#include "common.hpp"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

namespace easymsa_prep {
namespace {

std::vector<unsigned char> read_file_bytes(const fs::path& path, std::uintmax_t limit) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw PrepError("IO_ERROR", "input", "failed to open file: " + path.string());
  std::vector<unsigned char> bytes;
  std::array<char, 64 * 1024> buffer{};
  std::uintmax_t total = 0;
  while (in) {
    in.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const auto got = in.gcount();
    if (got <= 0) break;
    total += static_cast<std::uintmax_t>(got);
    if (total > limit) throw PrepError("LIMIT_EXCEEDED", "input", "file exceeds configured limit: " + path.string());
    bytes.insert(bytes.end(), buffer.begin(), buffer.begin() + got);
  }
  return bytes;
}

std::vector<unsigned char> decompress_xz_file(const fs::path& path, const Limits& limits) {
  if (fs::file_size(path) > limits.max_compressed_size) {
    throw PrepError("LIMIT_EXCEEDED", "archive", "xz input exceeds max_compressed_size: " + path.string());
  }
  std::ifstream in(path, std::ios::binary);
  if (!in) throw PrepError("IO_ERROR", "input", "failed to open xz file: " + path.string());

  lzma_stream stream = LZMA_STREAM_INIT;
  lzma_ret ret = lzma_stream_decoder(&stream, UINT64_MAX, LZMA_CONCATENATED);
  if (ret != LZMA_OK) throw PrepError("INVALID_ARCHIVE", "archive", "failed to initialize xz decoder");

  std::array<unsigned char, 64 * 1024> inbuf{};
  std::array<unsigned char, 64 * 1024> outbuf{};
  std::vector<unsigned char> output;
  lzma_action action = LZMA_RUN;

  stream.next_in = nullptr;
  stream.avail_in = 0;
  do {
    if (stream.avail_in == 0 && !in.eof()) {
      in.read(reinterpret_cast<char*>(inbuf.data()), static_cast<std::streamsize>(inbuf.size()));
      stream.next_in = inbuf.data();
      stream.avail_in = static_cast<std::size_t>(in.gcount());
      if (in.eof()) action = LZMA_FINISH;
    }

    stream.next_out = outbuf.data();
    stream.avail_out = outbuf.size();
    ret = lzma_code(&stream, action);
    const std::size_t produced = outbuf.size() - stream.avail_out;
    if (produced) {
      output.insert(output.end(), outbuf.begin(), outbuf.begin() + static_cast<std::ptrdiff_t>(produced));
      if (output.size() > limits.max_uncompressed_size) {
        lzma_end(&stream);
        throw PrepError("LIMIT_EXCEEDED", "archive", "xz uncompressed size exceeds limit: " + path.string());
      }
    }
    if (ret != LZMA_OK && ret != LZMA_STREAM_END) {
      lzma_end(&stream);
      throw PrepError("INVALID_ARCHIVE", "archive", "failed to decompress xz file: " + path.string());
    }
  } while (ret != LZMA_STREAM_END);

  lzma_end(&stream);
  return output;
}

std::vector<unsigned char> decompress_gzip_file(const fs::path& path, const Limits& limits) {
  if (fs::file_size(path) > limits.max_compressed_size) {
    throw PrepError("LIMIT_EXCEEDED", "archive", "gzip input exceeds max_compressed_size: " + path.string());
  }
  gzFile file = gzopen(path.string().c_str(), "rb");
  if (!file) throw PrepError("INVALID_ARCHIVE", "archive", "failed to open gzip file: " + path.string());
  std::vector<unsigned char> output;
  std::array<unsigned char, 64 * 1024> buffer{};
  while (true) {
    const int got = gzread(file, buffer.data(), static_cast<unsigned int>(buffer.size()));
    if (got < 0) {
      int err = 0;
      const char* msg = gzerror(file, &err);
      gzclose(file);
      throw PrepError("INVALID_ARCHIVE", "archive", msg ? msg : "gzip read error");
    }
    if (got == 0) break;
    output.insert(output.end(), buffer.begin(), buffer.begin() + got);
    if (output.size() > limits.max_uncompressed_size) {
      gzclose(file);
      throw PrepError("LIMIT_EXCEEDED", "archive", "gzip uncompressed size exceeds limit: " + path.string());
    }
  }
  gzclose(file);
  return output;
}

std::uint64_t parse_tar_octal(const unsigned char* data, std::size_t size) {
  std::uint64_t value = 0;
  std::size_t i = 0;
  while (i < size && (data[i] == ' ' || data[i] == '\0')) ++i;
  for (; i < size && data[i] >= '0' && data[i] <= '7'; ++i) {
    value = (value << 3) + static_cast<std::uint64_t>(data[i] - '0');
  }
  return value;
}

bool tar_block_is_zero(const unsigned char* block) {
  for (std::size_t i = 0; i < 512; ++i) {
    if (block[i] != 0) return false;
  }
  return true;
}

std::string tar_string(const unsigned char* data, std::size_t size) {
  std::size_t len = 0;
  while (len < size && data[len] != '\0') ++len;
  return std::string(reinterpret_cast<const char*>(data), len);
}

std::string normalize_tar_member_name(std::string name) {
  while (name == "." || name.rfind("./", 0) == 0) {
    if (name == ".") return "";
    name.erase(0, 2);
  }
  return name;
}

fs::path materialize_xz_fasta(const fs::path& path,
                              const Config& config,
                              std::vector<std::unique_ptr<TempDir>>& temp_dirs);

fs::path materialize_xz_bytes(const std::vector<unsigned char>& bytes,
                              const fs::path& source_name,
                              const Config& config,
                              std::vector<std::unique_ptr<TempDir>>& temp_dirs);

void extract_tar_bytes(const std::vector<unsigned char>& bytes,
                       const fs::path& archive_path,
                       InputBundle& bundle,
                       const Config& config,
                       std::vector<std::unique_ptr<TempDir>>& temp_dirs) {
  auto temp = std::make_unique<TempDir>();
  std::uintmax_t total_uncompressed = 0;
  std::size_t file_count = 0;

  std::size_t offset = 0;
  while (offset + 512 <= bytes.size()) {
    const unsigned char* header = bytes.data() + offset;
    if (tar_block_is_zero(header)) break;

    std::string name = tar_string(header, 100);
    const std::string prefix = tar_string(header + 345, 155);
    if (!prefix.empty()) name = prefix + "/" + name;
    name = normalize_tar_member_name(name);
    const char typeflag = static_cast<char>(header[156] ? header[156] : '0');
    const std::uint64_t size = parse_tar_octal(header + 124, 12);
    offset += 512;
    if (offset + size > bytes.size()) {
      throw PrepError("INVALID_ARCHIVE", "archive", "truncated tar archive: " + archive_path.string());
    }

    InputFileInfo info;
    info.path = archive_path.string() + "!" + name;
    info.kind = "tar_member";
    info.uncompressed_size = size;

    const fs::path member(name);
    if (typeflag == '5') {
      info.status = "ignored";
      info.note = "directory";
      bundle.files.push_back(info);
    } else if (typeflag == '0' || typeflag == '\0') {
      if (++file_count > config.limits.max_file_count) {
        throw PrepError("LIMIT_EXCEEDED", "archive", "tar file count exceeds limit: " + archive_path.string());
      }
      if (!is_safe_archive_member(member)) {
        info.status = "ignored";
        info.note = "unsafe_or_unsupported_member";
        bundle.files.push_back(info);
      } else if (!is_fasta_kind(path_kind(member))) {
        info.status = "ignored";
        info.note = "non_fasta_member";
        ++bundle.ignored_file_count;
        bundle.files.push_back(info);
      } else {
        total_uncompressed += size;
        if (total_uncompressed > config.limits.max_uncompressed_size) {
          throw PrepError("LIMIT_EXCEEDED", "archive", "tar uncompressed size exceeds limit: " + archive_path.string());
        }
        fs::path out = temp->path() / member;
        if (path_kind(member) == "fasta_xz") {
          std::vector<unsigned char> member_bytes(bytes.begin() + static_cast<std::ptrdiff_t>(offset),
                                                  bytes.begin() + static_cast<std::ptrdiff_t>(offset + size));
          out = materialize_xz_bytes(member_bytes, member, config, temp_dirs);
        } else {
          ensure_parent(out);
          std::ofstream file(out, std::ios::binary);
          if (!file) throw PrepError("IO_ERROR", "archive", "failed to extract tar member: " + name);
          file.write(reinterpret_cast<const char*>(bytes.data() + offset), static_cast<std::streamsize>(size));
        }
        info.status = "accepted";
        info.path = out;
        bundle.files.push_back(info);
      }
    }

    offset += static_cast<std::size_t>((size + 511) / 512) * 512;
  }

  bundle.uncompressed_size += total_uncompressed;
  temp_dirs.push_back(std::move(temp));
}

void extract_zip_archive(const fs::path& archive_path,
                         InputBundle& bundle,
                         const Config& config,
                         std::vector<std::unique_ptr<TempDir>>& temp_dirs) {
  if (fs::file_size(archive_path) > config.limits.max_compressed_size) {
    throw PrepError("LIMIT_EXCEEDED", "archive", "zip input exceeds max_compressed_size: " + archive_path.string());
  }

  mz_zip_archive archive{};
  if (!mz_zip_reader_init_file(&archive, archive_path.string().c_str(), 0)) {
    throw PrepError("INVALID_ARCHIVE", "archive", "failed to open zip: " + archive_path.string());
  }

  auto cleanup = [&archive]() { mz_zip_reader_end(&archive); };
  auto temp = std::make_unique<TempDir>();
  std::uintmax_t total_uncompressed = 0;
  std::size_t accepted_files = 0;

  const mz_uint count = mz_zip_reader_get_num_files(&archive);
  if (count > config.limits.max_file_count) {
    cleanup();
    throw PrepError("LIMIT_EXCEEDED", "archive", "zip contains too many files: " + archive_path.string());
  }

  std::vector<std::pair<std::string, mz_uint>> members;
  members.reserve(count);
  for (mz_uint i = 0; i < count; ++i) {
    mz_zip_archive_file_stat stat{};
    if (!mz_zip_reader_file_stat(&archive, i, &stat)) continue;
    members.emplace_back(stat.m_filename, i);
  }
  std::sort(members.begin(), members.end());

  for (const auto& [name_text, index] : members) {
    mz_zip_archive_file_stat stat{};
    if (!mz_zip_reader_file_stat(&archive, index, &stat)) continue;
    fs::path name(name_text);
    InputFileInfo info;
    info.path = archive_path.string() + "!" + name.generic_string();
    info.kind = "zip_member";
    info.compressed_size = stat.m_comp_size;
    info.uncompressed_size = stat.m_uncomp_size;

    if (mz_zip_reader_is_file_a_directory(&archive, index)) {
      info.status = "ignored";
      info.note = "directory";
      bundle.files.push_back(info);
      continue;
    }
    if (!is_safe_archive_member(name)) {
      info.status = "ignored";
      info.note = "unsafe_or_unsupported_member";
      bundle.files.push_back(info);
      continue;
    }
    if (!is_fasta_kind(path_kind(name))) {
      info.status = "ignored";
      info.note = "non_fasta_member";
      ++bundle.ignored_file_count;
      bundle.files.push_back(info);
      continue;
    }
    total_uncompressed += stat.m_uncomp_size;
    if (total_uncompressed > config.limits.max_uncompressed_size) {
      cleanup();
      throw PrepError("LIMIT_EXCEEDED", "archive", "zip uncompressed size exceeds limit: " + archive_path.string());
    }
    if (++accepted_files > config.limits.max_fasta_file_count) {
      cleanup();
      throw PrepError("LIMIT_EXCEEDED", "archive", "too many FASTA files in zip: " + archive_path.string());
    }
    fs::path out = temp->path() / name;
    ensure_parent(out);
    if (!mz_zip_reader_extract_to_file(&archive, index, out.string().c_str(), 0)) {
      cleanup();
      throw PrepError("INVALID_ARCHIVE", "archive", "failed to extract zip member: " + name.generic_string());
    }
    if (path_kind(name) == "fasta_xz") {
      out = materialize_xz_fasta(out, config, temp_dirs);
    }
    info.status = "accepted";
    info.path = out;
    bundle.files.push_back(info);
  }

  bundle.uncompressed_size += total_uncompressed;
  cleanup();
  temp_dirs.push_back(std::move(temp));
}

fs::path materialize_xz_fasta(const fs::path& path,
                              const Config& config,
                              std::vector<std::unique_ptr<TempDir>>& temp_dirs) {
  auto temp = std::make_unique<TempDir>();
  auto bytes = decompress_xz_file(path, config.limits);
  const fs::path out = temp->path() / path.stem();
  ensure_parent(out);
  std::ofstream file(out, std::ios::binary);
  if (!file) throw PrepError("IO_ERROR", "input", "failed to materialize xz FASTA: " + path.string());
  file.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
  temp_dirs.push_back(std::move(temp));
  return out;
}

fs::path materialize_xz_bytes(const std::vector<unsigned char>& bytes,
                              const fs::path& source_name,
                              const Config& config,
                              std::vector<std::unique_ptr<TempDir>>& temp_dirs) {
  auto temp = std::make_unique<TempDir>();
  const fs::path compressed = temp->path() / source_name;
  ensure_parent(compressed);
  {
    std::ofstream file(compressed, std::ios::binary);
    if (!file) throw PrepError("IO_ERROR", "archive", "failed to write xz archive member: " + source_name.string());
    file.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
  }
  auto decompressed = decompress_xz_file(compressed, config.limits);
  const fs::path out = temp->path() / source_name.stem();
  {
    std::ofstream file(out, std::ios::binary);
    if (!file) throw PrepError("IO_ERROR", "archive", "failed to materialize xz archive member: " + source_name.string());
    file.write(reinterpret_cast<const char*>(decompressed.data()), static_cast<std::streamsize>(decompressed.size()));
  }
  temp_dirs.push_back(std::move(temp));
  return out;
}

}  // namespace

void extract_archive(const fs::path& archive_path,
                     InputBundle& bundle,
                     const Config& config,
                     std::vector<std::unique_ptr<TempDir>>& temp_dirs) {
  const auto kind = path_kind(archive_path);
  bundle.compressed_size += fs::exists(archive_path) ? fs::file_size(archive_path) : 0;
  if (kind == "zip") {
    extract_zip_archive(archive_path, bundle, config, temp_dirs);
  } else if (kind == "tar_gz") {
    const auto bytes = decompress_gzip_file(archive_path, config.limits);
    extract_tar_bytes(bytes, archive_path, bundle, config, temp_dirs);
  } else if (kind == "tar_xz") {
    const auto bytes = decompress_xz_file(archive_path, config.limits);
    extract_tar_bytes(bytes, archive_path, bundle, config, temp_dirs);
  } else {
    throw PrepError("INVALID_ARCHIVE", "archive", "unsupported archive: " + archive_path.string());
  }
}

std::vector<fs::path> discover_fasta_paths(const fs::path& input,
                                           InputBundle& bundle,
                                           const Config& config,
                                           std::vector<std::unique_ptr<TempDir>>& temp_dirs) {
  if (!fs::exists(input)) throw PrepError("IO_ERROR", "input", "input does not exist: " + input.string());
  bundle.input_path = input;
  bundle.input_type = infer_input_type(input);

  std::vector<fs::path> paths;
  auto add_file = [&](const fs::path& path) {
    const auto kind = path_kind(path);
    InputFileInfo info;
    info.path = path;
    info.kind = kind;
    if (fs::exists(path) && fs::is_regular_file(path)) info.compressed_size = fs::file_size(path);

    if (is_archive_kind(kind)) {
      if (is_nested_archive(path) && path != input) {
        info.status = "ignored";
        info.note = "nested_archive";
        bundle.files.push_back(info);
        return;
      }
      info.status = "accepted";
      bundle.files.push_back(info);
      const auto before = bundle.files.size();
      extract_archive(path, bundle, config, temp_dirs);
      for (std::size_t i = before; i < bundle.files.size(); ++i) {
        if (bundle.files[i].status == "accepted" && is_fasta_kind(path_kind(bundle.files[i].path))) {
          paths.push_back(bundle.files[i].path);
        }
      }
      return;
    }

    if (!is_fasta_kind(kind)) {
      info.status = "ignored";
      info.note = "unsupported_extension";
      ++bundle.ignored_file_count;
      bundle.files.push_back(info);
      return;
    }

    if (kind == "fasta_xz") {
      info.status = "accepted";
      bundle.compressed_size += info.compressed_size;
      bundle.files.push_back(info);
      paths.push_back(materialize_xz_fasta(path, config, temp_dirs));
    } else {
      info.status = "accepted";
      bundle.compressed_size += info.compressed_size;
      bundle.files.push_back(info);
      paths.push_back(path);
    }
  };

  if (fs::is_directory(input)) {
    std::vector<fs::path> discovered;
    for (const auto& entry : fs::recursive_directory_iterator(input)) {
      if (entry.is_regular_file() && !is_hidden_or_system(entry.path().lexically_relative(input))) {
        discovered.push_back(entry.path());
      }
    }
    std::sort(discovered.begin(), discovered.end());
    for (const auto& path : discovered) add_file(path);
  } else {
    add_file(input);
  }

  std::sort(paths.begin(), paths.end());
  paths.erase(std::unique(paths.begin(), paths.end()), paths.end());
  if (paths.empty()) throw PrepError("IO_ERROR", "input", "no FASTA inputs found");
  if (paths.size() > config.limits.max_fasta_file_count) {
    throw PrepError("LIMIT_EXCEEDED", "input", "too many FASTA files");
  }
  return paths;
}

// kseq separates names from comments and discards their delimiter. Retain
// complete FASTA headers separately so tabs/spaces survive result publication.
std::vector<std::string> read_original_headers(const fs::path& path) {
  gzFile input = gzopen(path.string().c_str(), "rb");
  if (!input) throw PrepError("IO_ERROR", "input", "failed to read FASTA headers");
  std::vector<std::string> headers;
  std::array<char, 8192> buffer{};
  std::string header;
  bool start = true, capture = false, first = true;
  while (gzgets(input, buffer.data(), static_cast<int>(buffer.size()))) {
    std::string part(buffer.data());
    if (first && part.starts_with("\xef\xbb\xbf")) part.erase(0, 3);
    first = false;
    if (start) { capture = !part.empty() && part[0] == '>'; header.clear(); }
    if (capture) header += part;
    start = !part.empty() && part.back() == '\n';
    if (start && capture) {
      while (!header.empty() && (header.back() == '\n' || header.back() == '\r')) header.pop_back();
      headers.push_back(header.substr(1));
    }
  }
  if (!start && capture) {
    if (!header.empty() && header.back() == '\r') header.pop_back();
    headers.push_back(header.substr(1));
  }
  int error = Z_OK; gzerror(input, &error); gzclose(input);
  if (error != Z_OK && error != Z_STREAM_END) throw PrepError("IO_ERROR", "input", "failed to read FASTA headers");
  return headers;
}

std::vector<Record> read_fasta_records(const std::vector<fs::path>& files,
                                       InputBundle& bundle,
                                       const Config& config) {
  std::vector<Record> all;
  std::size_t next_index = 0;
  std::uint64_t total_bases = 0;

  for (const auto& path : files) {
    const auto headers = read_original_headers(path);
    gzFile fp = gzopen(path.string().c_str(), "rb");
    if (!fp) throw PrepError("INVALID_FASTA", "fasta_validation", "failed to open FASTA: " + path.string());
    kseq_t* seq = kseq_init(fp);
    std::vector<Record> records;
    std::uint64_t bases = 0;
    int status = 0;
    while ((status = kseq_read(seq)) >= 0) {
      Record record;
      record.record_index = next_index++;
      record.internal_id = format_internal_id(record.record_index + 1);
      record.original_id = seq->name.s ? seq->name.s : "";
      if (records.size() >= headers.size()) {
        kseq_destroy(seq); gzclose(fp);
        throw PrepError("INVALID_FASTA", "fasta_validation", "FASTA header count mismatch");
      }
      record.original_header = headers[records.size()];
      const auto separator = record.original_header.find_first_of(" \t");
      const auto description_start = separator == std::string::npos ? std::string::npos
          : record.original_header.find_first_not_of(" \t", separator);
      record.description = description_start == std::string::npos ? "" : record.original_header.substr(description_start);
      record.source_file = path;
      record.sequence = seq->seq.s ? seq->seq.s : "";
      if (record.original_id.empty()) {
        kseq_destroy(seq);
        gzclose(fp);
        throw PrepError("INVALID_FASTA", "fasta_validation", "empty FASTA header in " + path.string());
      }
      if (record.sequence.empty()) {
        kseq_destroy(seq);
        gzclose(fp);
        throw PrepError("INVALID_FASTA", "fasta_validation", "record has header but no sequence in " + path.string());
      }
      for (unsigned char c : record.sequence) {
        if (std::iscntrl(c) && !std::isspace(c)) {
          kseq_destroy(seq);
          gzclose(fp);
          throw PrepError("INVALID_FASTA", "fasta_validation", "control character in sequence: " + path.string());
        }
      }
      bases += record.sequence.size();
      total_bases += record.sequence.size();
      if (record.sequence.size() > config.limits.max_sequence_length) {
        kseq_destroy(seq);
        gzclose(fp);
        throw PrepError("LIMIT_EXCEEDED", "fasta_validation", "sequence exceeds max_sequence_length in " + path.string());
      }
      if (total_bases > config.limits.max_total_bases) {
        kseq_destroy(seq);
        gzclose(fp);
        throw PrepError("LIMIT_EXCEEDED", "fasta_validation", "input exceeds max_total_bases");
      }
      records.push_back(std::move(record));
    }
    if (status < -1) {
      kseq_destroy(seq);
      gzclose(fp);
      throw PrepError("INVALID_FASTA", "fasta_validation", "kseq parse error in " + path.string());
    }
    kseq_destroy(seq);
    gzclose(fp);
    if (records.empty()) throw PrepError("INVALID_FASTA", "fasta_validation", "empty FASTA file: " + path.string());

    for (auto& info : bundle.files) {
      const bool exact_match = info.path == path;
      const bool materialized_xz_match =
          path_kind(info.path) == "fasta_xz" &&
          info.path.filename().stem() == path.filename();
      if (exact_match || materialized_xz_match) {
        info.sequence_count = records.size();
        info.bases = bases;
        break;
      }
    }
    all.insert(all.end(), std::make_move_iterator(records.begin()), std::make_move_iterator(records.end()));
    if (all.size() > config.limits.max_raw_sequence_count) {
      throw PrepError("LIMIT_EXCEEDED", "fasta_validation", "input exceeds max_raw_sequence_count");
    }
  }

  return all;
}

InputBundle read_input_bundle(const fs::path& input, const Config& config) {
  InputBundle bundle;
  std::vector<std::unique_ptr<TempDir>> temp_dirs;
  const auto files = discover_fasta_paths(input, bundle, config, temp_dirs);
  bundle.records = read_fasta_records(files, bundle, config);
  return bundle;
}

}  // namespace easymsa_prep
