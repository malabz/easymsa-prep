#include "common.hpp"

namespace easymsa_prep {

PrepError::PrepError(std::string code, std::string stage, std::string message)
    : std::runtime_error(message), code_(std::move(code)), stage_(std::move(stage)) {}

const std::string& PrepError::code() const noexcept { return code_; }
const std::string& PrepError::stage() const noexcept { return stage_; }

std::string version() { return EASYMSA_PREP_VERSION; }

int default_threads() { return std::min(omp_get_max_threads(), 8); }

void configure_threads(int requested_threads) {
  const int threads = requested_threads > 0 ? requested_threads : default_threads();
  omp_set_num_threads(std::max(1, threads));
}

bool ends_with(const std::string& text, const std::string& suffix) {
  return text.size() >= suffix.size() &&
         text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string lower(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

std::string path_kind(const fs::path& path) {
  const std::string name = lower(path.filename().string());
  if (ends_with(name, ".tar.gz")) return "tar_gz";
  if (ends_with(name, ".tgz")) return "tar_gz";
  if (ends_with(name, ".tar.xz")) return "tar_xz";
  if (ends_with(name, ".zip")) return "zip";
  if (ends_with(name, ".gz")) return "fasta_gz";
  if (ends_with(name, ".xz")) return "fasta_xz";
  if (ends_with(name, ".fasta") || ends_with(name, ".fa") ||
      ends_with(name, ".fna") || ends_with(name, ".ffn") ||
      ends_with(name, ".txt")) {
    return "fasta";
  }
  return "ignored";
}

bool is_fasta_kind(const std::string& kind) {
  return kind == "fasta" || kind == "fasta_gz" || kind == "fasta_xz";
}

bool is_archive_kind(const std::string& kind) {
  return kind == "zip" || kind == "tar_gz" || kind == "tar_xz";
}

bool is_nested_archive(const fs::path& path) {
  const std::string kind = path_kind(path);
  return is_archive_kind(kind);
}

bool is_hidden_or_system(const fs::path& path) {
  for (const auto& part : path) {
    const std::string text = part.string();
    if (text.empty()) continue;
    if (text == "__MACOSX" || text == ".DS_Store") return true;
    if (text[0] == '.') return true;
  }
  return false;
}

bool is_safe_archive_member(const fs::path& path) {
  const std::string text = path.generic_string();
  return !text.empty() && !path.is_absolute() && text.find("..") == std::string::npos &&
         !is_hidden_or_system(path) && !is_nested_archive(path);
}

void ensure_parent(const fs::path& path) {
  if (!path.parent_path().empty()) fs::create_directories(path.parent_path());
}

void ensure_output_dir(const fs::path& path) {
  if (path.empty()) throw PrepError("INVALID_CLI", "cli", "output directory is empty");
  fs::create_directories(path);
}

void check_output_overwrite(const fs::path& output_dir, const OutputConfig& output, bool force) {
  if (force) return;
  const std::array<fs::path, 3> generated = {
      output_dir / output.clean_fasta_name,
      output_dir / "result.json",
      output_dir / "summary.txt",
  };
  for (const auto& path : generated) {
    if (fs::exists(path)) {
      throw PrepError("IO_ERROR", "output", "output file exists; use --force: " + path.string());
    }
  }
}

std::ofstream open_output(const fs::path& path, std::ios::openmode mode) {
  ensure_parent(path);
  std::ofstream out(path, mode);
  if (!out) throw PrepError("IO_ERROR", "output", "failed to open output file: " + path.string());
  return out;
}

std::string join(const std::vector<std::string>& values, const std::string& sep) {
  std::ostringstream out;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) out << sep;
    out << values[i];
  }
  return out.str();
}

std::string wrap_fasta(const std::string& sequence, std::size_t width) {
  std::ostringstream out;
  width = std::max<std::size_t>(1, width);
  for (std::size_t i = 0; i < sequence.size(); i += width) {
    out << sequence.substr(i, width) << '\n';
  }
  return out.str();
}

std::uint64_t xxh64(const std::string& text) {
  return XXH64(text.data(), text.size(), 0);
}

bool is_legal_base(char c, SequenceType type) {
  c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  if (type == SequenceType::protein) {
    switch (c) {
      case 'A': case 'C': case 'D': case 'E': case 'F':
      case 'G': case 'H': case 'I': case 'K': case 'L':
      case 'M': case 'N': case 'P': case 'Q': case 'R':
      case 'S': case 'T': case 'V': case 'W': case 'Y':
      case 'B': case 'Z': case 'J': case 'X': case '*':
        return true;
      default:
        return false;
    }
  }
  switch (c) {
    case 'A':
    case 'C':
    case 'G':
    case 'N':
    case 'R':
    case 'Y':
    case 'S':
    case 'W':
    case 'K':
    case 'M':
    case 'B':
    case 'D':
    case 'H':
    case 'V':
      return true;
    case 'T':
      return type == SequenceType::dna;
    case 'U':
      return type == SequenceType::rna || type == SequenceType::dna;
    default:
      return false;
  }
}

char replacement_char(SequenceType type) {
  return type == SequenceType::protein ? 'X' : 'N';
}

std::size_t count_n_like(const std::string& sequence, SequenceType type) {
  const char unknown = type == SequenceType::protein ? 'X' : 'N';
  return static_cast<std::size_t>(std::count(sequence.begin(), sequence.end(), unknown));
}

double median_length(const std::vector<RecordResult>& results) {
  if (results.empty()) return 0.0;
  std::vector<std::size_t> lengths;
  lengths.reserve(results.size());
  for (const auto& result : results) lengths.push_back(result.original_length);
  std::sort(lengths.begin(), lengths.end());
  const auto mid = lengths.size() / 2;
  if (lengths.size() % 2 == 0) {
    return (static_cast<double>(lengths[mid - 1]) + static_cast<double>(lengths[mid])) / 2.0;
  }
  return static_cast<double>(lengths[mid]);
}

std::string infer_input_type(const fs::path& input) {
  if (fs::is_directory(input)) return "directory";
  const auto kind = path_kind(input);
  if (is_archive_kind(kind)) return "archive";
  if (is_fasta_kind(kind)) return "fasta";
  return "unknown";
}

std::string format_internal_id(std::size_t one_based_index) {
  std::ostringstream out;
  out << "seq_" << std::setw(6) << std::setfill('0') << one_based_index;
  return out.str();
}

TempDir::TempDir() {
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  path_ = fs::temp_directory_path() /
          ("easymsa-prep-" + std::to_string(now) + "-" + std::to_string(std::rand()));
  fs::create_directories(path_);
}

TempDir::~TempDir() {
  std::error_code ec;
  fs::remove_all(path_, ec);
}

const fs::path& TempDir::path() const { return path_; }

json thresholds_to_json(const Thresholds& thresholds) {
  return {
      {"max_n_ratio", thresholds.max_n_ratio},
      {"max_illegal_char_ratio", thresholds.max_illegal_char_ratio},
      {"min_length_ratio", thresholds.min_length_ratio},
      {"max_length_ratio", thresholds.max_length_ratio},
      {"min_similarity_to_reference", thresholds.min_similarity_to_reference},
      {"max_low_complexity_ratio", thresholds.max_low_complexity_ratio},
  };
}

json rules_to_json(const Rules& rules) {
  return {
      {"duplicate_ids", rules.duplicate_ids},
      {"lowercase", rules.lowercase},
      {"u_to_t", rules.u_to_t},
      {"gaps", rules.gaps},
      {"illegal_characters", rules.illegal_characters},
      {"high_n_ratio", rules.high_n_ratio},
      {"high_illegal_char_ratio", rules.high_illegal_char_ratio},
      {"too_short", rules.too_short},
      {"too_long", rules.too_long},
      {"low_similarity_outlier", rules.low_similarity_outlier},
      {"low_complexity", rules.low_complexity},
      {"identical_sequences", rules.identical_sequences},
  };
}

json limits_to_json(const Limits& limits) {
  return {
      {"max_compressed_size", limits.max_compressed_size},
      {"max_uncompressed_size", limits.max_uncompressed_size},
      {"max_file_count", limits.max_file_count},
      {"max_fasta_file_count", limits.max_fasta_file_count},
      {"max_raw_sequence_count", limits.max_raw_sequence_count},
      {"max_unique_sequence_count", limits.max_unique_sequence_count},
      {"max_total_bases", limits.max_total_bases},
      {"max_sequence_length", limits.max_sequence_length},
  };
}

json output_to_json(const OutputConfig& output) {
  return {{"clean_fasta_name", output.clean_fasta_name}, {"wrap_width", output.wrap_width}};
}

}  // namespace easymsa_prep
