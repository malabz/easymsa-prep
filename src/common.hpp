#pragma once

#include "easymsa_prep.hpp"
#include "miniz.h"
#include "xxhash.h"

#include <lzma.h>
#include <omp.h>
#include <zlib.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace easymsa_prep {

constexpr std::size_t kMaxIdLength = 200;
constexpr std::size_t kReferenceCandidateLimit = 1000;

class TempDir {
 public:
  TempDir();
  ~TempDir();
  const fs::path& path() const;

 private:
  fs::path path_;
};

bool ends_with(const std::string& text, const std::string& suffix);
std::string lower(std::string text);
std::string path_kind(const fs::path& path);
bool is_fasta_kind(const std::string& kind);
bool is_archive_kind(const std::string& kind);
bool is_nested_archive(const fs::path& path);
bool is_hidden_or_system(const fs::path& path);
bool is_safe_archive_member(const fs::path& path);
void ensure_parent(const fs::path& path);
void ensure_output_dir(const fs::path& path);
void check_output_overwrite(const fs::path& output_dir, const OutputConfig& output, bool force);
std::ofstream open_output(const fs::path& path, std::ios::openmode mode = std::ios::out);
std::string join(const std::vector<std::string>& values, const std::string& sep);
std::string wrap_fasta(const std::string& sequence, std::size_t width);
std::uint64_t xxh64(const std::string& text);

bool is_legal_base(char c, SequenceType type);
char replacement_char(SequenceType type);
std::size_t count_n_like(const std::string& sequence, SequenceType type);
double median_length(const std::vector<RecordResult>& results);
std::string infer_input_type(const fs::path& input);

std::vector<fs::path> discover_fasta_paths(const fs::path& input,
                                           InputBundle& bundle,
                                           const Config& config,
                                           std::vector<std::unique_ptr<TempDir>>& temp_dirs);
std::vector<Record> read_fasta_records(const std::vector<fs::path>& files,
                                       InputBundle& bundle,
                                       const Config& config);
void extract_archive(const fs::path& archive_path,
                     InputBundle& bundle,
                     const Config& config,
                     std::vector<std::unique_ptr<TempDir>>& temp_dirs);

json thresholds_to_json(const Thresholds& thresholds);
json rules_to_json(const Rules& rules);
json limits_to_json(const Limits& limits);
json output_to_json(const OutputConfig& output);

}  // namespace easymsa_prep
