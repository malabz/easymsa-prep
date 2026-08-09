#pragma once

#include "json.hpp"

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace easymsa_prep {

namespace fs = std::filesystem;
using json = nlohmann::json;

enum class SequenceType { dna, rna, protein };

struct Limits {
  std::uintmax_t max_compressed_size = 1024ULL * 1024ULL * 1024ULL;
  std::uintmax_t max_uncompressed_size = 10ULL * 1024ULL * 1024ULL * 1024ULL;
  std::size_t max_file_count = 10000;
  std::size_t max_fasta_file_count = 1000;
  std::size_t max_raw_sequence_count = 10000000;
  std::size_t max_unique_sequence_count = 10000000;
  std::uint64_t max_total_bases = 50000000000ULL;
  std::size_t max_sequence_length = 10000000;
};

struct Thresholds {
  double max_n_ratio = 0.2;
  double max_illegal_char_ratio = 0.01;
  double min_length_ratio = 0.7;
  double max_length_ratio = 1.3;
  double min_similarity_to_reference = 0.8;
  double max_low_complexity_ratio = 0.8;
};

struct Rules {
  std::string duplicate_ids = "rename_if_needed";
  std::string lowercase = "to_uppercase";
  bool u_to_t = true;
  std::string gaps = "remove";
  std::string illegal_characters = "replace_with_N";
  std::string high_n_ratio = "warn_only";
  std::string high_illegal_char_ratio = "warn_only";
  std::string too_short = "warn_only";
  std::string too_long = "warn_only";
  std::string low_similarity_outlier = "warn_only";
  std::string low_complexity = "warn_only";
  std::string identical_sequences = "keep";
};

struct OutputConfig {
  std::string clean_fasta_name = "clean.fasta";
  std::size_t wrap_width = 80;
};

struct Config {
  std::string mode = "audit";
  std::string strictness = "normal";
  SequenceType sequence_type = SequenceType::dna;
  Thresholds thresholds;
  Rules rules;
  Limits limits;
  OutputConfig output;
};

struct CliOptions {
  fs::path input;
  fs::path output_dir;
  fs::path params_path;
  std::string mode = "audit";
  std::string strictness = "normal";
  int threads = 0;
  int batch_size = 10000;
  bool force = false;
  bool quiet = false;
  std::vector<std::string> command_line;
};

struct Record {
  std::size_t record_index = 0;
  std::string internal_id;
  std::string original_id;
  fs::path source_file;
  std::string sequence;
};

struct InputFileInfo {
  fs::path path;
  std::string kind;
  std::string status;
  std::uintmax_t compressed_size = 0;
  std::uintmax_t uncompressed_size = 0;
  std::size_t sequence_count = 0;
  std::uint64_t bases = 0;
  std::string note;
};

struct InputBundle {
  fs::path input_path;
  std::string input_type = "unknown";
  std::vector<InputFileInfo> files;
  std::vector<Record> records;
  std::size_t ignored_file_count = 0;
  std::uintmax_t compressed_size = 0;
  std::uintmax_t uncompressed_size = 0;
};

struct ReferenceInfo {
  std::string mode = "auto_first_1000";
  std::string selected_id;
  std::string selection_reason;
  bool available = false;
};

struct RecordResult {
  std::size_t record_index = 0;
  std::string internal_id;
  std::string original_id;
  std::string output_id;
  fs::path source_file;
  std::size_t original_length = 0;
  std::size_t clean_length = 0;
  double n_ratio = 0.0;
  double illegal_char_ratio = 0.0;
  std::size_t gap_count = 0;
  bool invalid_id = false;
  bool duplicate_id = false;
  bool renamed = false;
  bool low_complexity = false;
  double similarity_to_reference = 1.0;
  std::string clean_sequence;
  std::vector<std::string> warnings;
  std::vector<std::string> remove_reasons;
  bool removed = false;
  bool representative = true;
  std::string representative_id;
  std::uint64_t hash = 0;
};

struct SummaryCounts {
  std::size_t raw_sequence_count = 0;
  std::uint64_t raw_total_bases = 0;
  std::size_t clean_sequence_count = 0;
  std::uint64_t clean_total_bases = 0;
  std::size_t removed_sequence_count = 0;
  std::size_t collapsed_duplicate_count = 0;
};

struct DatasetFeatureValues {
  std::size_t n_sequences = 0;
  std::uint64_t total_bases = 0;
  double length_min = 0.0;
  double length_max = 0.0;
  double length_mean = 0.0;
  double length_median = 0.0;
  double length_std = 0.0;
  double length_iqr = 0.0;
  double gc_mean = 0.0;
  double gc_std = 0.0;
  double ambiguous_mean = 0.0;
  std::size_t kmer_k = 5;
  std::size_t kmer_sample_size = 0;
  double kmer_sample_fraction = 0.0;
  std::size_t kmer_distinct = 0;
  double kmer_diversity = 0.0;
  double similarity_mean = 1.0;
  double similarity_median = 1.0;
  double similarity_std = 0.0;
  std::size_t similarity_pair_count = 0;
  bool low_memory_mode = false;
};

struct DatasetFeatures {
  std::string schema = "easymsa-dataset-features/v1";
  std::string status = "not_computed";
  std::string profile = "uniform-jaccard-2000-v1";
  std::string source = "retained_representatives";
  std::string error_code;
  std::string error_message;
  DatasetFeatureValues values;
  std::chrono::milliseconds runtime{0};
};

struct RunResult {
  bool success = true;
  int exit_code = 0;
  std::string status = "success";
  std::string error_code;
  std::string error_message;
  std::string failed_stage;
  Config config;
  CliOptions options;
  InputBundle input;
  ReferenceInfo reference;
  DatasetFeatures dataset_features;
  SummaryCounts summary;
  std::vector<RecordResult> records;
  std::map<std::string, std::size_t> qc_counts;
  std::map<std::string, std::size_t> cleaning_counts;
  std::map<std::string, std::size_t> removal_counts;
  std::map<std::string, std::size_t> warning_counts;
  std::size_t duplicate_group_count = 0;
  std::chrono::milliseconds runtime{0};
};

class PrepError : public std::runtime_error {
 public:
  PrepError(std::string code, std::string stage, std::string message);
  const std::string& code() const noexcept;
  const std::string& stage() const noexcept;

 private:
  std::string code_;
  std::string stage_;
};

CliOptions parse_cli(int argc, char** argv);
Config load_config(const CliOptions& options);
DatasetFeatures compute_dataset_features(const std::vector<RecordResult>& records,
                                         SequenceType sequence_type);
RunResult run_pipeline(const CliOptions& options, const Config& config);
void write_success_outputs(const RunResult& result);
void write_error_outputs(const CliOptions& options,
                         const Config* config,
                         const PrepError& error,
                         int exit_code);

int default_threads();
void configure_threads(int requested_threads);
std::string version();

std::string normalize_mode_name(std::string mode);
std::string normalize_strictness_name(std::string strictness);
std::string sequence_type_name(SequenceType type);
SequenceType parse_sequence_type(const std::string& value);
Config default_config_for_mode(std::string mode, std::string strictness);
void validate_config(const Config& config);

InputBundle read_input_bundle(const fs::path& input, const Config& config);

std::string normalize_sequence(const std::string& sequence,
                               SequenceType type,
                               const Rules& rules,
                               std::size_t* gap_count = nullptr,
                               std::size_t* illegal_count = nullptr,
                               std::size_t* whitespace_count = nullptr);
bool valid_id(const std::string& id);
bool is_low_complexity(const std::string& sequence, double max_ratio);
double kmer_similarity(const std::string& lhs, const std::string& rhs, std::size_t k = 7);
std::string format_internal_id(std::size_t one_based_index);

}  // namespace easymsa_prep
