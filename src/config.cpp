#include "common.hpp"

namespace easymsa_prep {
namespace {

void require_object(const json& root, const std::string& name) {
  if (!root.is_object()) throw PrepError("INVALID_PARAMS_JSON", "params", name + " must be an object");
}

void check_keys(const json& object, const std::set<std::string>& allowed, const std::string& where) {
  require_object(object, where);
  for (const auto& [key, value] : object.items()) {
    (void)value;
    if (allowed.find(key) == allowed.end()) {
      throw PrepError("INVALID_PARAMS_JSON", "params", "unknown field in " + where + ": " + key);
    }
  }
}

template <class T>
void assign_number(const json& object, const std::string& key, T& target) {
  if (!object.contains(key)) return;
  if (!object.at(key).is_number()) {
    throw PrepError("INVALID_PARAMS_JSON", "params", key + " must be numeric");
  }
  target = object.at(key).get<T>();
}

void assign_string(const json& object, const std::string& key, std::string& target) {
  if (!object.contains(key)) return;
  if (!object.at(key).is_string()) {
    throw PrepError("INVALID_PARAMS_JSON", "params", key + " must be a string");
  }
  target = object.at(key).get<std::string>();
}

void assign_bool(const json& object, const std::string& key, bool& target) {
  if (!object.contains(key)) return;
  if (!object.at(key).is_boolean()) {
    throw PrepError("INVALID_PARAMS_JSON", "params", key + " must be a boolean");
  }
  target = object.at(key).get<bool>();
}

void validate_choice(const std::string& value,
                     const std::set<std::string>& choices,
                     const std::string& field) {
  if (choices.find(value) == choices.end()) {
    throw PrepError("INVALID_PARAMS_JSON", "params", "invalid " + field + ": " + value);
  }
}

void merge_params(Config& config, const json& root) {
  check_keys(root,
             {"mode", "strictness", "sequence_type", "thresholds", "rules", "limits", "output"},
             "params");

  if (root.contains("mode") || root.contains("strictness")) {
    std::string mode = config.mode;
    std::string strictness = config.strictness;
    if (root.contains("mode")) {
      if (!root.at("mode").is_string()) {
        throw PrepError("INVALID_PARAMS_JSON", "params", "mode must be a string");
      }
      mode = root.at("mode").get<std::string>();
    }
    if (root.contains("strictness")) {
      if (!root.at("strictness").is_string()) {
        throw PrepError("INVALID_PARAMS_JSON", "params", "strictness must be a string");
      }
      strictness = root.at("strictness").get<std::string>();
    }
    config = default_config_for_mode(mode, strictness);
  }
  if (root.contains("sequence_type")) {
    if (!root.at("sequence_type").is_string()) {
      throw PrepError("INVALID_PARAMS_JSON", "params", "sequence_type must be a string");
    }
    config.sequence_type = parse_sequence_type(root.at("sequence_type").get<std::string>());
  }

  if (root.contains("thresholds")) {
    const auto& object = root.at("thresholds");
    check_keys(object,
               {"max_n_ratio", "max_illegal_char_ratio", "min_length_ratio", "max_length_ratio",
                "min_similarity_to_reference", "max_low_complexity_ratio"},
               "thresholds");
    assign_number(object, "max_n_ratio", config.thresholds.max_n_ratio);
    assign_number(object, "max_illegal_char_ratio", config.thresholds.max_illegal_char_ratio);
    assign_number(object, "min_length_ratio", config.thresholds.min_length_ratio);
    assign_number(object, "max_length_ratio", config.thresholds.max_length_ratio);
    assign_number(object, "min_similarity_to_reference", config.thresholds.min_similarity_to_reference);
    assign_number(object, "max_low_complexity_ratio", config.thresholds.max_low_complexity_ratio);
  }

  if (root.contains("rules")) {
    const auto& object = root.at("rules");
    check_keys(object,
               {"duplicate_ids", "lowercase", "u_to_t", "gaps", "illegal_characters",
                "high_n_ratio", "high_illegal_char_ratio", "too_short", "too_long",
                "low_similarity_outlier", "low_complexity", "identical_sequences"},
               "rules");
    assign_string(object, "duplicate_ids", config.rules.duplicate_ids);
    assign_string(object, "lowercase", config.rules.lowercase);
    assign_bool(object, "u_to_t", config.rules.u_to_t);
    assign_string(object, "gaps", config.rules.gaps);
    assign_string(object, "illegal_characters", config.rules.illegal_characters);
    assign_string(object, "high_n_ratio", config.rules.high_n_ratio);
    assign_string(object, "high_illegal_char_ratio", config.rules.high_illegal_char_ratio);
    assign_string(object, "too_short", config.rules.too_short);
    assign_string(object, "too_long", config.rules.too_long);
    assign_string(object, "low_similarity_outlier", config.rules.low_similarity_outlier);
    assign_string(object, "low_complexity", config.rules.low_complexity);
    assign_string(object, "identical_sequences", config.rules.identical_sequences);
  }

  if (root.contains("limits")) {
    const auto& object = root.at("limits");
    check_keys(object,
               {"max_compressed_size", "max_uncompressed_size", "max_file_count",
                "max_fasta_file_count", "max_raw_sequence_count", "max_unique_sequence_count",
                "max_total_bases", "max_sequence_length"},
               "limits");
    assign_number(object, "max_compressed_size", config.limits.max_compressed_size);
    assign_number(object, "max_uncompressed_size", config.limits.max_uncompressed_size);
    assign_number(object, "max_file_count", config.limits.max_file_count);
    assign_number(object, "max_fasta_file_count", config.limits.max_fasta_file_count);
    assign_number(object, "max_raw_sequence_count", config.limits.max_raw_sequence_count);
    assign_number(object, "max_unique_sequence_count", config.limits.max_unique_sequence_count);
    assign_number(object, "max_total_bases", config.limits.max_total_bases);
    assign_number(object, "max_sequence_length", config.limits.max_sequence_length);
  }

  if (root.contains("output")) {
    const auto& object = root.at("output");
    check_keys(object, {"clean_fasta_name", "wrap_width"}, "output");
    assign_string(object, "clean_fasta_name", config.output.clean_fasta_name);
    assign_number(object, "wrap_width", config.output.wrap_width);
  }
}

}  // namespace

std::string normalize_mode_name(std::string mode) {
  mode = lower(std::move(mode));
  if (mode == "audit" || mode == "filter") return mode;
  throw PrepError("INVALID_CLI", "cli", "invalid mode: " + mode);
}

std::string normalize_strictness_name(std::string strictness) {
  strictness = lower(std::move(strictness));
  if (strictness == "strict" || strictness == "normal" || strictness == "lenient") return strictness;
  throw PrepError("INVALID_CLI", "cli", "invalid strictness: " + strictness);
}

Config default_config_for_mode(std::string mode, std::string strictness) {
  Config config;
  config.mode = normalize_mode_name(std::move(mode));
  config.strictness = normalize_strictness_name(std::move(strictness));
  config.sequence_type = SequenceType::dna;
  if (config.strictness == "strict") {
    config.thresholds.max_n_ratio = 0.05;
    config.thresholds.max_illegal_char_ratio = 0.001;
    config.thresholds.min_length_ratio = 0.9;
    config.thresholds.max_length_ratio = 1.1;
    config.thresholds.min_similarity_to_reference = 0.95;
    config.thresholds.max_low_complexity_ratio = 0.6;
  } else if (config.strictness == "lenient") {
    config.thresholds.max_n_ratio = 0.5;
    config.thresholds.max_illegal_char_ratio = 0.05;
    config.thresholds.min_length_ratio = 0.3;
    config.thresholds.max_length_ratio = 2.0;
    config.thresholds.min_similarity_to_reference = 0.5;
    config.thresholds.max_low_complexity_ratio = 0.95;
  }

  if (config.mode == "filter") {
    config.rules.high_n_ratio = "remove";
    config.rules.high_illegal_char_ratio = "remove";
    config.rules.too_short = "remove";
    config.rules.too_long = "remove";
    config.rules.low_similarity_outlier = "remove";
    config.rules.low_complexity = "remove";
    config.rules.identical_sequences = "collapse";
  }
  return config;
}

void validate_config(const Config& config) {
  if (config.thresholds.max_n_ratio < 0.0 || config.thresholds.max_n_ratio > 1.0) {
    throw PrepError("INVALID_PARAMS_JSON", "params", "max_n_ratio must be 0..1");
  }
  if (config.thresholds.max_illegal_char_ratio < 0.0 ||
      config.thresholds.max_illegal_char_ratio > 1.0) {
    throw PrepError("INVALID_PARAMS_JSON", "params", "max_illegal_char_ratio must be 0..1");
  }
  if (config.thresholds.min_length_ratio < 0.0 ||
      config.thresholds.max_length_ratio < config.thresholds.min_length_ratio) {
    throw PrepError("INVALID_PARAMS_JSON", "params", "invalid length ratio thresholds");
  }
  if (config.thresholds.min_similarity_to_reference < 0.0 ||
      config.thresholds.min_similarity_to_reference > 1.0) {
    throw PrepError("INVALID_PARAMS_JSON", "params", "min_similarity_to_reference must be 0..1");
  }
  if (config.thresholds.max_low_complexity_ratio < 0.0 ||
      config.thresholds.max_low_complexity_ratio > 1.0) {
    throw PrepError("INVALID_PARAMS_JSON", "params", "max_low_complexity_ratio must be 0..1");
  }
  validate_choice(config.rules.duplicate_ids, {"rename_if_needed"}, "duplicate_ids");
  validate_choice(config.rules.lowercase, {"to_uppercase", "keep"}, "lowercase");
  validate_choice(config.rules.gaps, {"remove", "keep"}, "gaps");
  validate_choice(config.rules.illegal_characters,
                  {"replace_with_N", "remove_sequence", "warn_only"},
                  "illegal_characters");
  validate_choice(config.rules.high_n_ratio, {"remove", "warn_only"}, "high_n_ratio");
  validate_choice(config.rules.high_illegal_char_ratio,
                  {"remove", "warn_only"},
                  "high_illegal_char_ratio");
  validate_choice(config.rules.too_short, {"remove", "warn_only"}, "too_short");
  validate_choice(config.rules.too_long, {"remove", "warn_only"}, "too_long");
  validate_choice(config.rules.low_similarity_outlier,
                  {"remove", "warn_only"},
                  "low_similarity_outlier");
  validate_choice(config.rules.low_complexity, {"remove", "warn_only"}, "low_complexity");
  validate_choice(config.rules.identical_sequences, {"collapse", "keep"}, "identical_sequences");
  if (config.output.clean_fasta_name.empty() ||
      fs::path(config.output.clean_fasta_name).is_absolute() ||
      config.output.clean_fasta_name.find('/') != std::string::npos ||
      config.output.clean_fasta_name.find('\\') != std::string::npos) {
    throw PrepError("INVALID_PARAMS_JSON", "params", "clean_fasta_name must be a simple filename");
  }
  if (config.output.wrap_width == 0 || config.output.wrap_width > 1000) {
    throw PrepError("INVALID_PARAMS_JSON", "params", "wrap_width must be 1..1000");
  }
}

Config load_config(const CliOptions& options) {
  Config config = default_config_for_mode(options.mode, options.strictness);
  if (!options.params_path.empty()) {
    std::ifstream in(options.params_path);
    if (!in) {
      throw PrepError("INVALID_PARAMS_JSON", "params", "failed to open params JSON: " + options.params_path.string());
    }
    json root;
    try {
      in >> root;
    } catch (const std::exception& ex) {
      throw PrepError("INVALID_PARAMS_JSON", "params", std::string("failed to parse params JSON: ") + ex.what());
    }
    merge_params(config, root);
  }
  validate_config(config);
  return config;
}

}  // namespace easymsa_prep
