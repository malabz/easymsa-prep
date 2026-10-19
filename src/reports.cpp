#include "common.hpp"

namespace easymsa_prep {
namespace {

json input_file_to_json(const InputFileInfo& info) {
  return {
      {"path", info.path.string()},
      {"kind", info.kind},
      {"status", info.status},
      {"compressed_size", info.compressed_size},
      {"uncompressed_size", info.uncompressed_size},
      {"sequence_count", info.sequence_count},
      {"bases", info.bases},
      {"note", info.note},
  };
}

json map_to_json(const std::map<std::string, std::size_t>& values) {
  json out = json::object();
  for (const auto& [key, value] : values) out[key] = value;
  return out;
}

json command_line_to_json(const std::vector<std::string>& command_line) {
  json out = json::array();
  for (const auto& part : command_line) out.push_back(part);
  return out;
}

json dataset_features_to_json(const DatasetFeatures& features) {
  json values = json::object();
  if (features.status == "complete") {
    const auto& item = features.values;
    values = {
        {"n_sequences", item.n_sequences},
        {"total_bases", item.total_bases},
        {"length_min", item.length_min},
        {"length_max", item.length_max},
        {"length_mean", item.length_mean},
        {"length_median", item.length_median},
        {"length_std", item.length_std},
        {"length_iqr", item.length_iqr},
        {"gc_mean", item.gc_mean},
        {"gc_std", item.gc_std},
        {"ambiguous_mean", item.ambiguous_mean},
        {"kmer_k", item.kmer_k},
        {"kmer_sample_size", item.kmer_sample_size},
        {"kmer_sample_fraction", item.kmer_sample_fraction},
        {"kmer_distinct", item.kmer_distinct},
        {"kmer_diversity", item.kmer_diversity},
        {"similarity_mean", item.similarity_mean},
        {"similarity_median", item.similarity_median},
        {"similarity_std", item.similarity_std},
        {"similarity_pair_count", item.similarity_pair_count},
        {"low_memory_mode", item.low_memory_mode ? 1 : 0},
    };
  }

  json output = {
      {"schema", features.schema},
      {"status", features.status},
      {"profile", features.profile},
      {"source", features.source},
      {"parameters",
       {
           {"alphabet", "ACGT"},
           {"kmer_k", 5},
           {"similarity_sampling", "uniform_without_replacement"},
           {"similarity_max_pairs", 2000},
           {"seed", 42},
           {"standard_deviation", "sample"},
       }},
      {"values", values},
      {"timing_ms", features.runtime.count()},
  };
  if (!features.error_code.empty()) output["error_code"] = features.error_code;
  if (!features.error_message.empty()) output["error_message"] = features.error_message;
  return output;
}

json success_result_json(const RunResult& result) {
  json input_files = json::array();
  std::size_t fasta_file_count = 0;
  for (const auto& file : result.input.files) {
    input_files.push_back(input_file_to_json(file));
    if (file.status == "accepted" && is_fasta_kind(path_kind(file.path))) {
      ++fasta_file_count;
    }
  }

  return {
      {"status", "success"},
      {"success", true},
      {"exit_code", 0},
      {"program_version", version()},
      {"command_line", command_line_to_json(result.options.command_line)},
      {"input",
       {
           {"path", result.options.input.string()},
           {"type", result.input.input_type},
           {"fasta_file_count", fasta_file_count},
           {"ignored_file_count", result.input.ignored_file_count},
           {"compressed_size", result.input.compressed_size},
           {"uncompressed_size", result.input.uncompressed_size},
       }},
      {"params",
       {
           {"params_json", result.options.params_path.empty() ? nullptr : json(result.options.params_path.string())},
           {"batch_size", result.options.batch_size},
           {"threads", result.options.threads > 0 ? result.options.threads : default_threads()},
      }},
      {"mode", result.config.mode},
      {"strictness", result.config.strictness},
      {"sequence_type", sequence_type_name(result.config.sequence_type)},
      {"thresholds", thresholds_to_json(result.config.thresholds)},
      {"rules", rules_to_json(result.config.rules)},
      {"limits", limits_to_json(result.config.limits)},
      {"reference",
       {
           {"mode", result.reference.mode},
           {"available", result.reference.available},
           {"selected_id", result.reference.selected_id},
           {"selection_reason", result.reference.selection_reason},
       }},
      {"summary_counts",
       {
           {"raw_sequence_count", result.summary.raw_sequence_count},
           {"raw_total_bases", result.summary.raw_total_bases},
           {"clean_sequence_count", result.summary.clean_sequence_count},
           {"clean_total_bases", result.summary.clean_total_bases},
           {"removed_sequence_count", result.summary.removed_sequence_count},
           {"collapsed_duplicate_count", result.summary.collapsed_duplicate_count},
       }},
      {"input_files", input_files},
      {"dataset_features", dataset_features_to_json(result.dataset_features)},
      {"qc_counts", map_to_json(result.qc_counts)},
      {"cleaning_counts", map_to_json(result.cleaning_counts)},
      {"removal_counts", map_to_json(result.removal_counts)},
      {"dedup_summary",
       {
           {"enabled", result.config.rules.identical_sequences == "collapse"},
           {"duplicate_groups", result.duplicate_group_count},
           {"collapsed_count", result.summary.collapsed_duplicate_count},
       }},
      {"warnings", map_to_json(result.warning_counts)},
      {"outputs",
       {
           {"clean_fasta", result.config.output.clean_fasta_name},
           {"sequence_map", "sequence_map.jsonl"},
           {"sequence_map_tsv", "sequence_map.tsv"},
           {"result_json", "result.json"},
           {"summary_txt", "summary.txt"},
       }},
      {"runtime", {{"milliseconds", result.runtime.count()}}},
  };
}

void write_summary_txt(const RunResult& result, const fs::path& path) {
  auto out = open_output(path);
  out << "easymsa-prep " << version() << "\n";
  out << "Status: success\n";
  out << "Input: " << result.options.input << "\n";
  out << "Input type: " << result.input.input_type << "\n";
  out << "Mode: " << result.config.mode << "\n";
  out << "Strictness: " << result.config.strictness << "\n";
  out << "Sequence type: " << sequence_type_name(result.config.sequence_type) << "\n";
  out << "Params JSON: " << (result.options.params_path.empty() ? "not used" : result.options.params_path.string()) << "\n";
  out << "Threads: " << (result.options.threads > 0 ? result.options.threads : default_threads()) << "\n\n";

  out << "Thresholds:\n";
  out << "- max_n_ratio: " << result.config.thresholds.max_n_ratio << "\n";
  out << "- max_illegal_char_ratio: " << result.config.thresholds.max_illegal_char_ratio << "\n";
  out << "- min_length_ratio: " << result.config.thresholds.min_length_ratio << "\n";
  out << "- max_length_ratio: " << result.config.thresholds.max_length_ratio << "\n";
  out << "- min_similarity_to_reference: " << result.config.thresholds.min_similarity_to_reference << "\n\n";

  out << "Input summary:\n";
  out << "- FASTA-related files: " << result.input.files.size() << "\n";
  out << "- Raw sequences: " << result.summary.raw_sequence_count << "\n";
  out << "- Raw bases: " << result.summary.raw_total_bases << "\n\n";

  out << "Output summary:\n";
  out << "- Clean sequences: " << result.summary.clean_sequence_count << "\n";
  out << "- Clean bases: " << result.summary.clean_total_bases << "\n";
  out << "- Removed sequences: " << result.summary.removed_sequence_count << "\n";
  out << "- Collapsed duplicates: " << result.summary.collapsed_duplicate_count << "\n";
  out << "- Runtime ms: " << result.runtime.count() << "\n\n";

  out << "Reference:\n";
  out << "- Mode: " << result.reference.mode << "\n";
  out << "- Selected ID: " << (result.reference.selected_id.empty() ? "none" : result.reference.selected_id) << "\n";
  out << "- Reason: " << result.reference.selection_reason << "\n\n";

  out << "Removal reasons:\n";
  if (result.removal_counts.empty()) {
    out << "- none\n";
  } else {
    for (const auto& [reason, count] : result.removal_counts) {
      out << "- " << reason << ": " << count << "\n";
    }
  }
  out << "\nQC counts:\n";
  if (result.qc_counts.empty()) {
    out << "- none\n";
  } else {
    for (const auto& [name, count] : result.qc_counts) {
      out << "- " << name << ": " << count << "\n";
    }
  }
  out << "\nWarnings:\n";
  if (result.warning_counts.empty()) {
    out << "- none\n";
  } else {
    for (const auto& [name, count] : result.warning_counts) {
      out << "- " << name << ": " << count << "\n";
    }
  }
  out << "\nOutputs:\n";
  out << "- clean.fasta: " << (result.options.output_dir / result.config.output.clean_fasta_name) << "\n";
  out << "- result.json: " << (result.options.output_dir / "result.json") << "\n";
  out << "- summary.txt: " << (result.options.output_dir / "summary.txt") << "\n";
}

json error_result_json(const CliOptions& options,
                       const Config* config,
                       const PrepError& error,
                       int exit_code) {
  return {
      {"status", "error"},
      {"success", false},
      {"exit_code", exit_code},
      {"program_version", version()},
      {"command_line", command_line_to_json(options.command_line)},
      {"error_code", error.code()},
      {"error_message", error.what()},
      {"failed_stage", error.stage()},
      {"input", {{"path", options.input.string()}}},
      {"mode", config ? json(config->mode) : nullptr},
      {"strictness", config ? json(config->strictness) : nullptr},
      {"outputs",
       {
           {"clean_fasta", nullptr},
           {"result_json", "result.json"},
           {"summary_txt", "summary.txt"},
       }},
  };
}

void write_error_summary(const CliOptions& options,
                         const PrepError& error,
                         const fs::path& path) {
  auto out = open_output(path);
  out << "easymsa-prep " << version() << "\n";
  out << "Status: error\n";
  out << "Failed stage: " << error.stage() << "\n";
  out << "Error code: " << error.code() << "\n";
  out << "Error message: " << error.what() << "\n";
  if (!options.input.empty()) out << "Input: " << options.input << "\n";
}

}  // namespace

void write_success_outputs(const RunResult& result) {
  ensure_output_dir(result.options.output_dir);
  {
    auto out = open_output(result.options.output_dir / result.config.output.clean_fasta_name);
    for (const auto& record : result.records) {
      if (!record.removed && record.representative) {
        out << '>' << record.output_id << '\n'
            << wrap_fasta(record.clean_sequence, result.config.output.wrap_width);
      }
    }
  }
  {
    auto machine = open_output(result.options.output_dir / "sequence_map.jsonl");
    auto table = open_output(result.options.output_dir / "sequence_map.tsv");
    table << "input_index\toriginal_id\toriginal_header\toutput_id\tstatus\trepresentative_id\treasons\n";
    auto quote = [](const std::string& value) {
      std::string escaped = "\"";
      for (char c : value) { if (c == '\"') escaped += '\"'; escaped += c; }
      return escaped + "\"";
    };
    std::map<std::string, std::string> output_ids;
    for (const auto& record : result.records) output_ids[record.internal_id] = record.output_id;
    for (const auto& record : result.records) {
      const std::string status = record.removed ? "removed" : (record.representative ? "retained" : "collapsed");
      std::string representative;
      if (!record.representative && !record.removed) {
        representative = output_ids.at(record.representative_id);
      }
      auto reasons = record.remove_reasons;
      if (record.renamed) reasons.push_back(record.duplicate_id ? "duplicate_id" : "invalid_id");
      if (status == "collapsed") reasons.push_back("identical_sequence");
      const json entry = {{"input_index", record.record_index + 1}, {"original_id", record.original_id},
          {"original_header", record.original_header}, {"description", record.description},
          {"output_id", record.output_id}, {"status", status}, {"representative_id", representative}, {"reasons", reasons}};
      machine << entry.dump() << '\n';
      table << record.record_index + 1 << '\t' << quote(record.original_id) << '\t' << quote(record.original_header)
            << '\t' << quote(record.output_id) << '\t' << status << '\t' << quote(representative)
            << '\t' << quote(join(reasons, ";")) << '\n';
    }
  }
  {
    auto out = open_output(result.options.output_dir / "result.json");
    out << success_result_json(result).dump(2) << '\n';
  }
  write_summary_txt(result, result.options.output_dir / "summary.txt");
}

void write_error_outputs(const CliOptions& options,
                         const Config* config,
                         const PrepError& error,
                         int exit_code) {
  if (options.output_dir.empty()) return;
  ensure_output_dir(options.output_dir);
  {
    auto out = open_output(options.output_dir / "result.json");
    out << error_result_json(options, config, error, exit_code).dump(2) << '\n';
  }
  write_error_summary(options, error, options.output_dir / "summary.txt");
}

}  // namespace easymsa_prep
