#include "common.hpp"

namespace easymsa_prep {
namespace {

double clean_median_length(const std::vector<RecordResult>& results) {
  std::vector<std::size_t> lengths;
  lengths.reserve(results.size());
  for (const auto& result : results) {
    if (result.clean_length > 0) lengths.push_back(result.clean_length);
  }
  if (lengths.empty()) return 0.0;
  std::sort(lengths.begin(), lengths.end());
  const auto mid = lengths.size() / 2;
  if (lengths.size() % 2 == 0) {
    return (static_cast<double>(lengths[mid - 1]) + static_cast<double>(lengths[mid])) / 2.0;
  }
  return static_cast<double>(lengths[mid]);
}

std::size_t select_reference(const std::vector<RecordResult>& results,
                             double median,
                             const Thresholds& thresholds) {
  std::size_t best = results.size();
  double best_score = -1.0;
  const std::size_t limit = std::min<std::size_t>(results.size(), kReferenceCandidateLimit);
  for (std::size_t i = 0; i < limit; ++i) {
    const auto& item = results[i];
    if (item.clean_length == 0 || item.low_complexity ||
        item.n_ratio > thresholds.max_n_ratio ||
        item.illegal_char_ratio > thresholds.max_illegal_char_ratio) {
      continue;
    }
    if (median > 0.0 &&
        (item.clean_length < median * thresholds.min_length_ratio ||
         item.clean_length > median * thresholds.max_length_ratio)) {
      continue;
    }
    const double len_score = median > 0.0
                                 ? 1.0 - std::min(1.0, std::abs(static_cast<double>(item.clean_length) - median) / median)
                                 : 1.0;
    const double score = len_score - item.n_ratio - item.illegal_char_ratio;
    if (score > best_score) {
      best = i;
      best_score = score;
    }
  }
  return best;
}

void add_reason_or_warning(RecordResult& result,
                           const std::string& reason,
                           const std::string& rule) {
  if (rule == "remove") {
    result.removed = true;
    result.remove_reasons.push_back(reason);
  } else if (rule == "warn_only") {
    result.warnings.push_back(reason);
  }
}

}  // namespace

RunResult run_pipeline(const CliOptions& options, const Config& config) {
  const auto start = std::chrono::steady_clock::now();
  configure_threads(options.threads);
  ensure_output_dir(options.output_dir);
  check_output_overwrite(options.output_dir, config.output, options.force);
  if (options.force) {
    std::error_code ec;
    fs::remove(options.output_dir / config.output.clean_fasta_name, ec);
  }

  RunResult result;
  result.options = options;
  result.config = config;
  result.input = read_input_bundle(options.input, config);
  result.summary.raw_sequence_count = result.input.records.size();
  for (const auto& record : result.input.records) {
    result.summary.raw_total_bases += record.sequence.size();
  }

  std::unordered_map<std::string, std::size_t> id_counts;
  for (const auto& record : result.input.records) ++id_counts[record.original_id];

  result.records.resize(result.input.records.size());
#pragma omp parallel for schedule(static)
  for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(result.input.records.size()); ++i) {
    const auto& record = result.input.records[static_cast<std::size_t>(i)];
    RecordResult item;
    item.record_index = record.record_index;
    item.internal_id = record.internal_id;
    item.original_id = record.original_id;
    item.original_header = record.original_header;
    item.description = record.description;
    item.output_id = record.original_id;
    item.source_file = record.source_file;
    item.original_length = record.sequence.size();

    std::size_t whitespace_count = 0;
    std::size_t illegal_count = 0;
    item.clean_sequence = normalize_sequence(record.sequence,
                                             config.sequence_type,
                                             config.rules,
                                             &item.gap_count,
                                             &illegal_count,
                                             &whitespace_count);
    item.clean_length = item.clean_sequence.size();
    item.illegal_char_ratio = item.original_length
                                  ? static_cast<double>(illegal_count) / static_cast<double>(item.original_length)
                                  : 1.0;
    const auto n_count = count_n_like(item.clean_sequence, config.sequence_type);
    item.n_ratio = item.clean_length ? static_cast<double>(n_count) / static_cast<double>(item.clean_length) : 1.0;
    item.low_complexity = is_low_complexity(item.clean_sequence, config.thresholds.max_low_complexity_ratio);
    item.invalid_id = !valid_id(item.original_id);
    item.hash = xxh64(item.clean_sequence);
    if (whitespace_count > 0) item.warnings.push_back("whitespace_removed");
    if (illegal_count > 0 && config.rules.illegal_characters == "remove_sequence") {
      item.removed = true;
      item.remove_reasons.push_back("illegal_characters");
    }
    if (item.clean_sequence.empty()) {
      item.removed = true;
      item.remove_reasons.push_back("empty_after_cleaning");
    }
    result.records[static_cast<std::size_t>(i)] = std::move(item);
  }

  std::unordered_set<std::string> reserved_ids;
  for (const auto& item : result.records) reserved_ids.insert(item.original_id);
  for (auto& item : result.records) {
    item.duplicate_id = id_counts[item.original_id] > 1;
    if (item.invalid_id || item.duplicate_id) {
      item.output_id = item.internal_id;
      std::size_t suffix = 1;
      while (reserved_ids.count(item.output_id)) {
        item.output_id = item.internal_id + "_" + std::to_string(suffix++);
      }
      reserved_ids.insert(item.output_id);
      item.renamed = true;
    }
  }

  const double median = clean_median_length(result.records);
  const std::size_t reference_index = select_reference(result.records, median, config.thresholds);
  std::string reference_sequence;
  if (reference_index < result.records.size()) {
    result.reference.available = true;
    result.reference.selected_id = result.records[reference_index].internal_id;
    result.reference.selection_reason = "close_to_median_length_low_N_low_illegal_chars";
    reference_sequence = result.records[reference_index].clean_sequence;
  } else {
    result.reference.available = false;
    result.reference.selection_reason = "no_suitable_reference";
    result.warning_counts["reference_unavailable"] = 1;
  }

#pragma omp parallel for schedule(static)
  for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(result.records.size()); ++i) {
    auto& item = result.records[static_cast<std::size_t>(i)];
    if (result.reference.available) {
      item.similarity_to_reference = kmer_similarity(item.clean_sequence, reference_sequence);
    }
  }

  std::unordered_set<std::string> used_output_ids;
  for (auto& item : result.records) {
    if (item.renamed) ++result.qc_counts["renamed_ids"];
    if (item.invalid_id) ++result.qc_counts["invalid_ids"];
    if (item.duplicate_id) ++result.qc_counts["duplicate_ids"];
    if (item.gap_count > 0) ++result.qc_counts["sequences_with_gaps"];
    if (item.illegal_char_ratio > 0.0) ++result.qc_counts["sequences_with_illegal_chars"];
    if (item.n_ratio >= 1.0 && item.clean_length > 0) {
      item.removed = true;
      item.remove_reasons.push_back("all_n");
      ++result.qc_counts["all_n"];
    }
    if (item.n_ratio > config.thresholds.max_n_ratio) {
      ++result.qc_counts["high_n_ratio"];
      add_reason_or_warning(item, "high_n_ratio", config.rules.high_n_ratio);
    }
    if (item.illegal_char_ratio > config.thresholds.max_illegal_char_ratio) {
      ++result.qc_counts["high_illegal_char_ratio"];
      add_reason_or_warning(item, "high_illegal_char_ratio", config.rules.high_illegal_char_ratio);
    }
    if (median > 0.0 && item.clean_length < median * config.thresholds.min_length_ratio) {
      ++result.qc_counts["too_short"];
      add_reason_or_warning(item, "too_short", config.rules.too_short);
    }
    if (median > 0.0 && item.clean_length > median * config.thresholds.max_length_ratio) {
      ++result.qc_counts["too_long"];
      add_reason_or_warning(item, "too_long", config.rules.too_long);
    }
    if (item.low_complexity) {
      ++result.qc_counts["low_complexity"];
      add_reason_or_warning(item, "low_complexity", config.rules.low_complexity);
    }
    if (result.reference.available &&
        item.record_index != result.records[reference_index].record_index &&
        item.similarity_to_reference < config.thresholds.min_similarity_to_reference) {
      ++result.qc_counts["low_similarity_outlier"];
      add_reason_or_warning(item, "low_similarity_outlier", config.rules.low_similarity_outlier);
    }
    if (used_output_ids.find(item.output_id) != used_output_ids.end()) {
      item.output_id = item.internal_id;
      item.renamed = true;
      ++result.qc_counts["renamed_ids"];
    }
    used_output_ids.insert(item.output_id);
  }

  if (config.rules.identical_sequences == "collapse") {
    std::unordered_map<std::string, std::size_t> representatives;
    representatives.reserve(result.records.size());
    for (auto& item : result.records) {
      if (item.removed) continue;
      const std::string key = std::to_string(item.hash) + ":" + item.clean_sequence;
      const auto [it, inserted] = representatives.emplace(key, item.record_index);
      if (!inserted) {
        item.representative = false;
        item.representative_id = format_internal_id(it->second + 1);
        ++result.summary.collapsed_duplicate_count;
      }
    }
    result.duplicate_group_count = representatives.size() < result.records.size()
                                       ? result.summary.collapsed_duplicate_count
                                       : 0;
    if (representatives.size() > config.limits.max_unique_sequence_count) {
      throw PrepError("LIMIT_EXCEEDED", "dedup", "unique sequence count exceeds max_unique_sequence_count");
    }
  }

  for (const auto& item : result.records) {
    for (const auto& warning : item.warnings) ++result.warning_counts[warning];
    if (item.gap_count > 0 && config.rules.gaps == "remove") {
      ++result.cleaning_counts["gaps_removed"];
    }
    if (item.illegal_char_ratio > 0.0 && config.rules.illegal_characters == "replace_with_N") {
      ++result.cleaning_counts["illegal_characters_replaced"];
    }
    if (item.removed) {
      ++result.summary.removed_sequence_count;
      if (item.remove_reasons.empty()) {
        ++result.removal_counts["removed"];
      } else {
        for (const auto& reason : item.remove_reasons) ++result.removal_counts[reason];
      }
      continue;
    }
    if (!item.representative) continue;
    ++result.summary.clean_sequence_count;
    result.summary.clean_total_bases += item.clean_sequence.size();
  }

  result.dataset_features =
      compute_dataset_features(result.records, config.sequence_type);
  result.runtime = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start);
  write_success_outputs(result);
  return result;
}

}  // namespace easymsa_prep
