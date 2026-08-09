#include "easymsa_prep.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

namespace easymsa_prep {
namespace {

constexpr std::size_t kKmerK = 5;
constexpr std::size_t kKmerUniverse = 1U << (2U * kKmerK);
constexpr std::size_t kKmerWords = kKmerUniverse / 64U;
constexpr std::size_t kLowMemoryThreshold = 10000;
constexpr std::size_t kLowMemorySampleSize = 2000;
constexpr std::size_t kMaximumSimilarityPairs = 2000;
constexpr std::uint64_t kFeatureSeed = 42;

using KmerBits = std::array<std::uint64_t, kKmerWords>;

class StableRandom {
 public:
  explicit StableRandom(std::uint64_t seed) : state_(seed) {}

  std::uint64_t next() {
    std::uint64_t value = (state_ += 0x9E3779B97F4A7C15ULL);
    value = (value ^ (value >> 30U)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27U)) * 0x94D049BB133111EBULL;
    return value ^ (value >> 31U);
  }

  std::uint64_t bounded(std::uint64_t bound) {
    if (bound == 0) throw std::invalid_argument("random bound must be positive");
    const std::uint64_t threshold = static_cast<std::uint64_t>(-bound) % bound;
    while (true) {
      const auto value = next();
      if (value >= threshold) return value % bound;
    }
  }

 private:
  std::uint64_t state_;
};

int dna_code(char base) {
  switch (base) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    default:
      return -1;
  }
}

KmerBits kmer_bits(const std::string& sequence) {
  KmerBits bits{};
  std::uint16_t rolling = 0;
  std::size_t valid = 0;
  constexpr std::uint16_t mask = static_cast<std::uint16_t>(kKmerUniverse - 1U);
  for (char base : sequence) {
    const int code = dna_code(base);
    if (code < 0) {
      rolling = 0;
      valid = 0;
      continue;
    }
    rolling = static_cast<std::uint16_t>(((rolling << 2U) | code) & mask);
    ++valid;
    if (valid >= kKmerK) {
      bits[rolling / 64U] |= 1ULL << (rolling % 64U);
    }
  }
  return bits;
}

std::size_t bit_count(const KmerBits& bits) {
  std::size_t count = 0;
  for (const auto word : bits) count += static_cast<std::size_t>(std::popcount(word));
  return count;
}

double jaccard(const KmerBits& left, const KmerBits& right) {
  std::size_t intersection = 0;
  std::size_t union_count = 0;
  for (std::size_t i = 0; i < kKmerWords; ++i) {
    intersection += static_cast<std::size_t>(std::popcount(left[i] & right[i]));
    union_count += static_cast<std::size_t>(std::popcount(left[i] | right[i]));
  }
  if (union_count == 0) return 1.0;
  return static_cast<double>(intersection) / static_cast<double>(union_count);
}

double mean(const std::vector<double>& values, double empty_value = 0.0) {
  if (values.empty()) return empty_value;
  long double total = 0.0L;
  for (const double value : values) total += static_cast<long double>(value);
  return static_cast<double>(total / static_cast<long double>(values.size()));
}

double sample_std(const std::vector<double>& values) {
  if (values.size() < 2) return 0.0;
  long double running_mean = 0.0L;
  long double m2 = 0.0L;
  std::size_t count = 0;
  for (const double raw : values) {
    ++count;
    const long double value = static_cast<long double>(raw);
    const long double delta = value - running_mean;
    running_mean += delta / static_cast<long double>(count);
    m2 += delta * (value - running_mean);
  }
  return static_cast<double>(
      std::sqrt(m2 / static_cast<long double>(values.size() - 1)));
}

double median(std::vector<double> values, double empty_value = 0.0) {
  if (values.empty()) return empty_value;
  std::sort(values.begin(), values.end());
  const auto midpoint = values.size() / 2;
  if (values.size() % 2 == 0) {
    return (values[midpoint - 1] + values[midpoint]) / 2.0;
  }
  return values[midpoint];
}

double iqr(std::vector<double> values) {
  if (values.size() < 2) return 0.0;
  std::sort(values.begin(), values.end());
  const auto midpoint = values.size() / 2;
  std::vector<double> lower(values.begin(), values.begin() + midpoint);
  std::vector<double> upper(values.end() - midpoint, values.end());
  return median(std::move(upper)) - median(std::move(lower));
}

std::vector<std::size_t> sample_record_indices(std::size_t sequence_count) {
  const std::size_t sample_size =
      sequence_count <= kLowMemoryThreshold
          ? sequence_count
          : std::min(sequence_count, kLowMemorySampleSize);
  std::vector<std::size_t> indices(sample_size);
  std::iota(indices.begin(), indices.end(), 0);
  if (sample_size == sequence_count) return indices;

  StableRandom rng(kFeatureSeed);
  for (std::size_t index = sample_size; index < sequence_count; ++index) {
    const auto replacement = static_cast<std::size_t>(rng.bounded(index + 1U));
    if (replacement < sample_size) indices[replacement] = index;
  }
  std::sort(indices.begin(), indices.end());
  return indices;
}

std::vector<std::uint64_t> sample_pair_ranks(std::size_t sequence_count) {
  if (sequence_count < 2) return {};
  const auto total_pairs =
      static_cast<std::uint64_t>(sequence_count) *
      static_cast<std::uint64_t>(sequence_count - 1U) / 2U;
  const auto sample_size = std::min<std::uint64_t>(
      total_pairs, static_cast<std::uint64_t>(kMaximumSimilarityPairs));
  std::vector<std::uint64_t> ranks;
  ranks.reserve(static_cast<std::size_t>(sample_size));
  if (total_pairs <= sample_size) {
    for (std::uint64_t rank = 0; rank < total_pairs; ++rank) ranks.push_back(rank);
    return ranks;
  }

  StableRandom rng(kFeatureSeed);
  std::unordered_set<std::uint64_t> selected;
  selected.reserve(static_cast<std::size_t>(sample_size * 2U));
  for (std::uint64_t current = total_pairs - sample_size; current < total_pairs;
       ++current) {
    const auto candidate = rng.bounded(current + 1U);
    if (!selected.insert(candidate).second) selected.insert(current);
  }
  ranks.assign(selected.begin(), selected.end());
  std::sort(ranks.begin(), ranks.end());
  return ranks;
}

std::pair<std::size_t, std::size_t> unrank_pair(std::uint64_t rank,
                                                std::size_t sequence_count) {
  for (std::size_t left = 0; left + 1U < sequence_count; ++left) {
    const auto row_size = static_cast<std::uint64_t>(sequence_count - left - 1U);
    if (rank < row_size) {
      return {left, left + 1U + static_cast<std::size_t>(rank)};
    }
    rank -= row_size;
  }
  throw std::out_of_range("similarity pair rank is outside the pair universe");
}

}  // namespace

DatasetFeatures compute_dataset_features(const std::vector<RecordResult>& records,
                                         SequenceType sequence_type) {
  const auto started = std::chrono::steady_clock::now();
  DatasetFeatures result;
  if (sequence_type == SequenceType::protein) {
    result.status = "unsupported";
    result.error_code = "UNSUPPORTED_SEQUENCE_TYPE";
    result.error_message =
        "The adaptive dataset feature profile supports DNA and RNA only.";
    result.runtime = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - started);
    return result;
  }

  try {
    std::vector<const RecordResult*> retained;
    retained.reserve(records.size());
    for (const auto& record : records) {
      if (!record.removed && record.representative) retained.push_back(&record);
    }
    if (retained.empty()) {
      throw std::runtime_error("No retained sequences are available for feature extraction.");
    }

    auto& values = result.values;
    values.n_sequences = retained.size();
    std::vector<double> lengths(retained.size());
    std::vector<double> gc_values(retained.size());
    std::vector<double> ambiguous_values(retained.size());

#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t raw_index = 0;
         raw_index < static_cast<std::ptrdiff_t>(retained.size()); ++raw_index) {
      const auto index = static_cast<std::size_t>(raw_index);
      const auto& sequence = retained[index]->clean_sequence;
      lengths[index] = static_cast<double>(sequence.size());
      std::size_t valid_bases = 0;
      std::size_t gc_bases = 0;
      std::size_t ambiguous_bases = 0;
      for (const char base : sequence) {
        const int code = dna_code(base);
        if (code < 0) {
          ++ambiguous_bases;
          continue;
        }
        ++valid_bases;
        if (base == 'G' || base == 'C') ++gc_bases;
      }
      gc_values[index] =
          valid_bases == 0
              ? 0.0
              : static_cast<double>(gc_bases) / static_cast<double>(valid_bases);
      ambiguous_values[index] =
          sequence.empty()
              ? 0.0
              : static_cast<double>(ambiguous_bases) /
                    static_cast<double>(sequence.size());
    }

    for (const double length : lengths) {
      values.total_bases += static_cast<std::uint64_t>(length);
    }
    values.length_min = *std::min_element(lengths.begin(), lengths.end());
    values.length_max = *std::max_element(lengths.begin(), lengths.end());
    values.length_mean = mean(lengths);
    values.length_median = median(lengths);
    values.length_std = sample_std(lengths);
    values.length_iqr = iqr(lengths);
    values.gc_mean = mean(gc_values);
    values.gc_std = sample_std(gc_values);
    values.ambiguous_mean = mean(ambiguous_values);

    const auto sampled_indices = sample_record_indices(retained.size());
    values.low_memory_mode = sampled_indices.size() < retained.size();
    values.kmer_sample_size = sampled_indices.size();
    values.kmer_sample_fraction =
        static_cast<double>(sampled_indices.size()) /
        static_cast<double>(retained.size());

    std::vector<KmerBits> bitsets(sampled_indices.size());
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t raw_index = 0;
         raw_index < static_cast<std::ptrdiff_t>(sampled_indices.size());
         ++raw_index) {
      const auto index = static_cast<std::size_t>(raw_index);
      bitsets[index] =
          kmer_bits(retained[sampled_indices[index]]->clean_sequence);
    }

    KmerBits all_kmers{};
    std::size_t total_per_sequence_distinct = 0;
    for (const auto& bits : bitsets) {
      total_per_sequence_distinct += bit_count(bits);
      for (std::size_t word = 0; word < kKmerWords; ++word) {
        all_kmers[word] |= bits[word];
      }
    }
    values.kmer_distinct = bit_count(all_kmers);
    values.kmer_diversity =
        total_per_sequence_distinct == 0
            ? 0.0
            : static_cast<double>(values.kmer_distinct) /
                  static_cast<double>(total_per_sequence_distinct);

    const auto pair_ranks = sample_pair_ranks(bitsets.size());
    std::vector<double> similarities(pair_ranks.size());
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t raw_index = 0;
         raw_index < static_cast<std::ptrdiff_t>(pair_ranks.size());
         ++raw_index) {
      const auto index = static_cast<std::size_t>(raw_index);
      const auto [left, right] =
          unrank_pair(pair_ranks[index], bitsets.size());
      similarities[index] = jaccard(bitsets[left], bitsets[right]);
    }
    values.similarity_pair_count = similarities.size();
    values.similarity_mean = mean(similarities, 1.0);
    values.similarity_median = median(similarities, 1.0);
    values.similarity_std = sample_std(similarities);

    result.status = "complete";
  } catch (const std::exception& error) {
    result.status = "error";
    result.error_code = "DATASET_FEATURE_EXTRACTION_FAILED";
    result.error_message = error.what();
  }
  result.runtime = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - started);
  return result;
}

}  // namespace easymsa_prep
