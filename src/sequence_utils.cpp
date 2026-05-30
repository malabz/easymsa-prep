#include "common.hpp"

namespace easymsa_prep {

std::string normalize_sequence(const std::string& sequence,
                               SequenceType type,
                               const Rules& rules,
                               std::size_t* gap_count,
                               std::size_t* illegal_count,
                               std::size_t* whitespace_count) {
  if (gap_count) *gap_count = 0;
  if (illegal_count) *illegal_count = 0;
  if (whitespace_count) *whitespace_count = 0;

  std::string out;
  out.reserve(sequence.size());
  for (unsigned char raw : sequence) {
    if (std::isspace(raw)) {
      if (whitespace_count) ++*whitespace_count;
      continue;
    }
    char c = static_cast<char>(raw);
    if (rules.lowercase == "to_uppercase") {
      c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    if (type == SequenceType::dna && rules.u_to_t && c == 'U') c = 'T';
    if (type == SequenceType::rna && rules.u_to_t && c == 'T') c = 'U';

    if (c == '-') {
      if (gap_count) ++*gap_count;
      if (rules.gaps == "keep") out.push_back(c);
      continue;
    }

    if (!is_legal_base(c, type)) {
      if (illegal_count) ++*illegal_count;
      if (rules.illegal_characters == "replace_with_N") {
        out.push_back(replacement_char(type));
      } else if (rules.illegal_characters == "warn_only") {
        out.push_back(c);
      }
      continue;
    }
    out.push_back(c);
  }
  return out;
}

bool valid_id(const std::string& id) {
  if (id.empty() || id.size() > kMaxIdLength) return false;
  return std::all_of(id.begin(), id.end(), [](unsigned char c) {
    return std::isalnum(c) || c == '_' || c == '-' || c == '.' || c == ':' || c == '|';
  });
}

bool is_low_complexity(const std::string& sequence, double max_ratio) {
  if (sequence.empty()) return true;
  std::array<std::size_t, 256> counts{};
  std::size_t longest_run = 0;
  std::size_t current_run = 0;
  char previous = '\0';
  for (char c : sequence) {
    ++counts[static_cast<unsigned char>(c)];
    current_run = (c == previous) ? current_run + 1 : 1;
    longest_run = std::max(longest_run, current_run);
    previous = c;
  }
  const auto max_count = *std::max_element(counts.begin(), counts.end());
  if (static_cast<double>(max_count) / static_cast<double>(sequence.size()) >= max_ratio) return true;
  if (longest_run >= std::max<std::size_t>(20, sequence.size() / 2)) return true;

  double entropy = 0.0;
  for (std::size_t count : counts) {
    if (!count) continue;
    const double p = static_cast<double>(count) / static_cast<double>(sequence.size());
    entropy -= p * std::log2(p);
  }
  return entropy < 1.0;
}

double kmer_similarity(const std::string& lhs, const std::string& rhs, std::size_t k) {
  if (lhs.empty() || rhs.empty()) return 0.0;
  if (lhs.size() < k || rhs.size() < k) return lhs == rhs ? 1.0 : 0.0;

  std::unordered_set<std::uint64_t> left;
  left.reserve(lhs.size());
  for (std::size_t i = 0; i + k <= lhs.size(); ++i) {
    left.insert(XXH64(lhs.data() + i, k, 0));
  }
  if (left.empty()) return 0.0;

  std::unordered_set<std::uint64_t> seen;
  seen.reserve(rhs.size());
  std::size_t shared = 0;
  for (std::size_t i = 0; i + k <= rhs.size(); ++i) {
    const auto h = XXH64(rhs.data() + i, k, 0);
    if (left.find(h) != left.end() && seen.insert(h).second) ++shared;
  }
  return static_cast<double>(shared) / static_cast<double>(left.size());
}

std::string sequence_type_name(SequenceType type) {
  switch (type) {
    case SequenceType::dna:
      return "dna";
    case SequenceType::rna:
      return "rna";
    case SequenceType::protein:
      return "protein";
  }
  return "dna";
}

SequenceType parse_sequence_type(const std::string& value) {
  const auto text = lower(value);
  if (text == "dna") return SequenceType::dna;
  if (text == "rna") return SequenceType::rna;
  if (text == "protein") return SequenceType::protein;
  throw PrepError("INVALID_PARAMS_JSON", "params", "invalid sequence_type: " + value);
}

}  // namespace easymsa_prep
