#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "easymsa_prep.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

namespace {

void write_text(const fs::path& path, const std::string& text) {
  fs::create_directories(path.parent_path());
  std::ofstream out(path);
  out << text;
}

std::string read_text(const fs::path& path) {
  std::ifstream in(path);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

easymsa_prep::CliOptions options_for(const fs::path& input, const fs::path& output, int threads = 1) {
  easymsa_prep::CliOptions options;
  options.input = input;
  options.output_dir = output;
  options.mode = "audit";
  options.strictness = "normal";
  options.threads = threads;
  options.force = true;
  options.batch_size = 10000;
  options.command_line = {"easymsa-prep", input.string(), "-o", output.string()};
  return options;
}

}  // namespace

TEST_CASE("mode and strictness defaults") {
  CHECK(easymsa_prep::normalize_mode_name("audit") == "audit");
  CHECK(easymsa_prep::normalize_mode_name("filter") == "filter");
  CHECK(easymsa_prep::normalize_strictness_name("strict") == "strict");
  CHECK(easymsa_prep::normalize_strictness_name("normal") == "normal");
  CHECK(easymsa_prep::normalize_strictness_name("lenient") == "lenient");
  CHECK_THROWS_AS(easymsa_prep::normalize_strictness_name("balanced"), easymsa_prep::PrepError);
  CHECK_THROWS_AS(easymsa_prep::normalize_strictness_name("balance"), easymsa_prep::PrepError);

  const auto filter = easymsa_prep::default_config_for_mode("filter", "normal");
  const auto audit = easymsa_prep::default_config_for_mode("audit", "normal");
  CHECK(filter.mode == "filter");
  CHECK(filter.strictness == "normal");
  CHECK(filter.rules.high_n_ratio == "remove");
  CHECK(filter.rules.identical_sequences == "collapse");
  CHECK(audit.rules.high_n_ratio == "warn_only");
  CHECK(audit.rules.low_similarity_outlier == "warn_only");
  CHECK(audit.rules.identical_sequences == "keep");
  CHECK(easymsa_prep::default_config_for_mode("audit", "lenient").thresholds.max_n_ratio == doctest::Approx(0.5));
  CHECK(easymsa_prep::default_config_for_mode("audit", "strict").thresholds.max_n_ratio == doctest::Approx(0.05));
}

TEST_CASE("DNA and RNA normalization") {
  auto config = easymsa_prep::default_config_for_mode("audit", "normal");
  std::size_t gaps = 0;
  std::size_t illegal = 0;
  std::size_t spaces = 0;
  CHECK(easymsa_prep::normalize_sequence("acgu-xx n", easymsa_prep::SequenceType::dna,
                                         config.rules, &gaps, &illegal, &spaces) == "ACGTNNN");
  CHECK(gaps == 1);
  CHECK(illegal == 2);
  CHECK(spaces == 1);

  config.sequence_type = easymsa_prep::SequenceType::rna;
  CHECK(easymsa_prep::normalize_sequence("acgt-u", easymsa_prep::SequenceType::rna,
                                         config.rules, nullptr, nullptr, nullptr) == "ACGUU");
}

TEST_CASE("kmer similarity and low complexity") {
  CHECK(easymsa_prep::kmer_similarity("ACGTACGTACGT", "ACGTACGTACGT", 3) == doctest::Approx(1.0));
  CHECK(easymsa_prep::kmer_similarity("ACGTACGTACGT", "TTTTTTTTTTTT", 3) < 0.5);
  CHECK(easymsa_prep::is_low_complexity("AAAAAAAAAAAAAAAAAAAA", 0.8));
  CHECK_FALSE(easymsa_prep::is_low_complexity("ACGTACGTGCAATGCA", 0.8));
}

TEST_CASE("params JSON rejects unknown fields") {
  const auto root = fs::temp_directory_path() / "easymsa-prep-params-test";
  fs::remove_all(root);
  fs::create_directories(root);
  const auto params = root / "params.json";
  write_text(params, R"({"mode":"audit","strictness":"normal","unknown":true})");

  auto options = options_for(root / "input.fasta", root / "out");
  options.params_path = params;
  CHECK_THROWS_AS(easymsa_prep::load_config(options), easymsa_prep::PrepError);

  const auto legacy_params = root / "legacy_profile.json";
  write_text(legacy_params, R"({"profile":"strict"})");
  options.params_path = legacy_params;
  CHECK_THROWS_AS(easymsa_prep::load_config(options), easymsa_prep::PrepError);
  fs::remove_all(root);
}

TEST_CASE("pipeline writes simplified outputs and is deterministic") {
  const auto root = fs::temp_directory_path() / "easymsa-prep-pipeline-test";
  fs::remove_all(root);
  fs::create_directories(root);
  const auto input = root / "input.fasta";
  write_text(input,
             ">seqA\nACGTACGTACGTACGT\n"
             ">seqB\nacgu-acgtxxxx\n"
             ">seqC\nNNNNNNNNNNNNNNNN\n"
             ">seqA\nACGTACGTACGTACGT\n");

  auto config = easymsa_prep::default_config_for_mode("audit", "normal");
  auto one = options_for(input, root / "out1", 1);
  auto result1 = easymsa_prep::run_pipeline(one, config);
  CHECK(result1.success);
  CHECK(fs::exists(one.output_dir / "clean.fasta"));
  CHECK(fs::exists(one.output_dir / "result.json"));
  CHECK(fs::exists(one.output_dir / "summary.txt"));
  CHECK_FALSE(fs::exists(one.output_dir / "recommended_strategy.json"));

  auto eight = options_for(input, root / "out8", 8);
  auto result8 = easymsa_prep::run_pipeline(eight, config);
  CHECK(result8.success);
  CHECK(read_text(one.output_dir / "clean.fasta") == read_text(eight.output_dir / "clean.fasta"));

  fs::remove_all(root);
}

TEST_CASE("invalid FASTA surfaces as PrepError") {
  const auto root = fs::temp_directory_path() / "easymsa-prep-invalid-test";
  fs::remove_all(root);
  fs::create_directories(root);
  const auto input = root / "bad.fasta";
  write_text(input, ">bad\n");

  auto config = easymsa_prep::default_config_for_mode("audit", "normal");
  auto options = options_for(input, root / "out", 1);
  CHECK_THROWS_AS(easymsa_prep::run_pipeline(options, config), easymsa_prep::PrepError);

  fs::remove_all(root);
}
