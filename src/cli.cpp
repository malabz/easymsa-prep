#include "common.hpp"

#include <CLI/CLI.hpp>

namespace easymsa_prep {

CliOptions parse_cli(int argc, char** argv) {
  CLI::App app{"EasyMSA local FASTA preprocessing tool"};
  app.set_version_flag("--version", version());

  CliOptions options;
  std::string mode = "audit";
  std::string strictness = "normal";
  std::string input;

  app.add_option("input", input, "Input FASTA/GZ/XZ/ZIP/TAR archive file or directory")->required();
  app.add_option("-o,--output", options.output_dir, "Output directory")->required();
  app.add_option("--params", options.params_path, "Optional params JSON");
  app.add_option("--mode", mode, "audit or filter");
  app.add_option("--strictness", strictness, "strict, normal, or lenient");
  app.add_option("--threads", options.threads, "OpenMP worker threads");
  app.add_option("--batch-size", options.batch_size, "Records per analysis batch");
  app.add_flag("--force", options.force, "Overwrite generated output files");
  app.add_flag("--quiet", options.quiet, "Reduce terminal output");

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& ex) {
    const int code = app.exit(ex);
    if (code == 0) std::exit(0);
    throw PrepError("INVALID_CLI", "cli", ex.what());
  }

  options.input = input;
  options.mode = normalize_mode_name(mode);
  options.strictness = normalize_strictness_name(strictness);
  if (options.batch_size <= 0) {
    throw PrepError("INVALID_CLI", "cli", "--batch-size must be positive");
  }
  if (options.threads < 0) {
    throw PrepError("INVALID_CLI", "cli", "--threads must be non-negative");
  }
  for (int i = 0; i < argc; ++i) options.command_line.emplace_back(argv[i]);
  return options;
}

}  // namespace easymsa_prep
