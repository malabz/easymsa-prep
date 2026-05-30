#include "easymsa_prep.hpp"

#include <exception>
#include <iostream>

int main(int argc, char** argv) {
  easymsa_prep::CliOptions options;
  std::optional<easymsa_prep::Config> config;
  try {
    options = easymsa_prep::parse_cli(argc, argv);
    config = easymsa_prep::load_config(options);
    const auto result = easymsa_prep::run_pipeline(options, *config);
    if (!options.quiet) {
      std::cout << "Wrote " << result.summary.clean_sequence_count
                << " clean sequences to " << options.output_dir << "\n";
    }
    return result.exit_code;
  } catch (const easymsa_prep::PrepError& error) {
    int exit_code = error.code() == "INVALID_CLI" ? 2 : 1;
    try {
      if (options.command_line.empty()) {
        for (int i = 0; i < argc; ++i) options.command_line.emplace_back(argv[i]);
      }
      easymsa_prep::write_error_outputs(options, config ? &*config : nullptr, error, exit_code);
    } catch (...) {
    }
    std::cerr << "error: " << error.what() << "\n";
    return exit_code;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 1;
  }
}
