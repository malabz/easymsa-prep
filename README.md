# easymsa-prep

`easymsa-prep` cleans FASTA input before EasyMSA alignment. It accepts a FASTA
file, compressed FASTA, directory, or archive, then writes:

```text
clean.fasta
result.json
summary.txt
```

No subcommands are needed.

## Quick Install

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install -y cmake g++ libomp-dev zlib1g-dev liblzma-dev
```

Build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Optional local install:

```bash
cmake --install build --prefix "$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
```

Check:

```bash
easymsa-prep --help
```

If you do not install it, use `build/easymsa-prep` in the commands below.

## Quick Run

Run with defaults:

```bash
easymsa-prep input.fasta.gz -o result
```

The default is `--mode audit --strictness normal`: it normalizes sequence text,
renames unsafe IDs, removes empty/all-unknown records, and reports quality
outliers as warnings instead of deleting them. Use `--mode filter` when the
backend should automatically drop high-N, high-illegal-character, length
outlier, low-complexity, and low-similarity records. Use `--strictness strict`
or `--strictness lenient` to change QC thresholds.

Overwrite an existing result directory:

```bash
easymsa-prep input.fasta.gz -o result --force
```

Use more threads:

```bash
easymsa-prep input.fasta.gz -o result --threads 8 --force
```

Run the bundled examples:

```bash
bash examples/run_examples.sh
```

## Inputs

Supported inputs:

- FASTA: `.fasta`, `.fa`, `.fna`, `.ffn`, `.txt`
- Compressed FASTA: `.gz`, `.xz`
- Archives: `.zip`, `.tar.gz`, `.tgz`, `.tar.xz`
- Directories containing supported FASTA files

Archive members are checked for path traversal, absolute paths, hidden/system
files, nested archives, excessive file count, and excessive uncompressed size.

## Outputs

Success:

```text
result/
  clean.fasta
  result.json
  summary.txt
```

Failure:

```text
result/
  result.json
  summary.txt
```

`clean.fasta` is not produced on failure.

Backend workers should pass only `clean.fasta` to the MSA stage and treat
`result.json` as the machine-readable source of truth.

Successful DNA/RNA runs also expose a public `dataset_features` block in
`result.json`. It is computed from the final retained representatives and
contains length, GC, ambiguous-base, 5-mer, and uniformly sampled Jaccard
statistics. Mash remains an adaptive-selector concern and is not run by the
preprocessor. The feature profile is versioned as
`uniform-jaccard-2000-v1`; protein inputs report the feature block as
`unsupported`.

## Common Options

```bash
easymsa-prep <input> -o <output_dir> [options]
```

- `--mode audit|filter`: default `audit`.
- `--strictness strict|normal|lenient`: default `normal`.
- `--params params.json`: override thresholds, rules, limits, or output settings.
- `--threads N`: OpenMP threads; default is `min(OpenMP max threads, 8)`.
- `--batch-size N`: default `10000`.
- `--force`: overwrite generated output files.
- `--quiet`: reduce terminal output.
- `--version`, `--help`.

Mode summary:

| Mode | Default behavior |
| --- | --- |
| `audit` | Check and report suspicious records without dropping them |
| `filter` | Check, remove records that exceed QC thresholds, and collapse exact duplicates |

Strictness summary:

| Strictness | Threshold behavior |
| --- | --- |
| `strict` | Tight thresholds for cleaner, more uniform inputs |
| `normal` | Default middle setting for most runs |
| `lenient` | Wider thresholds for noisy or diverse inputs |

## Custom Filtering

Use `--params` to customize filtering:

```bash
easymsa-prep input.fasta.gz -o result --params params.json --force
```

Minimal example:

```json
{
  "mode": "filter",
  "strictness": "normal",
  "thresholds": {
    "max_n_ratio": 0.25,
    "max_illegal_char_ratio": 0.02
  },
  "rules": {
    "gaps": "remove",
    "illegal_characters": "replace_with_N",
    "high_n_ratio": "remove",
    "too_short": "remove",
    "too_long": "warn_only",
    "identical_sequences": "collapse"
  }
}
```

Detailed documentation:

- [Filtering Rules](docs/filtering-rules.md)

## Exit Codes

- `0`: success; all three output files exist.
- `1`: input, params, archive, FASTA, limit, or runtime error.
- `2`: command-line error.

When the output directory is known, errors write `result.json` and `summary.txt`.

## Development

Build and test:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Source files are intentionally flat under `src/`:

- `cli.cpp`: command line parsing
- `config.cpp`: mode, strictness, and params JSON
- `input.cpp`: discovery, archive extraction, FASTA reading
- `engine.cpp`: QC, cleaning, reference selection, dedup
- `reports.cpp`: `result.json` and `summary.txt`
- `sequence_utils.cpp`, `common.cpp`: shared helpers

Vendored dependencies: `CLI11`, `miniz`, `doctest`, `kseq.h`, `json.hpp`,
- `dataset_features.cpp`: deterministic non-Mash dataset statistics
`xxhash`. System dependencies: OpenMP, zlib, liblzma.
