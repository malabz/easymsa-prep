# Filtering Rules

This document explains how to customize `easymsa-prep` with `--params`.

Run with a params file:

```bash
easymsa-prep input.fasta.gz -o result --params params.json --force
```

Unknown fields are rejected. This is intentional: spelling mistakes should fail
early instead of silently changing preprocessing behavior.

## Default Policy

The default policy is `mode: audit` and `strictness: normal`:

- normalize FASTA for downstream MSA: uppercase, DNA `U -> T`, remove gap
  characters, replace illegal DNA/RNA characters with `N`
- rename invalid or duplicated sequence IDs to stable IDs
- remove records that become empty or contain only unknown bases
- report high `N`, illegal-character ratio, length outliers, low similarity,
  and low complexity as warnings
- keep identical cleaned sequences unless you explicitly request collapse

This follows the behavior of common FASTA utilities: format transformation and
duplicate removal are explicit operations, while quality thresholds usually
need to be selected by the pipeline owner. Use `mode: filter` when you want
automatic removal, and use `strictness` to choose how tight the QC thresholds
should be.

Related tool behavior:

- SeqKit exposes gap removal, uppercasing, validation, and duplicate removal as
  explicit commands/options rather than one implicit default cleanup pipeline:
  <https://bioinf.shenwei.me/seqkit/usage/>
- VSEARCH quality filtering removes records only when users run filtering
  commands and select thresholds; several defaults are non-removing, for
  example minimum quality defaults to no discard:
  <https://vcru.wisc.edu/simonlab/bioinformatics/programs/vsearch/vsearch_manual.pdf>
- MAFFT generally tolerates ambiguous/unusual symbols through explicit modes,
  and its ambiguous-character removal threshold is described as an option rather
  than a general command-line default:
  <https://mafft.cbrc.jp/alignment/software/textcomparison.html>
  <https://mafft.cbrc.jp/alignment/software/changelog.html>

## Minimal Params

```json
{
  "mode": "audit",
  "strictness": "normal",
  "thresholds": {
    "max_n_ratio": 0.25
  },
  "rules": {
    "gaps": "remove",
    "identical_sequences": "collapse"
  }
}
```

Any missing field keeps the selected mode and strictness defaults.

## Full Shape

```json
{
  "mode": "audit",
  "strictness": "normal",
  "sequence_type": "dna",
  "thresholds": {
    "max_n_ratio": 0.2,
    "max_illegal_char_ratio": 0.01,
    "min_length_ratio": 0.7,
    "max_length_ratio": 1.3,
    "min_similarity_to_reference": 0.8,
    "max_low_complexity_ratio": 0.8
  },
  "rules": {
    "duplicate_ids": "rename_if_needed",
    "lowercase": "to_uppercase",
    "u_to_t": true,
    "gaps": "remove",
    "illegal_characters": "replace_with_N",
    "high_n_ratio": "warn_only",
    "high_illegal_char_ratio": "warn_only",
    "too_short": "warn_only",
    "too_long": "warn_only",
    "low_similarity_outlier": "warn_only",
    "low_complexity": "warn_only",
    "identical_sequences": "keep"
  },
  "limits": {
    "max_compressed_size": 1073741824,
    "max_uncompressed_size": 10737418240,
    "max_file_count": 10000,
    "max_fasta_file_count": 1000,
    "max_raw_sequence_count": 10000000,
    "max_unique_sequence_count": 10000000,
    "max_total_bases": 50000000000,
    "max_sequence_length": 10000000
  },
  "output": {
    "clean_fasta_name": "clean.fasta",
    "wrap_width": 80
  }
}
```

## Mode And Strictness

`mode` decides what happens after QC detects a problem. `strictness` decides
which threshold set is used.

| Mode | Use case |
| --- | --- |
| `audit` | Default mode; clean obvious formatting issues and report suspicious records |
| `filter` | Removal mode; remove records over QC thresholds and collapse exact duplicates |

| Strictness | Use case |
| --- | --- |
| `strict` | Tight thresholds for cleaner, more uniform inputs |
| `normal` | Default middle setting for most EasyMSA jobs |
| `lenient` | Wider thresholds for noisy or diverse datasets |

Default thresholds:

| Threshold | strict | normal | lenient |
| --- | ---: | ---: | ---: |
| `max_n_ratio` | 0.05 | 0.2 | 0.5 |
| `max_illegal_char_ratio` | 0.001 | 0.01 | 0.05 |
| `min_length_ratio` | 0.9 | 0.7 | 0.3 |
| `max_length_ratio` | 1.1 | 1.3 | 2.0 |
| `min_similarity_to_reference` | 0.95 | 0.8 | 0.5 |
| `max_low_complexity_ratio` | 0.6 | 0.8 | 0.95 |

Default actions:

| Rule | filter mode | audit mode |
| --- | --- | --- |
| `high_n_ratio` | remove | warn_only |
| `high_illegal_char_ratio` | remove | warn_only |
| `too_short` | remove | warn_only |
| `too_long` | remove | warn_only |
| `low_similarity_outlier` | remove | warn_only |
| `low_complexity` | remove | warn_only |
| `identical_sequences` | collapse | keep |

## Sequence Type

```json
{
  "sequence_type": "dna"
}
```

Supported values:

| Value | Behavior |
| --- | --- |
| `dna` | DNA preprocessing; `U` is converted to `T` when `u_to_t` is true |
| `rna` | RNA preprocessing; `T` is converted to `U` when `u_to_t` is true |
| `protein` | Reserved extension point; MVP behavior is not optimized for protein QC |

## Thresholds

### `max_n_ratio`

Maximum allowed ratio of unknown bases after cleaning.

```json
{
  "thresholds": {
    "max_n_ratio": 0.2
  }
}
```

If the ratio is exceeded, `high_n_ratio` decides whether to remove or warn.
All-`N` sequences are removed.

### `max_illegal_char_ratio`

Maximum ratio of illegal characters in the original sequence.

```json
{
  "thresholds": {
    "max_illegal_char_ratio": 0.01
  }
}
```

If the ratio is exceeded, `high_illegal_char_ratio` decides whether to remove
or warn. Small numbers of illegal characters can be replaced with `N`.

### `min_length_ratio` and `max_length_ratio`

Length outlier thresholds relative to median cleaned length.

```text
too_short: clean_length < median_length * min_length_ratio
too_long:  clean_length > median_length * max_length_ratio
```

Example:

```json
{
  "thresholds": {
    "min_length_ratio": 0.5,
    "max_length_ratio": 1.8
  }
}
```

### `min_similarity_to_reference`

Minimum k-mer similarity to the automatically selected reference sequence.
Sequences below this value are considered low-similarity outliers.

```json
{
  "thresholds": {
    "min_similarity_to_reference": 0.8
  }
}
```

The reference is selected from the first 1000 usable records. If no suitable
reference is available, low-similarity filtering is skipped and a warning is
reported.

### `max_low_complexity_ratio`

Maximum allowed dominance ratio for a single base/residue before a sequence is
treated as low complexity.

```json
{
  "thresholds": {
    "max_low_complexity_ratio": 0.8
  }
}
```

Low complexity also considers long homopolymers and low entropy.

## Rules

### ID Rules

```json
{
  "rules": {
    "duplicate_ids": "rename_if_needed"
  }
}
```

Allowed value:

- `rename_if_needed`

Invalid, empty, too long, or duplicate IDs are replaced with stable IDs like
`seq_000001`.

### Character Rules

```json
{
  "rules": {
    "lowercase": "to_uppercase",
    "u_to_t": true,
    "gaps": "remove",
    "illegal_characters": "replace_with_N"
  }
}
```

Allowed values:

| Rule | Values |
| --- | --- |
| `lowercase` | `to_uppercase`, `keep` |
| `u_to_t` | `true`, `false` |
| `gaps` | `remove`, `keep` |
| `illegal_characters` | `replace_with_N`, `remove_sequence`, `warn_only` |

Notes:

- For DNA, `u_to_t: true` converts `U` to `T`.
- For RNA, `u_to_t: true` converts `T` to `U`.
- `replace_with_N` writes `N` for illegal DNA/RNA characters.
- `remove_sequence` removes any record containing illegal characters.
- `warn_only` keeps illegal characters and reports warnings/counts.

### Removal Rules

```json
{
  "rules": {
    "high_n_ratio": "remove",
    "high_illegal_char_ratio": "remove",
    "too_short": "remove",
    "too_long": "remove",
    "low_similarity_outlier": "remove",
    "low_complexity": "remove"
  }
}
```

Allowed values for each rule:

- `remove`
- `warn_only`

Use `warn_only` when you want `result.json` to report issues without dropping
the sequence.

### Duplicate Sequence Rule

```json
{
  "rules": {
    "identical_sequences": "collapse"
  }
}
```

Allowed values:

- `collapse`: keep the first representative of each identical cleaned sequence.
- `keep`: keep all cleaned sequences.

## Limits

Limits protect backend workers from oversized inputs and archive bombs.

```json
{
  "limits": {
    "max_compressed_size": 1073741824,
    "max_uncompressed_size": 10737418240,
    "max_file_count": 10000,
    "max_fasta_file_count": 1000,
    "max_raw_sequence_count": 10000000,
    "max_unique_sequence_count": 10000000,
    "max_total_bases": 50000000000,
    "max_sequence_length": 10000000
  }
}
```

When a limit is exceeded, `easymsa-prep` exits with code `1`, does not write
`clean.fasta`, and writes an error `result.json` when possible.

## Output Settings

```json
{
  "output": {
    "clean_fasta_name": "clean.fasta",
    "wrap_width": 80
  }
}
```

`clean_fasta_name` must be a simple filename, not a path. `wrap_width` must be
between 1 and 1000.

## Common Recipes

### Keep More Borderline Sequences

```json
{
  "mode": "audit",
  "strictness": "lenient",
  "rules": {
    "low_similarity_outlier": "warn_only",
    "low_complexity": "warn_only"
  }
}
```

The default `audit` mode already reports these issues without removal.
Use this recipe when you also want wider thresholds.

### Filter Backend Cleanup

```json
{
  "mode": "filter",
  "strictness": "normal",
  "rules": {
    "gaps": "remove",
    "illegal_characters": "replace_with_N",
    "identical_sequences": "collapse"
  }
}
```

`filter` mode already removes high-N, high-illegal-character, length-outlier,
low-similarity, and low-complexity records, and collapses exact duplicates.

### Preserve Gaps

```json
{
  "mode": "audit",
  "strictness": "normal",
  "rules": {
    "gaps": "keep"
  }
}
```

### RNA Input

```json
{
  "mode": "audit",
  "strictness": "normal",
  "sequence_type": "rna",
  "rules": {
    "u_to_t": true
  }
}
```

### Collapse Exact Duplicates

```json
{
  "mode": "audit",
  "strictness": "normal",
  "rules": {
    "identical_sequences": "collapse"
  }
}
```

The representative is the first record in stable input order.

## Debugging Params

If params are invalid:

- exit code is `1`
- `clean.fasta` is not written
- `result.json` contains `error_code: "INVALID_PARAMS_JSON"`
- unknown fields are reported as errors

Start from a small params file, run once, then inspect:

```bash
cat result/result.json
cat result/summary.txt
```
