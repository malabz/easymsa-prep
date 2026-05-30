# Examples

This folder contains small inputs for the simplified single-command
`easymsa-prep` workflow.

- `mixed_sequences.fasta`: plain FASTA with lowercase bases, U bases, gaps,
  illegal characters, high-N sequence, duplicate IDs, and duplicate sequences.
- `mixed_sequences.fasta.gz`: gzip-compressed FASTA.
- `mixed_sequences.zip`: zip archive input.
- `rna_sequences.fasta`: RNA input usable with a params JSON that sets
  `"sequence_type": "rna"`.

Run a single input:

```bash
../build/easymsa-prep mixed_sequences.fasta.gz -o output/gz --force --threads 4
```

Run with a params JSON:

```bash
../build/easymsa-prep rna_sequences.fasta -o output/rna --params rna_params.json --force
```

Check deterministic output:

```bash
../build/easymsa-prep mixed_sequences.fasta -o output/t1 --force --threads 1
../build/easymsa-prep mixed_sequences.fasta -o output/t8 --force --threads 8
cmp output/t1/clean.fasta output/t8/clean.fasta
```

Run all example scenarios from the repository root:

```bash
bash examples/run_examples.sh
```
