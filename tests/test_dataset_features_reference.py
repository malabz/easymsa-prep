#!/usr/bin/env python3
"""Cross-language reference test for easymsa-prep dataset features."""

from __future__ import annotations

import json
import math
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

MASK64 = (1 << 64) - 1
SEED = 42
MAX_PAIRS = 2000


class StableRandom:
    def __init__(self, seed: int) -> None:
        self.state = seed & MASK64

    def next(self) -> int:
        self.state = (self.state + 0x9E3779B97F4A7C15) & MASK64
        value = self.state
        value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
        value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
        return (value ^ (value >> 31)) & MASK64

    def bounded(self, bound: int) -> int:
        threshold = ((1 << 64) - bound) % bound
        while True:
            value = self.next()
            if value >= threshold:
                return value % bound


def pair_ranks(sequence_count: int) -> list[int]:
    total = sequence_count * (sequence_count - 1) // 2
    count = min(total, MAX_PAIRS)
    if total <= count:
        return list(range(total))
    rng = StableRandom(SEED)
    selected: set[int] = set()
    for current in range(total - count, total):
        candidate = rng.bounded(current + 1)
        if candidate in selected:
            selected.add(current)
        else:
            selected.add(candidate)
    return sorted(selected)


def unrank_pair(rank: int, sequence_count: int) -> tuple[int, int]:
    for left in range(sequence_count - 1):
        row_size = sequence_count - left - 1
        if rank < row_size:
            return left, left + 1 + rank
        rank -= row_size
    raise AssertionError("invalid pair rank")


def parse_fasta(path: Path) -> list[str]:
    sequences: list[str] = []
    current: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        else:
            current.append(raw.strip())
    if current:
        sequences.append("".join(current))
    return sequences


def iqr(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    ordered = sorted(values)
    midpoint = len(ordered) // 2
    return statistics.median(ordered[-midpoint:]) - statistics.median(
        ordered[:midpoint]
    )


def expected_values(sequences: list[str]) -> dict[str, float | int]:
    lengths = [float(len(sequence)) for sequence in sequences]
    gc_values: list[float] = []
    ambiguous_values: list[float] = []
    kmer_sets: list[set[str]] = []
    for sequence in sequences:
        valid = [base for base in sequence if base in "ACGT"]
        gc_values.append(
            sum(base in "GC" for base in valid) / len(valid) if valid else 0.0
        )
        ambiguous_values.append(
            sum(base not in "ACGT" for base in sequence) / len(sequence)
            if sequence
            else 0.0
        )
        kmer_sets.append(
            {
                sequence[index : index + 5]
                for index in range(max(0, len(sequence) - 4))
                if all(base in "ACGT" for base in sequence[index : index + 5])
            }
        )

    distinct = len(set().union(*kmer_sets))
    total_per_sequence = sum(len(items) for items in kmer_sets)
    similarities: list[float] = []
    for rank in pair_ranks(len(sequences)):
        left, right = unrank_pair(rank, len(sequences))
        union = kmer_sets[left] | kmer_sets[right]
        similarities.append(
            len(kmer_sets[left] & kmer_sets[right]) / len(union) if union else 1.0
        )

    return {
        "n_sequences": len(sequences),
        "total_bases": sum(len(sequence) for sequence in sequences),
        "length_min": min(lengths),
        "length_max": max(lengths),
        "length_mean": statistics.fmean(lengths),
        "length_median": statistics.median(lengths),
        "length_std": statistics.stdev(lengths) if len(lengths) > 1 else 0.0,
        "length_iqr": iqr(lengths),
        "gc_mean": statistics.fmean(gc_values),
        "gc_std": statistics.stdev(gc_values) if len(gc_values) > 1 else 0.0,
        "ambiguous_mean": statistics.fmean(ambiguous_values),
        "kmer_k": 5,
        "kmer_sample_size": len(sequences),
        "kmer_sample_fraction": 1.0,
        "kmer_distinct": distinct,
        "kmer_diversity": distinct / total_per_sequence if total_per_sequence else 0.0,
        "similarity_mean": statistics.fmean(similarities) if similarities else 1.0,
        "similarity_median": statistics.median(similarities) if similarities else 1.0,
        "similarity_std": (
            statistics.stdev(similarities) if len(similarities) > 1 else 0.0
        ),
        "similarity_pair_count": len(similarities),
        "low_memory_mode": 0,
    }


def make_sequences() -> list[str]:
    sequences: list[str] = []
    for sequence_index in range(80):
        state = sequence_index + 1
        chars: list[str] = []
        for position in range(48 + sequence_index % 19):
            state = (state * 1664525 + 1013904223) & 0xFFFFFFFF
            chars.append("ACGT"[(state >> 16) & 3])
            if position == 17 and sequence_index % 9 == 0:
                chars[-1] = "N"
        sequences.append("".join(chars))
    return sequences


def run(binary: Path, root: Path, threads: int) -> dict:
    input_path = root / "input.fasta"
    output_path = root / f"out-{threads}"
    completed = subprocess.run(
        [
            str(binary),
            str(input_path),
            "-o",
            str(output_path),
            "--threads",
            str(threads),
            "--force",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise AssertionError(completed.stdout + completed.stderr)
    return json.loads((output_path / "result.json").read_text(encoding="utf-8"))


def main() -> int:
    binary = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="easymsa-prep-reference-") as tmpdir:
        root = Path(tmpdir)
        sequences = make_sequences()
        (root / "input.fasta").write_text(
            "".join(
                f">seq_{index}\n{sequence}\n"
                for index, sequence in enumerate(sequences)
            ),
            encoding="utf-8",
        )
        one = run(binary, root, 1)
        eight = run(binary, root, 8)
        one_features = one["dataset_features"]
        eight_features = eight["dataset_features"]
        assert one_features["status"] == "complete"
        assert one_features["profile"] == "uniform-jaccard-2000-v1"
        assert one_features["values"] == eight_features["values"]
        observed = one_features["values"]
        expected = expected_values(parse_fasta(root / "out-1" / "clean.fasta"))
        for key, expected_value in expected.items():
            observed_value = observed[key]
            if isinstance(expected_value, int):
                assert observed_value == expected_value, (key, observed_value, expected_value)
            else:
                assert math.isclose(
                    float(observed_value),
                    expected_value,
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                ), (key, observed_value, expected_value)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

