#!/usr/bin/env python3
"""Reformat FLORA-seq FASTQs into Cell Ranger-compatible read pairs.

This module implements the sequence-coordinate logic of the processing
workflow while reading FASTQ files as streams. Sequencer R1 is written as Cell
Ranger R2. Cell Ranger R1 is reconstructed from the reverse-complemented
16-bp droplet barcode and the 10-bp UMI encoded in sequencer R2.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
from collections import Counter
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class IndexRecord:
    sample: str
    index: str
    index_length: int
    sample_number: str = "S1"
    primer: str = "both"


@dataclass(frozen=True)
class FastqRecord:
    header: str
    sequence: str
    plus: str
    quality: str


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return path.open(mode)


def read_fastq(path: Path) -> Iterator[FastqRecord]:
    with open_text(path, "r") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"Incomplete FASTQ record in {path}")
            record = FastqRecord(
                header.rstrip("\r\n"),
                sequence.rstrip("\r\n"),
                plus.rstrip("\r\n"),
                quality.rstrip("\r\n"),
            )
            if not record.header.startswith("@") or not record.plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record in {path}: {record.header}")
            if len(record.sequence) != len(record.quality):
                raise ValueError(f"Sequence/quality length mismatch in {path}: {record.header}")
            yield record


def read_name(header: str) -> str:
    return header.split()[0].removeprefix("@").removesuffix("/1").removesuffix("/2")


def paired_records(r1_path: Path, r2_path: Path) -> Iterator[tuple[FastqRecord, FastqRecord]]:
    r1_iter = read_fastq(r1_path)
    r2_iter = read_fastq(r2_path)
    position = 0
    while True:
        try:
            r1 = next(r1_iter)
        except StopIteration:
            r1 = None
        try:
            r2 = next(r2_iter)
        except StopIteration:
            r2 = None
        if r1 is None and r2 is None:
            return
        position += 1
        if r1 is None or r2 is None:
            raise ValueError("R1 and R2 contain different numbers of records")
        if read_name(r1.header) != read_name(r2.header):
            raise ValueError(f"R1/R2 name mismatch at record {position}: {r1.header} != {r2.header}")
        yield r1, r2


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def load_indices(path: Path) -> list[IndexRecord]:
    records: list[IndexRecord] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "index", "index_length"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: sample, index, index_length")
        for row in reader:
            index = row["index"].strip().upper()
            index_length = int(row["index_length"])
            if len(index) != index_length:
                raise ValueError(f"Index length mismatch for {row['sample']}: {index}")
            records.append(
                IndexRecord(
                    sample=row["sample"].strip(),
                    index=index,
                    index_length=index_length,
                    sample_number=(row.get("sample_number") or "S1").strip(),
                    primer=(row.get("primer") or "both").strip(),
                )
            )
            if records[-1].primer not in {"3", "r", "both"}:
                raise ValueError(f"primer must be 3, r or both for {records[-1].sample}")
    if not records:
        raise ValueError(f"No index records found in {path}")
    return records


def candidate_assignments(sequence: str, indices: Iterable[IndexRecord]) -> list[tuple[IndexRecord, str]]:
    """Apply exact index matching with 3-prime priority within each sample."""
    candidates: list[tuple[IndexRecord, str]] = []
    records_by_sample: dict[str, list[IndexRecord]] = {}
    for record in indices:
        records_by_sample.setdefault(record.sample, []).append(record)
    for sample_records in records_by_sample.values():
        three_prime_matches = [
            (record, "3")
            for record in sample_records
            if record.primer in {"3", "both"}
            and sequence[16 : 16 + record.index_length] == record.index
        ]
        if three_prime_matches:
            candidates.extend(three_prime_matches)
            continue
        candidates.extend(
            (record, "r")
            for record in sample_records
            if record.primer in {"r", "both"}
            and sequence[26 : 26 + record.index_length] == record.index
        )
    return candidates


def reconstruct_cellranger_r1(
    raw_r2: FastqRecord,
    index_length: int,
    primer: str,
    tail_mode: str,
    cell_barcode_orientation: str,
) -> FastqRecord:
    sequence = raw_r2.sequence
    quality = raw_r2.quality
    minimum = 16 + index_length + 10 if primer == "3" else 26 + index_length
    if len(sequence) < minimum:
        raise ValueError(f"Barcode read is shorter than required ({minimum} bp): {raw_r2.header}")

    if cell_barcode_orientation == "reverse-complement":
        cell_barcode = reverse_complement(sequence[:16])
        cell_quality = quality[:16][::-1]
    else:
        cell_barcode = sequence[:16]
        cell_quality = quality[:16]
    if primer == "3":
        umi_start = 16 + index_length
        umi_end = umi_start + 10
        umi = sequence[umi_start:umi_end]
        umi_quality = quality[umi_start:umi_end]
        if tail_mode == "keep":
            out_sequence = cell_barcode + sequence[umi_start:]
            out_quality = cell_quality + quality[umi_start:]
        else:
            out_sequence = cell_barcode + umi
            out_quality = cell_quality + umi_quality
    else:
        umi = sequence[16:26]
        umi_quality = quality[16:26]
        if tail_mode == "keep":
            out_sequence = cell_barcode + umi + sequence[26 + index_length :]
            out_quality = cell_quality + umi_quality + quality[26 + index_length :]
        else:
            out_sequence = cell_barcode + umi
            out_quality = cell_quality + umi_quality

    return FastqRecord(raw_r2.header, out_sequence, raw_r2.plus, out_quality)


def write_record(handle: TextIO, record: FastqRecord) -> None:
    handle.write(f"{record.header}\n{record.sequence}\n{record.plus}\n{record.quality}\n")


def output_paths(
    outdir: Path,
    index: IndexRecord,
    primer: str,
    lane: str,
    gzip_output: bool,
) -> tuple[Path, Path]:
    library = f"{index.sample}_{primer}"
    suffix = ".fastq.gz" if gzip_output else ".fastq"
    library_dir = outdir / library
    r1 = library_dir / f"{library}_{index.sample_number}_{lane}_R1_001{suffix}"
    r2 = library_dir / f"{library}_{index.sample_number}_{lane}_R2_001{suffix}"
    return r1, r2


def run_reformat(args: argparse.Namespace) -> None:
    indices = load_indices(args.indices)
    args.outdir.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    output_counts: Counter[str] = Counter()

    with ExitStack() as stack:
        handles: dict[tuple[str, str], tuple[TextIO, TextIO]] = {}
        output_records = {
            (record.sample, primer): record
            for record in indices
            for primer in ("3", "r")
            if record.primer in {primer, "both"}
            and (args.only_primer is None or primer == args.only_primer)
        }
        for (sample, primer), index in output_records.items():
            r1_path, r2_path = output_paths(
                args.outdir, index, primer, args.lane, args.gzip_output
            )
            r1_path.parent.mkdir(parents=True, exist_ok=True)
            handles[(sample, primer)] = (
                stack.enter_context(open_text(r1_path, "w")),
                stack.enter_context(open_text(r2_path, "w")),
            )

        for raw_r1, raw_r2 in paired_records(args.raw_r1, args.raw_r2):
            counts["input_pairs"] += 1
            candidates = candidate_assignments(raw_r2.sequence, indices)
            if not candidates:
                counts["unassigned_pairs"] += 1
                continue
            if len(candidates) > 1:
                counts["multi_assigned_input_pairs"] += 1
            selected = [
                candidate
                for candidate in candidates
                if args.only_primer is None or candidate[1] == args.only_primer
            ]
            if not selected:
                counts["excluded_other_primer_pairs"] += 1
                continue
            for index, primer in selected:
                tail_mode = (
                    args.three_prime_tail_mode
                    if primer == "3" and args.three_prime_tail_mode
                    else args.random_tail_mode
                    if primer == "r" and args.random_tail_mode
                    else args.tail_mode
                )
                cellranger_r1 = reconstruct_cellranger_r1(
                    raw_r2,
                    index.index_length,
                    primer,
                    tail_mode,
                    args.cell_barcode_orientation,
                )
                r1_handle, r2_handle = handles[(index.sample, primer)]
                write_record(r1_handle, cellranger_r1)
                write_record(r2_handle, raw_r1)
                output_counts[f"{index.sample}_{primer}"] += 1
                counts["written_pairs"] += 1

    report = {
        "raw_r1": str(args.raw_r1),
        "raw_r2": str(args.raw_r2),
        "indices": str(args.indices),
        "tail_mode": args.tail_mode,
        "three_prime_tail_mode": args.three_prime_tail_mode or args.tail_mode,
        "random_tail_mode": args.random_tail_mode or args.tail_mode,
        "cell_barcode_orientation": args.cell_barcode_orientation,
        "counts": dict(sorted(counts.items())),
        "output_pairs": dict(sorted(output_counts.items())),
    }
    report_path = args.outdir / "reformat_summary.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


def raw_output_paths(
    outdir: Path,
    index: IndexRecord,
    primer: str,
    gzip_output: bool,
) -> tuple[Path, Path]:
    library = f"{index.sample}_{primer}"
    suffix = ".fastq.gz" if gzip_output else ".fastq"
    library_dir = outdir / library
    return (
        library_dir / f"{library}_sequencer_R1{suffix}",
        library_dir / f"{library}_sequencer_R2{suffix}",
    )


def run_split_raw(args: argparse.Namespace) -> None:
    indices = load_indices(args.indices)
    args.outdir.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    output_counts: Counter[str] = Counter()

    with ExitStack() as stack:
        handles: dict[tuple[str, str], tuple[TextIO, TextIO]] = {}
        output_records = {
            (record.sample, primer): record
            for record in indices
            for primer in ("3", "r")
            if record.primer in {primer, "both"}
            and (args.only_primer is None or primer == args.only_primer)
        }
        for (sample, primer), index in output_records.items():
            r1_path, r2_path = raw_output_paths(
                args.outdir, index, primer, args.gzip_output
            )
            r1_path.parent.mkdir(parents=True, exist_ok=True)
            handles[(sample, primer)] = (
                stack.enter_context(open_text(r1_path, "w")),
                stack.enter_context(open_text(r2_path, "w")),
            )

        for raw_r1, raw_r2 in paired_records(args.raw_r1, args.raw_r2):
            counts["input_pairs"] += 1
            candidates = candidate_assignments(raw_r2.sequence, indices)
            if not candidates:
                counts["unassigned_pairs"] += 1
                continue
            if len(candidates) > 1:
                counts["multi_assigned_input_pairs"] += 1
            selected = [
                candidate
                for candidate in candidates
                if args.only_primer is None or candidate[1] == args.only_primer
            ]
            if not selected:
                counts["excluded_other_primer_pairs"] += 1
                continue
            for index, primer in selected:
                r1_handle, r2_handle = handles[(index.sample, primer)]
                write_record(r1_handle, raw_r1)
                write_record(r2_handle, raw_r2)
                output_counts[f"{index.sample}_{primer}"] += 1
                counts["written_pairs"] += 1

    report = {
        "raw_r1": str(args.raw_r1),
        "raw_r2": str(args.raw_r2),
        "indices": str(args.indices),
        "counts": dict(sorted(counts.items())),
        "output_pairs": dict(sorted(output_counts.items())),
    }
    report_path = args.outdir / "split_raw_summary.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


def stream_copy_fastq(inputs: list[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(output, "w") as out_handle:
        for path in inputs:
            with open_text(path, "r") as in_handle:
                shutil.copyfileobj(in_handle, out_handle, length=1024 * 1024)


def run_merge(args: argparse.Namespace) -> None:
    if len(args.r1) != len(args.r2):
        raise ValueError("The number of R1 and R2 inputs must match")
    stream_copy_fastq(args.r1, args.output_r1)
    stream_copy_fastq(args.r2, args.output_r2)
    print(f"Merged {len(args.r1)} FASTQ pairs into {args.output_r1} and {args.output_r2}")


def run_truncate(args: argparse.Namespace) -> None:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.output, "w") as out_handle:
        count = 0
        for record in read_fastq(args.input):
            if len(record.sequence) < args.length:
                raise ValueError(f"Read shorter than {args.length} bp: {record.header}")
            write_record(
                out_handle,
                FastqRecord(
                    record.header,
                    record.sequence[: args.length],
                    record.plus,
                    record.quality[: args.length],
                ),
            )
            count += 1
    print(f"Truncated {count} reads to {args.length} bp: {args.output}")


def run_reverse_barcode(args: argparse.Namespace) -> None:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.output, "w") as out_handle:
        count = 0
        for record in read_fastq(args.input):
            if len(record.sequence) < 16:
                raise ValueError(f"Read shorter than 16 bp: {record.header}")
            write_record(
                out_handle,
                FastqRecord(
                    record.header,
                    reverse_complement(record.sequence[:16]) + record.sequence[16:],
                    record.plus,
                    record.quality[:16][::-1] + record.quality[16:],
                ),
            )
            count += 1
    print(f"Reverse-complemented 16-nt cell barcodes in {count} reads: {args.output}")


def length_summary(path: Path) -> tuple[int, Counter[int]]:
    counts: Counter[int] = Counter()
    total = 0
    for record in read_fastq(path):
        total += 1
        counts[len(record.sequence)] += 1
    return total, counts


def run_inspect(args: argparse.Namespace) -> None:
    rows = []
    for path in args.fastq:
        total, lengths = length_summary(path)
        rows.append(
            {
                "fastq": str(path),
                "reads": total,
                "length_counts": ",".join(f"{k}:{v}" for k, v in sorted(lengths.items())),
            }
        )
    writer = csv.DictWriter(
        __import__("sys").stdout,
        fieldnames=["fastq", "reads", "length_counts"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    reformat = subparsers.add_parser("reformat", help="Demultiplex and reconstruct Cell Ranger reads")
    reformat.add_argument("--raw-r1", type=Path, required=True)
    reformat.add_argument("--raw-r2", type=Path, required=True)
    reformat.add_argument("--indices", type=Path, required=True)
    reformat.add_argument("--outdir", type=Path, required=True)
    reformat.add_argument("--lane", default="L001")
    reformat.add_argument("--only-primer", choices=["3", "r"])
    reformat.add_argument(
        "--cell-barcode-orientation",
        choices=["raw", "reverse-complement"],
        default="reverse-complement",
    )
    reformat.add_argument("--tail-mode", choices=["keep", "trim"], default="trim")
    reformat.add_argument(
        "--three-prime-tail-mode",
        choices=["keep", "trim"],
        help="Override --tail-mode for oligo(dT)-primed reads",
    )
    reformat.add_argument(
        "--random-tail-mode",
        choices=["keep", "trim"],
        help="Override --tail-mode for random-primed reads",
    )
    reformat.add_argument(
        "--gzip-output", action="store_true"
    )
    reformat.set_defaults(func=run_reformat)

    split_raw = subparsers.add_parser(
        "split-raw", help="Demultiplex while retaining sequencer read orientation"
    )
    split_raw.add_argument("--raw-r1", type=Path, required=True)
    split_raw.add_argument("--raw-r2", type=Path, required=True)
    split_raw.add_argument("--indices", type=Path, required=True)
    split_raw.add_argument("--outdir", type=Path, required=True)
    split_raw.add_argument("--only-primer", choices=["3", "r"])
    split_raw.add_argument("--gzip-output", action="store_true")
    split_raw.set_defaults(func=run_split_raw)

    merge = subparsers.add_parser("merge", help="Merge corresponding FASTQ pairs")
    merge.add_argument("--r1", type=Path, nargs="+", required=True)
    merge.add_argument("--r2", type=Path, nargs="+", required=True)
    merge.add_argument("--output-r1", type=Path, required=True)
    merge.add_argument("--output-r2", type=Path, required=True)
    merge.set_defaults(func=run_merge)

    truncate = subparsers.add_parser("truncate", help="Truncate one FASTQ to a fixed read length")
    truncate.add_argument("--input", type=Path, required=True)
    truncate.add_argument("--output", type=Path, required=True)
    truncate.add_argument("--length", type=int, default=26)
    truncate.set_defaults(func=run_truncate)

    reverse_barcode = subparsers.add_parser(
        "reverse-barcode", help="Reverse-complement the first 16 nt of each read"
    )
    reverse_barcode.add_argument("--input", type=Path, required=True)
    reverse_barcode.add_argument("--output", type=Path, required=True)
    reverse_barcode.set_defaults(func=run_reverse_barcode)

    inspect = subparsers.add_parser("inspect", help="Report read counts and length distributions")
    inspect.add_argument("fastq", type=Path, nargs="+")
    inspect.set_defaults(func=run_inspect)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
