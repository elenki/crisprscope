"""
CRISPRScope Orchestrator: FASTQ Demultiplexer

This module contains high-performance functions to parse raw, multiplexed 
FASTQ files and split them into per-cell files using the dnaio library.
"""
import logging
from pathlib import Path
from typing import Set
import dnaio
from collections import defaultdict
import itertools
import gzip

logger = logging.getLogger(__name__)

def demultiplex_fastq(
    r1_path: Path, 
    r2_path: Path, 
    known_bc1s: Set[str], 
    known_bc2s: Set[str], 
    output_dir: Path,
    barcode_1_start: int,
    barcode_1_end: int,
    barcode_2_start: int,
    barcode_2_end: int,
    chunk_size: int = 1_000_000
):
    """
    Parses raw R1/R2 FASTQ files in chunks, extracts cell barcodes, and appends
    the R2 reads to cell-specific FASTQ files. This approach is scalable and
    avoids hitting open file limits.
    """
    logger.info(f"Starting high-performance demultiplexing of R1: {r1_path.name}")
    logger.info(f"Processing in chunks of {chunk_size:,} reads.")
    
    total_reads, assigned_reads = 0, 0
    
    with dnaio.open(r1_path, file2=r2_path, mode='r') as reader:
        while True:
            chunk = list(itertools.islice(reader, chunk_size))
            if not chunk:
                break

            total_reads += len(chunk)
            chunk_records_by_cell = defaultdict(list)

            # 1. Process the in-memory chunk
            for r1, r2 in chunk:
                # dnaio returns strings, so no decoding is needed.
                r1_seq = r1.sequence
                
                bc1 = r1_seq[barcode_1_start:barcode_1_end]
                bc2 = r1_seq[barcode_2_start:barcode_2_end]
                
                if bc1 in known_bc1s and bc2 in known_bc2s:
                    full_barcode = bc1 + bc2
                    # Construct the record string directly from the string attributes
                    record_str = f"@{r2.name}\n{r2.sequence}\n+\n{r2.qualities}\n"
                    chunk_records_by_cell[full_barcode].append(record_str)
                    assigned_reads += 1

            # 2. Append the chunk to the respective files
            logger.info(f"   ...processed {total_reads:,} reads. Appending chunk to files...")
            for barcode, records in chunk_records_by_cell.items():
                output_path = output_dir / f"{barcode}.fastq.gz"
                with gzip.open(output_path, 'at') as f_out:
                    f_out.writelines(records)
    
    output_files = list(output_dir.glob("*.fastq.gz"))
    percent_assigned = (assigned_reads / total_reads * 100) if total_reads > 0 else 0
    logger.info("Demultiplexing complete.")
    logger.info(f"   Total reads processed: {total_reads:,}")
    logger.info(f"   Reads assigned to known barcodes: {assigned_reads:,} ({percent_assigned:.2f}%)")
    logger.info(f"   Created {len(output_files)} cell-specific FASTQ files in {output_dir}")