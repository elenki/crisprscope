"""
CRISPRScope Orchestrator: FASTQ Demultiplexer

This module contains functions to parse raw, multiplexed FASTQ files
from Mission Bio single-cell experiments and split them into per-cell files.
"""
import gzip
import logging
from pathlib import Path
from collections import defaultdict

logger = logging.getLogger(__name__)

# --- CORRECTED BARCODE POSITIONS ---
# Based on exploratory debugging, we discovered the barcode locations.
# Assumed structure: {26bp spacer}{BC1}{15bp const}{BC2}
BARCODE_1_START = 26
BARCODE_1_END = 35
BARCODE_2_START = 50
BARCODE_2_END = 59

def read_fastq_records(file_handle):
    """Generator to yield records (header, sequence, quality) from a FASTQ file."""
    while True:
        header = file_handle.readline().strip()
        if not header:
            break
        sequence = file_handle.readline().strip()
        plus = file_handle.readline().strip()
        quality = file_handle.readline().strip()
        yield header, sequence, quality

def demultiplex_fastq(r1_path: Path, r2_path: Path, barcodes_path: Path, output_dir: Path):
    """
    Parses raw R1/R2 FASTQ files, extracts cell barcodes, and writes
    the R2 reads to cell-specific FASTQ files.

    Args:
        r1_path: Path to the R1 FASTQ file (contains barcodes).
        r2_path: Path to the R2 FASTQ file (contains biological sequence).
        barcodes_path: Path to the text file of known cell barcodes.
        output_dir: The directory where per-cell FASTQ files will be written.
    """
    logger.info(f"Starting demultiplexing of R1: {r1_path.name} and R2: {r2_path.name}")
    
    # Load the set of known barcodes for fast lookup
    with open(barcodes_path, 'r') as f:
        known_barcodes = {line.strip() for line in f}
    logger.info(f"Successfully loaded {len(known_barcodes)} known cell barcodes.")

    # Dictionary to hold open file handles for writing
    output_files = {}
    
    total_reads = 0
    assigned_reads = 0

    try:
        with gzip.open(r1_path, 'rt') as r1_file, gzip.open(r2_path, 'rt') as r2_file:
            r1_records = read_fastq_records(r1_file)
            r2_records = read_fastq_records(r2_file)

            for r1_header, r1_seq, r1_qual in r1_records:
                total_reads += 1
                r2_header, r2_seq, r2_qual = next(r2_records)

                # Extract the two 9bp barcodes from the R1 read sequence
                bc1 = r1_seq[BARCODE_1_START:BARCODE_1_END]
                bc2 = r1_seq[BARCODE_2_START:BARCODE_2_END]
                
                # Check if both parts are in the known list
                if bc1 in known_barcodes and bc2 in known_barcodes:
                    full_barcode = bc1 + bc2
                    assigned_reads += 1
                    
                    # If this is the first time we've seen this barcode, open a new file
                    if full_barcode not in output_files:
                        out_path = output_dir / f"{full_barcode}.fastq.gz"
                        output_files[full_barcode] = gzip.open(out_path, 'wt')
                    
                    # Write the R2 record to the correct file
                    output_files[full_barcode].write(f"{r2_header}\n{r2_seq}\n+\n{r2_qual}\n")

                if total_reads % 1000000 == 0:
                    logger.info(f"   ...processed {total_reads:,} reads.")

    finally:
        # Crucially, close all the file handles we opened
        for handle in output_files.values():
            handle.close()

    percent_assigned = (assigned_reads / total_reads * 100) if total_reads > 0 else 0
    logger.info("Demultiplexing complete.")
    logger.info(f"   Total reads processed: {total_reads:,}")
    logger.info(f"   Reads assigned to known barcodes: {assigned_reads:,} ({percent_assigned:.2f}%)")
    logger.info(f"   Created {len(output_files)} cell-specific FASTQ files in {output_dir}")