"""
CRISPRScope: Data loading functions for the Integrator module.

This module contains functions to parse and load the various output and
configuration files generated by the CRISPRScope orchestrator and CRISPResso2.
"""

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
import gzip
import os
import traceback
from collections import Counter
from multiprocessing import Pool, cpu_count
import logging
import tempfile

logger = logging.getLogger(__name__)

# --- Helper function for parallel processing ---

def _parse_and_write_parquet(input_fastq_path: Path, output_parquet_path: Path) -> Optional[str]:
    """
    Parses a single CRISPResso_output.fastq.gz file, converts the data to a
    pandas DataFrame, and writes it to a temporary Parquet file.

    This function is designed to be called by a multiprocessing worker and now
    accepts unpacked arguments directly from starmap.

    Args:
        input_fastq_path: Path to the gzipped FASTQ file to parse.
        output_parquet_path: Path to write the output Parquet file to.

    Returns:
        The path to the created Parquet file as a string, or None on failure.
    """
    # no more unpacking needed, arguments are passed directly by starmap
    amplicon_name = input_fastq_path.parent.name.replace("CRISPResso_on_", "", 1)
    
    allele_counts = Counter()
    
    try:
        with gzip.open(input_fastq_path, 'rt') as f:
            for header in f:
                sequence = next(f).strip()
                next(f) # Skip '+' line
                next(f) # Skip quality line

                cell_barcode = header.strip().split(':')[-2]
                allele_counts[(cell_barcode, sequence)] += 1
                
    except StopIteration:
        logger.warning(f"File {input_fastq_path} appears to be truncated or empty.")
        return None
    except Exception as e:
        logger.error(f"Error parsing file {input_fastq_path}: {e}\n{traceback.format_exc()}")
        return None

    if not allele_counts:
        return None

    results = [
        {
            'cell_barcode': cell_barcode,
            'amplicon_name': amplicon_name,
            'allele_sequence': sequence,
            'count': count
        }
        for (cell_barcode, sequence), count in allele_counts.items()
    ]
    
    try:
        df = pd.DataFrame(results)
        df['cell_barcode'] = df['cell_barcode'].astype('category')
        df['amplicon_name'] = df['amplicon_name'].astype('category')
        
        table = pa.Table.from_pandas(df, preserve_index=False)
        pq.write_table(table, output_parquet_path, compression='snappy')
        
        return str(output_parquet_path)
    except Exception as e:
        logger.error(f"Failed to write Parquet file {output_parquet_path}: {e}")
        return None

# --- Main Loader Functions ---

def load_settings(settings_path: Path) -> Dict[str, Any]:
    """Parses the settings.txt configuration file."""
    if not settings_path.is_file():
        logger.error(f"Settings file not found at: {settings_path}")
        raise FileNotFoundError(f"Settings file not found: {settings_path}")
    logger.info(f"Parsing settings from: {settings_path}")
    settings = {}
    with open(settings_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                key, value = parts
                settings[key] = value
            else:
                logger.warning(f"Could not parse line in settings file: '{line}'")
    logger.info(f"Successfully loaded {len(settings)} parameters from settings file.")
    return settings

def load_amplicons(amplicons_path: Path) -> pd.DataFrame:
    """Parses the amplicons.txt file into a DataFrame."""
    if not amplicons_path.is_file():
        logger.error(f"Amplicons file not found at: {amplicons_path}")
        raise FileNotFoundError(f"Amplicons file not found: {amplicons_path}")
    logger.info(f"Parsing amplicons from: {amplicons_path}")
    try:
        df = pd.read_csv(amplicons_path, sep='\t', header=None, index_col=0, names=['sequence', 'guide'])
        df.index.name = 'amplicon_name'
        logger.info(f"Successfully loaded {len(df)} amplicons.")
        return df
    except Exception as e:
        logger.error(f"Failed to parse amplicons file: {e}")
        raise

def load_editing_summary(summary_path: Path) -> pd.DataFrame:
    """Loads the main data matrix from settings.txt.filteredEditingSummary.txt."""
    if not summary_path.is_file():
        logger.error(f"Editing summary file not found at: {summary_path}")
        raise FileNotFoundError(f"Editing summary file not found: {summary_path}")
    logger.info(f"Parsing editing summary from: {summary_path}")
    try:
        df = pd.read_csv(summary_path, sep='\t', index_col=0)
        df.index.name = 'cell_barcode'
        logger.info(f"Successfully loaded editing data for {df.shape[0]} cells and {df.shape[1]} metrics.")
        return df
    except Exception as e:
        logger.error(f"Failed to parse editing summary file: {e}")
        raise

def load_quality_scores(scores_path: Path) -> pd.DataFrame:
    """Loads the cell quality metrics from settings.txt.amplicon_score.txt."""
    if not scores_path.is_file():
        logger.error(f"Amplicon scores file not found at: {scores_path}")
        raise FileNotFoundError(f"Amplicon scores file not found: {scores_path}")
    logger.info(f"Parsing quality scores from: {scores_path}")
    try:
        df = pd.read_csv(scores_path, sep='\t', index_col=0)
        df.index.name = 'cell_barcode'
        logger.info(f"Successfully loaded quality scores for {len(df)} cells.")
        return df
    except Exception as e:
        logger.error(f"Failed to parse quality scores file: {e}")
        raise

def load_crispresso_alleles(crispresso_dir: Path, n_processes: int = None) -> List[Path]:
    """
    Parses allele data from CRISPResso2 FASTQ outputs in parallel, writing
    intermediate results to temporary Parquet files to manage memory usage.

    Args:
        crispresso_dir: Path to the 'settings.txt.crispresso.filtered' directory.
        n_processes: Number of parallel processes to use. Defaults to all available CPUs.

    Returns:
        A list of Path objects pointing to the temporary Parquet files created.
    """
    if not crispresso_dir.is_dir():
        logger.error(f"CRISPResso output directory not found at: {crispresso_dir}")
        raise FileNotFoundError(f"Directory not found: {crispresso_dir}")

    logger.info(f"Scanning for allele data in: {crispresso_dir}")
    
    input_files_to_process = list(crispresso_dir.glob("CRISPResso_on_*/CRISPResso_output.fastq.gz"))
    
    if not input_files_to_process:
        logger.warning("No CRISPResso FASTQ files found to parse. Returning empty list.")
        return []

    logger.info(f"Found {len(input_files_to_process)} amplicon files to parse.")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        logger.info(f"Using temporary directory for intermediate files: {temp_path}")

        job_args = [
            (input_path, temp_path / f"{input_path.parent.name}.parquet")
            for input_path in input_files_to_process
        ]

        if n_processes is None:
            n_processes = cpu_count()
        logger.info(f"Starting parallel parsing with {n_processes} processes...")
        
        with Pool(processes=n_processes) as pool:
            results = pool.starmap(_parse_and_write_parquet, job_args)
            temp_parquet_paths = [Path(p) for p in results if p is not None]
        
        logger.info(f"Parallel processing complete. Generated {len(temp_parquet_paths)} intermediate files.")
        
        if not temp_parquet_paths:
            logger.warning("Parsing yielded no allele data. Returning empty list.")
            return []
            
        # IMPORTANT: The temporary directory is deleted when this 'with' block exits.
        # We must move the successful files to a persistent location before that.
        persistent_temp_dir = Path(tempfile.mkdtemp(prefix="crisprscope_alleles_persistent_"))
        final_paths = []
        for temp_file in temp_parquet_paths:
            final_path = persistent_temp_dir / temp_file.name
            temp_file.rename(final_path)
            final_paths.append(final_path)
            
        logger.info(f"Moved {len(final_paths)} intermediate files to persistent temp directory: {persistent_temp_dir}")

    return final_paths