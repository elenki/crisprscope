"""
Tests for the CRISPRScope Orchestrator module.
"""
import pytest
import gzip
from pathlib import Path

from crisprscope.orchestrator.orchestrator import CRISPRScopeOrchestrator
# Use the same test data definitions from the integrator test for consistency
from .test_integrator import N_CELLS, AMPLICON_SPECS, EDIT_SCENARIOS

def test_demultiplexing_logic(simulation_factory):
    """
    Tests the demultiplexing logic of the orchestrator from end-to-end.
    It verifies that raw FASTQs are correctly split into per-cell files.
    """
    # 1. SETUP: Generate a simulation directory WITH raw FASTQ files
    sim_path = simulation_factory(N_CELLS, AMPLICON_SPECS, EDIT_SCENARIOS, generate_raw_fastq=True)
    
    settings_path = sim_path / "settings.txt"
    config_path = sim_path / "config.yaml"

    # 2. EXECUTION: Run the demultiplexing part of the orchestrator
    orchestrator = CRISPRScopeOrchestrator(settings_path=settings_path, config_path=config_path)
    orchestrator.load_and_validate_settings()
    orchestrator.run_demultiplexing()

    # 3. ASSERTION: Verify the output files and their contents
    
    # Check that the output directory and correct number of files were created
    demux_dir = orchestrator.demux_dir
    assert demux_dir.is_dir()
    output_files = list(demux_dir.glob("*.fastq.gz"))
    # In our scenario, all 3 cells have reads
    assert len(output_files) == N_CELLS

    # --- Critical Assertion: Check the content of one cell's file ---
    
    # Get the barcode for cell_id 0 from the simulator's internal state
    # (A more advanced simulator might return this info directly)
    with open(sim_path / "barcodes.txt") as f:
        barcodes = [line.strip() for line in f]
    cell_0_barcode = barcodes[0]

    cell_0_fastq = demux_dir / f"{cell_0_barcode}.fastq.gz"
    assert cell_0_fastq.is_file()

    # Get the ground truth R2 sequences for cell 0 from our scenario
    ground_truth_seqs = []
    for scenario in EDIT_SCENARIOS:
        if scenario['cell_id'] == 0:
            for seq, count in scenario['alleles'].items():
                ground_truth_seqs.extend([seq] * count)
    
    # Read the sequences from the generated file
    output_seqs = []
    with gzip.open(cell_0_fastq, "rt") as f:
        for header in f:
            seq = next(f).strip()
            next(f) # skip +
            next(f) # skip qual
            output_seqs.append(seq)

    # Assert that the sets of sequences are identical
    assert len(output_seqs) == len(ground_truth_seqs)
    assert set(output_seqs) == set(ground_truth_seqs)

def test_demultiplexing_with_chunking(simulation_factory):
    """
    Tests that the chunk-and-append logic correctly processes all reads
    without loss when the data spans multiple chunks.
    """
    # 1. SETUP: Simulate a dataset with 101 total reads for one cell.
    # We use a prime number to ensure reads are split unevenly across chunks.
    n_reads_for_cell_0 = 101
    sim_path = simulation_factory(
        n_cells=1,
        amplicon_specs={'trp53': 'GATTACA'},
        edit_scenarios=[{
            'cell_id': 0, 
            'amp_name': 'trp53', 
            'alleles': {'GATTACA': n_reads_for_cell_0}
        }],
        generate_raw_fastq=True
    )
    
    settings_path = sim_path / "settings.txt"
    config_path = sim_path / "config.yaml"

    # 2. EXECUTION: Run the orchestrator, but override the chunk_size
    # to be very small, forcing multiple chunking iterations.
    orchestrator = CRISPRScopeOrchestrator(settings_path=settings_path, config_path=config_path)
    orchestrator.load_and_validate_settings()
    
    # Manually override the chunk size for this specific test
    orchestrator.config['performance_parameters']['demultiplexer_chunk_size'] = 10
    
    orchestrator.run_demultiplexing()

    # 3. ASSERTION: Verify the final output file contains all reads.
    demux_dir = orchestrator.demux_dir
    output_files = list(demux_dir.glob("*.fastq.gz"))
    assert len(output_files) == 1

    # Read back the generated file and count the records
    reads_in_output = 0
    with gzip.open(output_files[0], "rt") as f:
        for _ in f:
            reads_in_output += 1
    
    # A fastq record has 4 lines.
    assert (reads_in_output / 4) == n_reads_for_cell_0