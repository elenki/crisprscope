"""
CRISPRScope Integrator: Performance and Scalability Benchmark

This script measures the runtime and peak memory usage of the `build_anndata`
function across a range of dataset sizes (number of cells).
"""

import sys
import time
import tracemalloc
import argparse
import logging
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import tempfile

# Add src and tests directories to path
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root / 'src'))
sys.path.insert(0, str(project_root / 'tests'))

from crisprscope.api import build_anndata
from conftest import CRISPRScopeSimulator

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Benchmark Configuration ---
# Use a representative number of amplicons from the paper
N_AMPLICONS = 40
# Define a consistent, simple edit scenario for each cell
EDIT_COMPLEXITY_PER_CELL = 5 # Number of amplicons to have edits in each cell

def run_benchmark(n_cells: int, simulation_path: Path) -> dict:
    """
    Runs a single benchmark for a given number of cells.

    Args:
        n_cells: The number of cells to simulate.
        simulation_path: The temporary directory to build the simulation in.

    Returns:
        A dictionary containing the benchmark results.
    """
    logger.info(f"--- Running benchmark for {n_cells:,} cells ---")
    
    # 1. Setup: Simulate the dataset
    logger.info("Generating simulated dataset...")
    amplicon_specs = {f'amp_{i}': 'A'*200 for i in range(N_AMPLICONS)}
    
    # Create a semi-realistic but consistent set of edits
    edit_scenarios = []
    for i in range(n_cells):
        for j in range(EDIT_COMPLEXITY_PER_CELL):
            amp_name = f'amp_{(i + j) % N_AMPLICONS}'
            edit_scenarios.append({
                'cell_id': i, 'amp_name': amp_name,
                'alleles': {'A'*199 + '-': 20, 'A'*200: 10} # 2:1 het edit
            })

    simulator = CRISPRScopeSimulator(output_dir=simulation_path)
    simulator.generate(n_cells, amplicon_specs, edit_scenarios)
    logger.info("Simulation complete.")

    # 2. Execution & Profiling
    tracemalloc.start()
    start_time = time.perf_counter()

    # Run the main build function
    build_anndata(base_path=simulation_path)

    end_time = time.perf_counter()
    _, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    elapsed_time = end_time - start_time
    peak_mem_mb = peak_mem / 1024**2
    
    logger.info(f"Completed in {elapsed_time:.2f}s, Peak Memory: {peak_mem_mb:.2f} MB")
    
    return {'n_cells': n_cells, 'time_s': elapsed_time, 'mem_mb': peak_mem_mb}

def main():
    """Main function to run the benchmark suite."""
    parser = argparse.ArgumentParser(description="CRISPRScope Integrator Benchmark")
    parser.add_argument(
        '--cell-counts',
        nargs='+',
        type=int,
        default=[100, 1000, 5000, 10000],
        help="A list of cell counts to benchmark (e.g., 100 1000 10000)."
    )
    parser.add_argument(
        '--output-plot',
        type=str,
        default="benchmark_results.png",
        help="Path to save the output plot."
    )
    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    results = []
    for n_cells in args.cell_counts:
        with tempfile.TemporaryDirectory() as temp_dir:
            result = run_benchmark(n_cells, Path(temp_dir))
            results.append(result)

    # --- 3. Report & Plot Results ---
    results_df = pd.DataFrame(results).set_index('n_cells')
    print("\n--- Benchmark Results ---")
    print(results_df)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot Time vs. N_Cells
    results_df['time_s'].plot(kind='line', ax=ax1, marker='o')
    ax1.set_title('Runtime Scaling')
    ax1.set_xlabel('Number of Cells')
    ax1.set_ylabel('Time (seconds)')
    ax1.grid(True, linestyle='--', alpha=0.6)

    # Plot Memory vs. N_Cells
    results_df['mem_mb'].plot(kind='line', ax=ax2, marker='o', color='r')
    ax2.set_title('Peak Memory Usage Scaling')
    ax2.set_xlabel('Number of Cells')
    ax2.set_ylabel('Peak Memory (MB)')
    ax2.grid(True, linestyle='--', alpha=0.6)

    fig.suptitle('CRISPRScope Integrator Performance', fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    
    plt.savefig(args.output_plot)
    logger.info(f"Benchmark plot saved to {args.output_plot}")

if __name__ == "__main__":
    main()