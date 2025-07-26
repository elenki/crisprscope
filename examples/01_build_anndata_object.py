"""
CRISPRScope Example 1: Building and Saving the Final AnnData Object

This script demonstrates the final, user-facing API for the CRISPRScope
Integrator. It uses a single function call to perform all loading,
building, and saving steps.
"""
from pathlib import Path
import sys
import logging

# Add the 'src' directory to Python's path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

# Import the new high-level API function
from crisprscope.api import build_anndata

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    """Main function to run the final build process."""
    
    # --- Configuration ---
    base_data_path = "/uufs/chpc.utah.edu/common/home/u6046470/clement/projects/20230818_scCRISPR/analysis/01_run_on_20200804_BaF3_revision"
    # Define where we want to save the final output file
    output_file_path = "crisprscope_output.h5ad"

    print("---" * 15)
    print("--- Running CRISPRScope End-to-End Builder ---")
    print("---" * 15)

    try:
        # --- Run the entire pipeline with one function call ---
        adata = build_anndata(
            base_path=base_data_path,
            output_path=output_file_path
        )

        print("\n---" * 15)
        print("--- Pipeline Complete ---")
        print("---" * 15)
        
        print("\n✅ Success: AnnData object built and saved.")
        print(f"   - Final object located at: {output_file_path}")
        print(f"   - Object dimensions: {adata.n_obs} cells × {adata.n_vars} amplicons")
        print(f"   - Available layers: {list(adata.layers.keys())}")

    except Exception as e:
        print(f"\n❌ An error occurred during the build process: {e}")

if __name__ == "__main__":
    main()