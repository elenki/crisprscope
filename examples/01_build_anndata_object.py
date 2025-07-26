"""
CRISPRScope Example 1: Validating Core Data Loaders

This script serves as a validation test for the foundational data loading
functions in the `loaders.py` module.
"""
from pathlib import Path
import sys
import logging
import pandas as pd

# Ensure the src directory is in the Python path for module imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from crisprscope.integrator.loaders import (
    load_settings, 
    load_amplicons,
    load_editing_summary,
    load_quality_scores
)

# Configure logging to display informational messages
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def main():
    """Main function to run the validation."""
    
    # --- Configuration ---
    base_data_path = Path("/uufs/chpc.utah.edu/common/home/u6046470/clement/projects/20230818_scCRISPR/analysis/01_run_on_20200804_BaF3_revision")
    settings_file = base_data_path / "settings.txt"
    summary_file = base_data_path / "settings.txt.filteredEditingSummary.txt"
    scores_file = base_data_path / "settings.txt.amplicon_score.txt"

    print("---" * 15)
    print("--- CRISPRScope Loader Validation Script ---")
    print("---" * 15)

    # --- Step 1: Validate load_settings ---
    print("\n[Step 1/4] Validating the settings loader...")
    try:
        settings = load_settings(settings_file)
        print("✅ Success: `load_settings` executed.")
    except Exception as e:
        print(f"❌ Error loading settings: {e}")
        return

    # --- Step 2: Validate load_amplicons ---
    print("\n[Step 2/4] Validating the amplicons loader...")
    try:
        amplicon_file_path_str = settings.get('amplicons')
        amplicon_file = Path(amplicon_file_path_str)
        amplicons_df = load_amplicons(amplicon_file)
        print("✅ Success: `load_amplicons` executed.")
        print(f"   - Loaded {amplicons_df.shape[0]} amplicons with {amplicons_df.shape[1]} attributes.")
    except Exception as e:
        print(f"❌ Error loading amplicons: {e}")
        return
        
    # --- Step 3: Validate load_editing_summary ---
    print("\n[Step 3/4] Validating the editing summary loader...")
    try:
        summary_df = load_editing_summary(summary_file)
        print("✅ Success: `load_editing_summary` executed.")
        print(f"   - Loaded summary data for {summary_df.shape[0]} cells across {summary_df.shape[1]} metrics.")
        print("   - Head of the summary DataFrame:")
        print(summary_df.head().iloc[:, :5].to_string()) # Show first 5 columns for brevity
    except Exception as e:
        print(f"❌ Error loading editing summary: {e}")
        return

    # --- Step 4: Validate load_quality_scores ---
    print("\n[Step 4/4] Validating the quality scores loader...")
    try:
        scores_df = load_quality_scores(scores_file)
        print("✅ Success: `load_quality_scores` executed.")
        print(f"   - Loaded quality scores for {scores_df.shape[0]} cells.")
        print("   - Head of the scores DataFrame:")
        print(scores_df.head().to_string())
    except Exception as e:
        print(f"❌ Error loading quality scores: {e}")
        return

    print("\n---" * 15)
    print("--- Validation Complete ---")
    print("---" * 15)


if __name__ == "__main__":
    main()