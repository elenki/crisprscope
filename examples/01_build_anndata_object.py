"""
CRISPRScope Example 1: Validating Core Data Loaders

This script serves as a validation test for the foundational data loading
functions in the `loaders.py` module. It performs the following steps:

1.  Sets up the Python path to find the `crisprscope` source code.
2.  Defines the path to the experimental data directory.
3.  Calls `load_settings` to parse the main settings.txt file and prints a sample.
4.  Retrieves the path to the amplicons.txt file from the loaded settings.
5.  Calls `load_amplicons` to parse the amplicon data into a DataFrame and
    prints the shape and head of the resulting table.
"""
from pathlib import Path
import sys
import logging
import pandas as pd

# This is a common pattern in example scripts. It adds the 'src' directory
# to Python's path, allowing us to import our `crisprscope` package even
# before it's formally installed via pip.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from crisprscope.integrator.loaders import load_settings, load_amplicons

# Configure logging to display informational messages
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def main():
    """Main function to run the validation."""
    
    # --- Configuration ---
    # Define the absolute path to your data directory.
    # This should be the directory containing 'settings.txt' and other outputs.
    base_data_path = Path("/uufs/chpc.utah.edu/common/home/u6046470/clement/projects/20230818_scCRISPR/analysis/01_run_on_20200804_BaF3_revision")
    settings_file = base_data_path / "settings.txt"

    print("---" * 15)
    print("--- CRISPRScope Loader Validation Script ---")
    print("---" * 15)

    # --- Step 1: Validate load_settings ---
    print("\n[Step 1/2] Validating the settings loader...")
    try:
        settings = load_settings(settings_file)
        print("✅ Success: `load_settings` executed without errors.")
        print(f"   - Found {len(settings)} parameters in '{settings_file.name}'.")
        
        # Print a sample of the loaded settings for manual verification
        print("   - Sample parameters:")
        for i, (key, value) in enumerate(settings.items()):
            if i >= 5:
                break
            print(f"     - {key}: {value}")
        
    except Exception as e:
        print(f"❌ Error loading settings: {e}")
        # Exit if settings can't be loaded, as the next step depends on it.
        return

    # --- Step 2: Validate load_amplicons ---
    print("\n[Step 2/2] Validating the amplicons loader...")
    try:
        # Dynamically get the amplicons file path from the loaded settings
        amplicon_file_path_str = settings.get('amplicons')
        if not amplicon_file_path_str:
            print("❌ Error: 'amplicons' key not found in settings.txt. Cannot proceed.")
            return
        
        amplicon_file = Path(amplicon_file_path_str)
        print(f"   - Found amplicon file path in settings: {amplicon_file}")
        
        amplicons_df = load_amplicons(amplicon_file)
        print("✅ Success: `load_amplicons` executed without errors.")
        print(f"   - Loaded {amplicons_df.shape[0]} amplicons with {amplicons_df.shape[1]} attributes.")
        
        # Print the head of the DataFrame for verification
        print("   - Head of the amplicons DataFrame:")
        # Use pandas to_string() for better console formatting
        print(amplicons_df.head().to_string())

    except Exception as e:
        print(f"❌ Error loading amplicons: {e}")
        return
        
    print("\n---" * 15)
    print("--- Validation Complete ---")
    print("---" * 15)


if __name__ == "__main__":
    main()