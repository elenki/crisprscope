"""
CRISPRScope Example 1: Building and Validating the Full AnnData Object

This script demonstrates the complete data loading and AnnData construction
pipeline, including the parsing of allele data and the population of all
advanced data layers (.X, counts, zygosity, and alleles).
"""
from pathlib import Path
import sys
import logging
import pandas as pd

# Add the 'src' directory to Python's path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from crisprscope.integrator.loaders import (
    load_settings, 
    load_amplicons,
    load_editing_summary,
    load_quality_scores,
    load_crispresso_alleles
)
from crisprscope.integrator.builder import CRISPRScopeAnnDataBuilder

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def main():
    """Main function to run the validation."""
    
    # --- Configuration ---
    base_data_path = Path("/uufs/chpc.utah.edu/common/home/u6046470/clement/projects/20230818_scCRISPR/analysis/01_run_on_20200804_BaF3_revision")

    print("---" * 15)
    print("--- CRISPRScope AnnData Builder Test (Full) ---")
    print("---" * 15)

    # --- Step 1: Load all data sources ---
    print("\n[Step 1/2] Loading all data files...")
    try:
        settings = load_settings(base_data_path / "settings.txt")
        amplicon_file = Path(settings.get('amplicons'))
        amplicons_df = load_amplicons(amplicon_file)
        summary_df = load_editing_summary(base_data_path / "settings.txt.filteredEditingSummary.txt")
        scores_df = load_quality_scores(base_data_path / "settings.txt.amplicon_score.txt")
        
        # Load the new allele data
        alleles_data = load_crispresso_alleles(base_data_path / "settings.txt.crispresso.filtered")
        
        print("✅ Success: All data loaded.")
    except Exception as e:
        print(f"❌ Error during data loading: {e}")
        return

    # --- Step 2: Build the AnnData object ---
    print("\n[Step 2/2] Building AnnData object from loaded data...")
    try:
        # Initialize the builder, now including the alleles_data
        builder = CRISPRScopeAnnDataBuilder(
            settings=settings,
            amplicons=amplicons_df,
            editing_summary=summary_df,
            quality_scores=scores_df,
            alleles_data=alleles_data
        )
        
        # Build the object
        adata = builder.build()
        
        print("✅ Success: AnnData object built.")
        print(f"   - Object dimensions: {adata.n_obs} cells × {adata.n_vars} amplicons")
        
        print("\n--- Verifying .obs (Cell Metadata) ---")
        print(adata.obs.head().to_string())
        
        print("\n--- Verifying .var (Amplicon Metadata) ---")
        print(adata.var.head().to_string())
        
        print("\n--- Verifying .X (Modification Percentage) ---")
        print("   - Displaying top-left 5x5 slice:")
        print(pd.DataFrame(adata.X[:5, :5], 
                           index=adata.obs_names[:5], 
                           columns=adata.var_names[:5]))

        print("\n--- Verifying .layers['counts'] (Total Counts) ---")
        print("   - Displaying top-left 5x5 slice:")
        print(pd.DataFrame(adata.layers['counts'][:5, :5], 
                           index=adata.obs_names[:5], 
                           columns=adata.var_names[:5]))
                           
        print("\n--- Verifying .layers['zygosity'] (Zygosity Codes) ---")
        print("   - Zygosity Encoding:", adata.uns['zygosity_encoding'])
        print("   - Displaying top-left 5x5 slice:")
        print(pd.DataFrame(adata.layers['zygosity'][:5, :5], 
                           index=adata.obs_names[:5], 
                           columns=adata.var_names[:5]))

        print("\n--- Verifying .layers['alleles'] (Allele Strings) ---")
        print("   - Displaying top-left 5x5 slice:")
        # Increase pandas display width to see full allele strings
        with pd.option_context('display.max_colwidth', 150):
            print(pd.DataFrame(adata.layers['alleles'][:5, :5], 
                               index=adata.obs_names[:5], 
                               columns=adata.var_names[:5]))

    except Exception as e:
        print(f"❌ Error during AnnData construction: {e}")
        return

    print("\n---" * 15)
    print("--- Builder Test Complete ---")
    print("---" * 15)

if __name__ == "__main__":
    main()