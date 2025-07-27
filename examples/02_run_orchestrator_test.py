"""
CRISPRScope Orchestrator: Test Script

This script tests the initialization and the first step (demultiplexing)
of the CRISPRScopeOrchestrator pipeline.
"""
import logging
from pathlib import Path
import sys

# Add the src directory to the Python path to allow for direct imports
# This is a common pattern for running example scripts
project_root = Path(__file__).resolve().parent.parent
sys.path.append(str(project_root / "src"))

from crisprscope.orchestrator.orchestrator import CRISPRScopeOrchestrator

def main():
    """Main function to run the orchestrator test."""
    
    # Configure logging to see the output from the orchestrator
    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')
    
    logger = logging.getLogger(__name__)

    print("\n-------------------------------------------------")
    print("--- CRISPRScope Orchestrator Test (Demux) ---")
    print("-------------------------------------------------\n")

    # --- Configuration ---
    # Path to the main settings file for the run
    # This should be the only path you need to configure
    base_analysis_path = Path("/uufs/chpc.utah.edu/common/home/u6046470/clement/projects/20230818_scCRISPR/analysis/01_run_on_20200804_BaF3_revision")
    settings_file = base_analysis_path / "settings.txt"

    try:
        logger.info(f"Initializing orchestrator with settings file: {settings_file}")
        
        # 1. Initialize the Orchestrator
        orchestrator = CRISPRScopeOrchestrator(settings_path=settings_file)

        # 2. Execute the pipeline
        # NOTE: For a real run, this will be very time-consuming.
        # For this test, we can interrupt it (Ctrl+C) after we see it start
        # processing reads, or let it run on a small subset of data.
        orchestrator.run_pipeline()

    except FileNotFoundError as e:
        logger.error(f"❌ A required file was not found. Please check your paths.")
        logger.error(e)
    except Exception as e:
        logger.error(f"❌ An unexpected error occurred during the orchestrator run.")
        logger.error(e, exc_info=True) # exc_info=True prints the full traceback

if __name__ == "__main__":
    main()