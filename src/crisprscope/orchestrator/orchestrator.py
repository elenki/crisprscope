"""
CRISPRScope Orchestrator: Main Controller

This module contains the main orchestrator class that manages the entire
workflow from raw FASTQ files to CRISPResso2 outputs.
"""
import logging
from pathlib import Path
import sys

# Add the src directory to the Python path for imports
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from crisprscope.integrator.loaders import load_settings
from .demultiplex import demultiplex_fastq

logger = logging.getLogger(__name__)

class CRISPRScopeOrchestrator:
    """
    Orchestrates the running of the demultiplexing and CRISPResso2 pipeline.
    """
    def __init__(self, settings_path: str):
        """
        Initializes the orchestrator with the path to the settings file.

        Args:
            settings_path: Path to the main settings.txt configuration file.
        """
        self.settings_path = Path(settings_path)
        if not self.settings_path.is_file():
            raise FileNotFoundError(f"Settings file not found at: {self.settings_path}")
        
        self.settings = None
        self.r1_path = None
        self.r2_path = None
        self.barcodes_path = None
        self.output_dir = None
        self.demux_dir = None
        self.known_barcodes = None
        self.cell_fastq_paths = None

    def load_and_validate_settings(self):
        """Loads settings and validates required file paths."""
        logger.info("Loading and validating orchestrator settings...")
        self.settings = load_settings(self.settings_path)

        # Validate required paths from settings file
        self.r1_path = Path(self.settings.get("r1"))
        self.r2_path = Path(self.settings.get("r2"))
        self.barcodes_path = Path(self.settings.get("barcodes"))

        for p in [self.r1_path, self.r2_path, self.barcodes_path]:
            if not p.is_file():
                raise FileNotFoundError(f"Required data file not found: {p}")
        
        # Load known barcodes
        with open(self.barcodes_path, 'r') as f:
            self.known_barcodes = {line.strip() for line in f}
        logger.info(f"Successfully loaded {len(self.known_barcodes)} known cell barcodes.")


    def run_demultiplexing(self):
        """Runs the FASTQ demultiplexing step."""
        logger.info("--- Starting Step 1: Demultiplexing ---")
        self.demux_dir = self.output_dir / "demultiplexed_fastqs"
        self.demux_dir.mkdir(exist_ok=True)

        # Corrected function call: Pass the barcodes_path, not the loaded barcodes
        demultiplex_fastq(
            r1_path=self.r1_path,
            r2_path=self.r2_path,
            barcodes_path=self.barcodes_path, 
            output_dir=self.demux_dir
        )
        
        # Store the paths of the created files for the next step
        self.cell_fastq_paths = list(self.demux_dir.glob("*.fastq.gz"))
        logger.info(f"Demultiplexing complete. Created {len(self.cell_fastq_paths)} cell-specific FASTQ files.")


    def run_crispresso(self):
        """Runs CRISPResso2 on the demultiplexed files."""
        logger.info("--- Starting Step 2: Running CRISPResso2 (Placeholder) ---")
        if not self.cell_fastq_paths:
            logger.warning("No demultiplexed FASTQ files found to run CRISPResso2 on.")
            return
        # --- Future Implementation ---
        # This is where the logic to build and run CRISPResso2 commands in parallel will go.
        # For each file in self.cell_fastq_paths:
        #   - Construct a CRISPResso2 command
        #   - Submit the command to the cluster (e.g., using subprocess or a job scheduler)
        logger.warning("CRISPResso2 execution is not yet implemented.")

    def run_summarization(self):
        """Summarizes the results from all CRISPResso2 runs."""
        logger.info("--- Starting Step 3: Summarizing Results (Placeholder) ---")
        # --- Future Implementation ---
        # This is where the logic to parse all the individual CRISPResso2
        # output folders and create the summary files (e.g., filteredEditingSummary.txt) will go.
        logger.warning("Summarization is not yet implemented.")

    def run_pipeline(self):
        """
        Executes the full orchestrator pipeline from start to finish.
        """
        logger.info("=== Starting CRISPRScope Orchestrator Pipeline ===")
        
        # Define the main output directory, named after the settings file
        output_dir_name = self.settings_path.stem
        self.output_dir = self.settings_path.parent / f"{output_dir_name}.orchestrator_output"
        self.output_dir.mkdir(exist_ok=True)
        
        self.load_and_validate_settings()
        self.run_demultiplexing()
        self.run_crispresso()
        self.run_summarization()
        
        logger.info("=== CRISPRScope Orchestrator Pipeline Finished ===")