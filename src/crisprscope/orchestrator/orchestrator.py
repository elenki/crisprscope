"""
CRISPRScope Orchestrator: Main Controller

This module contains the main orchestrator class that manages the entire
workflow from raw FASTQ files to CRISPResso2 outputs.
"""
import logging
from pathlib import Path
import sys
import yaml
from importlib import resources

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
    def __init__(self, settings_path: str, config_path: str = None):
        """
        Initializes the orchestrator, creating the main output directory.
        """
        self.settings_path = Path(settings_path)
        self.config_path = config_path
        if not self.settings_path.is_file():
            raise FileNotFoundError(f"Settings file not found at: {self.settings_path}")
        
        # Define the main output directory upon instantiation.
        output_dir_name = self.settings_path.stem
        self.output_dir = self.settings_path.parent / f"{output_dir_name}.orchestrator_output"
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize all attributes the object will use to None or empty.
        self.settings = {}
        self.config = {}
        self.r1_path = None
        self.r2_path = None
        self.barcodes_path = None
        self.demux_dir = self.output_dir / "demultiplexed_fastqs"
        self.known_bc1s = set()
        self.known_bc2s = set()
        self.cell_fastq_paths = []

    def _load_config(self):
        """Loads a YAML config, using the package default if no path is provided."""
        if self.config_path:
            config_file = Path(self.config_path)
            logger.info(f"Loading user-provided config from: {config_file}")
            if not config_file.is_file():
                raise FileNotFoundError(f"Custom config file not found: {self.config_path}")
        else:
            logger.info("No custom config provided. Loading package default.")
            config_file = resources.files('crisprscope').joinpath('config.yaml')

        with open(config_file, 'r') as f:
            self.config = yaml.safe_load(f)

    def load_and_validate_settings(self):
        """Loads settings and config, then validates required file paths."""
        logger.info("Loading and validating orchestrator settings and config...")
        self._load_config()
        self.settings = load_settings(self.settings_path)

        # Ensure all required paths from the settings file are loaded.
        self.r1_path = Path(self.settings.get("r1"))
        self.r2_path = Path(self.settings.get("r2"))
        self.barcodes_path = Path(self.settings.get("barcodes"))

        for p in [self.r1_path, self.r2_path, self.barcodes_path]:
            if not p or not p.is_file(): # Added check for None
                raise FileNotFoundError(f"Required data file not found or specified in settings.txt: {p}")
        
        # Load and split barcodes
        with open(self.barcodes_path, 'r') as f:
            for line in f:
                bc = line.strip()
                if len(bc) == 18:
                    self.known_bc1s.add(bc[:9])
                    self.known_bc2s.add(bc[9:])
        logger.info(f"Successfully loaded and split {len(self.known_bc1s)} unique barcode parts.")


    def run_demultiplexing(self):
        """Runs the FASTQ demultiplexing step using parameters from the config."""
        logger.info("--- Starting Step 1: Demultiplexing ---")
        self.demux_dir = self.output_dir / "demultiplexed_fastqs"
        self.demux_dir.mkdir(exist_ok=True, parents=True)

        try:
            read_params = self.config['read_structure']['mission_bio_v2']
            perf_params = self.config.get('performance_parameters', {})
            chunk_size = perf_params.get('demultiplexer_chunk_size', 1_000_000) # Default to 1M
        except KeyError:
            logger.error("❌ 'read_structure: mission_bio_v2' not found in config.yaml.")
            raise
        
        # Pass the new chunk_size parameter to the demultiplexing function
        demultiplex_fastq(
            r1_path=self.r1_path,
            r2_path=self.r2_path,
            known_bc1s=self.known_bc1s,
            known_bc2s=self.known_bc2s,
            output_dir=self.demux_dir,
            barcode_1_start=read_params['barcode_1_start'],
            barcode_1_end=read_params['barcode_1_end'],
            barcode_2_start=read_params['barcode_2_start'],
            barcode_2_end=read_params['barcode_2_end'],
            chunk_size=chunk_size # <-- PASS IT HERE
        )
        
        self.cell_fastq_paths = list(self.demux_dir.glob("*.fastq.gz"))
        logger.info(f"Demultiplexing complete. Created {len(self.cell_fastq_paths)} cell-specific FASTQ files.")


    def run_crispresso(self):
        """Runs CRISPResso2 on the demultiplexed files."""
        logger.info("--- Starting Step 2: Running CRISPResso2 (Placeholder) ---")
        # --- Future Implementation ---
        logger.warning("CRISPResso2 execution is not yet implemented.")

    def run_summarization(self):
        """Summarizes the results from all CRISPResso2 runs."""
        logger.info("--- Starting Step 3: Summarizing Results (Placeholder) ---")
        # --- Future Implementation ---
        logger.warning("Summarization is not yet implemented.")

    def run_pipeline(self):
        """
        Executes the full orchestrator pipeline from start to finish.
        """
        logger.info("=== Starting CRISPRScope Orchestrator Pipeline ===")
        
        # The output directory is now created in __init__, so we just run the steps.
        self.load_and_validate_settings()
        self.run_demultiplexing()
        self.run_crispresso()
        self.run_summarization()
        
        logger.info("=== CRISPRScope Orchestrator Pipeline Finished ===")