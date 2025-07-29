"""
CRISPRScope: Command-Line Interface (CLI)

This module provides the main entry point for running CRISPRScope from the terminal.
"""

import argparse
import logging
from pathlib import Path

from .api import build_anndata

def main():
    """Main function for the CRISPRScope command-line tool."""
    
    # --- Main Parser ---
    parser = argparse.ArgumentParser(
        description="CRISPRScope: A toolkit for analyzing single-cell CRISPR screen data."
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s 0.1.0'  # Consider pulling this from __init__.py later
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    # --- Build Command Parser ---
    parser_build = subparsers.add_parser(
        "build",
        help="Build an AnnData object from CRISPResso2 output.",
        description="This command runs the Integrator module to process a CRISPResso2 output directory and create a single, analysis-ready .h5ad file."
    )
    parser_build.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the input directory (containing settings.txt and subdirectories)."
    )
    parser_build.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path for the output .h5ad file."
    )
    parser_build.add_argument(
        "--config",
        type=str,
        default=None,
        help="Optional path to a custom YAML config file for methodological parameters. If not provided, the package default is used."
    )
    parser_build.add_argument(
        "--verbose",
        action="store_true",
        help="Enable detailed logging output."
    )

    # --- Parse Arguments ---
    args = parser.parse_args()

    # --- Configure Logging ---
    if args.verbose:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(level=logging.WARNING, format='%(levelname)s: %(message)s')

    # --- Execute Command ---
    if args.command == "build":
        try:
            build_anndata(
                base_path=args.input, 
                output_path=args.output,
                config_path=args.config  # Pass the new argument to the API
            )
            print(f"\n✅ Success! AnnData object saved to: {args.output}")
        except Exception as e:
            logging.error(f"An error occurred during the build process: {e}", exc_info=args.verbose)
            # Exit with a non-zero status code to indicate failure
            exit(1)

if __name__ == '__main__':
    main()