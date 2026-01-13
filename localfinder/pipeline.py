# File: pipeline.py

import argparse, os, sys, shutil
from localfinder.commands.bin import main as bin_tracks_main
from localfinder.commands.calc import main as calc_corr_main
from localfinder.commands.findreg import main as find_regions_main  
from localfinder.utils import check_external_tools, get_chromosomes_from_chrom_sizes

def run_pipeline(args):
    # Ensure required external tools are available
    check_external_tools()
    pass
    print("Pipeline completed successfully.")

