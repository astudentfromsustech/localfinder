# File: commands/findreg.py

import os
from localfinder.utils import (                       # util helpers
    find_significantly_different_regions as util_find,  
    get_chromosomes_from_chrom_sizes
)

# ----------------------------------------------------------------------
# CLI entry-point (wired in __main__.py)
# ----------------------------------------------------------------------
def main(args):                                        
    # • pull CLI arguments
    track_ES_file   = args.track_E                     
    track_HMC_file = args.track_C                     
    output_dir      = args.output_dir
    p_thresh        = args.p_thresh                    
    binNum_thresh   = args.binNum_thresh               
    max_gap_bins    = args.max_gap_bins
    chroms          = args.chroms
    chrom_sizes     = args.chrom_sizes

    os.makedirs(output_dir, exist_ok=True)

    # • expand “all” → actual chromosome list
    if chroms == ['all'] or chroms is None:
        chroms = get_chromosomes_from_chrom_sizes(chrom_sizes)
        print(f"'chroms' set to all chromosomes: {chroms}")
    else:
        print(f"'chroms' set to specified chromosomes: {chroms}")

    # • call the utility function (aliased as util_find)
    util_find(
        track_ES_file   = track_ES_file,               
        track_HMC_file = track_HMC_file,            
        output_dir      = output_dir,
        p_thresh        = p_thresh,                    
        binNum_thresh   = binNum_thresh,               
        max_gap_bins    = max_gap_bins,
        chroms          = chroms,
        chrom_sizes     = chrom_sizes
    )