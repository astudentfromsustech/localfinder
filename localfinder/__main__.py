# File: __main__.py

import argparse, sys, importlib.metadata, textwrap, argcomplete, shutil   
from argparse import RawTextHelpFormatter, ArgumentDefaultsHelpFormatter  

from localfinder.commands.bin     import main as bin_tracks_main
from localfinder.commands.calc    import main as calc_corr_main
from localfinder.commands.findreg import main as find_regions_main
from localfinder.commands.viz     import main as visualize_main
from localfinder.pipeline         import run_pipeline
from localfinder.commands.simbg   import main as simbg_main  

# ---------------------------
# Pretty help formatter (wide)
# ---------------------------
### <<< NEW
def _wide_formatter():
    """
    Use a wide, raw-text formatter so newlines in epilog/description are preserved
    and option columns align nicely.
    """
    width = min(120, shutil.get_terminal_size((120, 20)).columns)
    class _F(RawTextHelpFormatter, ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super().__init__(prog, max_help_position=32, width=width)
    return _F


def main():
    # Retrieve package version
    try:
        version = importlib.metadata.version("localfinder")
    except importlib.metadata.PackageNotFoundError:
        version = "0.0.0"  # Fallback version

    parser = argparse.ArgumentParser(
        prog="localfinder",
        description=(
            "localfinder – calculate weighted local correlation (HMC) and enrichment significance (ES)\n"
            "between two genomic tracks, discover significantly different regions (SDRs), and visualize results.\n"
            "GitHub: https://github.com/astudentfromsustech/localfinder"
        ),
        formatter_class=_wide_formatter(),                                      
    )
    parser.add_argument('--version', '-V', action='version',
                        version=f'localfinder {version}',
                        help="Show program's version number and exit.")

    # Create subparsers for subcommands
    subparsers = parser.add_subparsers(
        dest="command", required=True, title="Sub-commands", metavar="{bin,calc,findreg,viz,pipeline,simbg}"  
    )

    # -----------------------
    # Subcommand: bin
    # -----------------------
    parser_bin = subparsers.add_parser(
        'bin',
        help='Bin genomic tracks into fixed-size windows (BedGraph output).',
        description='Bin genomic tracks into fixed-size bins and output BedGraph format.',
        epilog=textwrap.dedent('''\
            Usage Examples
            --------------
            1) Example 1:
               localfinder bin \\
                 --input_files track1.bw track2.bw \\
                 --output_dir ./binned_tracks \\
                 --bin_size 200 \\
                 --chrom_sizes mm10.chrom.sizes \\
                 --chroms chr1 chr2
        '''),
        formatter_class=_wide_formatter(),                                      
    )
    # Grouping for readability                                         
    g_io = parser_bin.add_argument_group("I/O")
    g_io.add_argument('--input_files', nargs='+', required=True, metavar='FILE',
                      help='Input files in BigWig/BedGraph/BAM format.')
    g_io.add_argument('--output_dir', required=True, metavar='DIR',
                      help='Output directory for binned data.')
    g_genome = parser_bin.add_argument_group("Genomics")
    g_genome.add_argument('--bin_size', type=int, default=10000, metavar='BP',
                          help='Size of each bin.')
    g_genome.add_argument('--chrom_sizes', type=str, required=True, metavar='FILE',
                          help='Path to the chromosome sizes file.')
    g_genome.add_argument('--chroms', nargs='+', default=['all'], metavar='CHR',
                          help="'all' or specific chromosomes (e.g. chr1 chr2).")
    g_perf = parser_bin.add_argument_group("Performance")
    g_perf.add_argument('--threads', '-t', type=int, metavar='N', default=1,
                        help='Number of worker processes to run in parallel.')
    parser_bin.set_defaults(func=bin_tracks_main)

    # -----------------------
    # Subcommand: calc
    # -----------------------
    parser_calc = subparsers.add_parser(
        'calc',
        help='Compute HMC & ES tracks from two binned BedGraphs.',
        description='Calculate hormonic-mean weighted local correlation (HMC) and enrichment significance (ES) between two BedGraph tracks.',
        epilog=textwrap.dedent('''\
            Usage Examples
            --------------
            1) Example 1:
               localfinder calc \\
                 --track1 track1.bedgraph --track2 track2.bedgraph \\
                 --output_dir ./results \\
                 --method locP_and_ES --FDR \\
                 --binNum_window 11 --step 1 \\
                 --percentile 90 --percentile_mode all \\
                 --binNum_peak 3 --FC_thresh 1.5 \\
                 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2 \\
                 --threads 4

            2) Example 2:
               localfinder calc \\
                 --track1 track1.bedgraph --track2 track2.bedgraph \\
                 --output_dir ./results \\
                 --percentile 99 --binNum_peak 2 \\
                 --chrom_sizes hg19.chrom.sizes
        '''),
        formatter_class=_wide_formatter(),                                      
    )
    g_io = parser_calc.add_argument_group("I/O")
    g_io.add_argument('--track1', required=True, metavar='BEDGRAPH',
                      help='First input BedGraph file.')
    g_io.add_argument('--track2', required=True, metavar='BEDGRAPH',
                      help='Second input BedGraph file.')
    g_io.add_argument('--output_dir', required=True, metavar='DIR',
                      help='Output directory for results.')

    g_method = parser_calc.add_argument_group("Method")
    g_method.add_argument('--method', choices=['locP_and_ES', 'locS_and_ES'],
                          default='locP_and_ES',
                          help='P=Pearson, S=Spearman for local correlation.')
    g_method.add_argument('--FDR', action='store_true',
                          help='Use Benjamini–Hochberg q-values instead of raw P-values.')
    g_method.add_argument('--binNum_window', type=int, default=11,
                          help='Number of bins in the sliding window.')
    g_method.add_argument('--step', type=int, default=1,
                          help='Step size for the sliding window.')
    g_method.add_argument('--percentile', type=float, default=90,
                          help='Percentile for floor correction of low-coverage bins.')
    g_method.add_argument('--percentile_mode', choices=['all', 'nonzero'], default='all',
                          help="Use all bins or only non-zero bins for percentile.")
    g_method.add_argument('--binNum_peak', type=int, default=3,
                          help='Number of bins in the ES peak window.')
    g_method.add_argument('--FC_thresh', type=float, default=1.5,
                          help='Fold-change threshold used as log base in enrichment.')
    g_method.add_argument('--norm_method', choices=['scale'], default='scale',
                          help='Normalisation strategy.')
    g_method.add_argument('--HMC_scale_pct', type=float, default=0.9995,
                          help='Quantile used to clip+rescale HMC into [0,1].')

    g_genome = parser_calc.add_argument_group("Genomics")
    g_genome.add_argument('--chrom_sizes', type=str, required=True, metavar='FILE',
                          help='Path to the chromosome sizes file.')
    g_genome.add_argument('--chroms', nargs='+', default=['all'], metavar='CHR',
                          help="Chromosomes to process (e.g. chr1 chr2) or 'all'.")

    g_perf = parser_calc.add_argument_group("Performance")
    g_perf.add_argument('--threads', '-t', type=int, metavar='N', default=1,
                        help='Worker processes per chromosome.')
    parser_calc.set_defaults(func=calc_corr_main)

    # -----------------------
    # Subcommand: findreg
    # -----------------------
    parser_find = subparsers.add_parser(
        "findreg",
        help="Find SDRs by ES and further classify SDRs by their HMC",
        description="Find SDRs by ES and further classify SDRs by their HMC",
        epilog=textwrap.dedent("""\
            Usage Examples
            --------------
            1) Example 1:
               localfinder findreg \\
                 --track_E track_ES.bedgraph --track_C track_HMC.bedgraph \\
                 --output_dir ./findreg_out \\
                 --p_thresh 0.05 --binNum_thresh 2 --max_gap_bins 0 \\
                 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

            2) Example 2:
               localfinder findreg \\
                 --track_E track_ES.bedgraph --track_C track_HMC.bedgraph \\
                 --output_dir ./findreg_out \\
                 --chrom_sizes hg19.chrom.sizes
        """),
        formatter_class=_wide_formatter(),                                       
    )
    g_io = parser_find.add_argument_group("I/O")
    g_io.add_argument("--track_E", required=True, metavar='BEDGRAPH', help="ES track (BedGraph).")
    g_io.add_argument("--track_C", required=True, metavar='BEDGRAPH', help="HMC track (BedGraph).")
    g_io.add_argument("--output_dir", required=True, metavar='DIR')

    g_thresh = parser_find.add_argument_group("Thresholds")
    g_thresh.add_argument("--p_thresh", type=float, default=0.05,   help="Bin-level P-value threshold.")
    g_thresh.add_argument("--binNum_thresh", type=int, default=2,   help="Minimum consecutive significant bins per region.")
    g_thresh.add_argument("--max_gap_bins", type=int, default=0,    
                          help="Merge two significant runs if the intervening gap ≤ this number of bins.")

    g_genome = parser_find.add_argument_group("Genomics")
    g_genome.add_argument("--chroms", nargs="+", default=["all"], metavar='CHR')
    g_genome.add_argument("--chrom_sizes", required=True, metavar='FILE')
    parser_find.set_defaults(func=find_regions_main)

    # -----------------------
    # Subcommand: viz
    # -----------------------
    parser_visualize = subparsers.add_parser(
        'viz',
        help='Visualize genomic tracks (pyGenomeTracks/plotly).',
        description='Visualize genomic tracks.',
        epilog=textwrap.dedent('''\
            Usage Examples
            --------------
            1) Example 1:
               localfinder viz \\
                 --input_files track1.bedgraph track2.bedgraph \\
                 --output_file output.html \\
                 --method plotly \\
                 --region chr1 1000000 2000000 \\
                 --colors blue red

            2) Example 2:
               localfinder viz \\
                 --input_files track1.bedgraph track2.bedgraph \\
                 --output_file output.png \\
                 --method pyGenomeTracks \\
                 --region chr1 1000000 2000000
        '''),
        formatter_class=_wide_formatter(),                                       
    )
    g_io = parser_visualize.add_argument_group("I/O")
    g_io.add_argument('--input_files', nargs='+', required=True, metavar='BEDGRAPH',
                      help='Input BedGraph files to visualize.')
    g_io.add_argument('--output_file', required=True, metavar='FILE',
                      help='Output visualization file (e.g. PNG, HTML).')

    g_opts = parser_visualize.add_argument_group("Options")
    g_opts.add_argument('--method', choices=['pyGenomeTracks', 'plotly'], required=True,
                        help='Visualization method to use.')
    g_opts.add_argument('--region', nargs=3, metavar=('CHROM', 'START', 'END'),
                        help='Region to visualize, e.g. chr20 1000000 2000000.')
    g_opts.add_argument('--colors', nargs='+',
                        help='Colors for the tracks (optional).')
    parser_visualize.set_defaults(func=visualize_main)

    # -----------------------
    # Subcommand: pipeline
    # -----------------------
    parser_pipe = subparsers.add_parser(
        'pipeline',
        help='Run bin → calc → findreg in one command.',
        description='Run localfinder sequentially: bin → calc → findreg.',
        epilog=textwrap.dedent('''\
            Usage Examples
            --------------
            1) Example 1:
               localfinder pipeline \\
                 --input_files track1.bedgraph track2.bedgraph \\
                 --output_dir ./results \\
                 --chrom_sizes hg19.chrom.sizes \\
                 --bin_size 200 \\
                 --method locP_and_ES --FDR \\
                 --binNum_window 11 --binNum_peak 3 --step 1 \\
                 --percentile 90 --percentile_mode all \\
                 --FC_thresh 1.5 --norm_method rpkm --HMC_scale_pct 0.9995 \\
                 --p_thresh 0.05 --binNum_thresh 2 --max_gap_bins 0 \\
                 --chroms chr1 chr2 --threads 4

            2) Example 2:
               localfinder pipeline \\
                 --input_files track1.bigwig track2.bigwig \\
                 --output_dir ./results \\
                 --chrom_sizes hg19.chrom.sizes
        '''),
        formatter_class=_wide_formatter(),                                       
    )

    g_io = parser_pipe.add_argument_group("I/O")
    g_io.add_argument("--input_files", nargs="+", required=True, metavar='FILE',
                      help="Input files in BigWig/BedGraph format.")
    g_io.add_argument("--output_dir", required=True, metavar='DIR',
                      help="Output directory for the pipeline results.")

    g_genome = parser_pipe.add_argument_group("Genomics")
    g_genome.add_argument("--chrom_sizes", required=True, metavar='FILE',
                          help="Path to the chromosome sizes file.")
    g_genome.add_argument("--bin_size", type=int, default=200, metavar='BP',
                          help="Size of each bin.")
    g_genome.add_argument('--chroms', nargs='+', default=['all'], metavar='CHR',
                          help='Chromosomes to process or "all".')

    g_calc = parser_pipe.add_argument_group("calc options")
    g_calc.add_argument("--method", choices=["locP_and_ES", "locS_and_ES"],
                        default="locP_and_ES",
                        help='P=Pearson, S=Spearman for local correlation.')
    g_calc.add_argument('--FDR', action='store_true',
                        help='Use Benjamini–Hochberg q-values.')
    g_calc.add_argument("--binNum_window", type=int, default=11,
                        help="Number of bins in the sliding window.")
    g_calc.add_argument("--binNum_peak", type=int, default=3,
                        help="Number of bins of the peak for ES.")
    g_calc.add_argument("--step", type=int, default=1,
                        help="Step size for the sliding window.")
    g_calc.add_argument("--percentile", type=float, default=90,
                        help="Percentile for floor correction of low-coverage bins.")
    g_calc.add_argument('--percentile_mode', choices=['all', 'nonzero'], default='all',
                        help='Use all bins or only non-zero bins for percentile.')
    g_calc.add_argument("--FC_thresh", type=float, default=1.5,
                        help="Fold-change threshold used as log base in enrichment.")
    g_calc.add_argument('--norm_method', choices=['scale']
                        help='Normalisation strategy.')
    g_calc.add_argument("--HMC_scale_pct", type=float, default=0.9995,
                        help="Quantile used to clip+rescale HMC into [0,1].")

    g_find = parser_pipe.add_argument_group("findreg options")
    g_find.add_argument("--p_thresh", type=float, default=0.05,
                        help="P-value threshold for merging significant bins into regions.")
    g_find.add_argument("--binNum_thresh", type=int, default=2,
                        help="Minimum consecutive significant bins per region.")
    g_find.add_argument("--max_gap_bins", type=int, default=0,
                        help="Merge two significant runs if the intervening gap ≤ this number of bins.")

    g_perf = parser_pipe.add_argument_group("Performance")
    g_perf.add_argument('--threads', '-t', type=int, default=1, metavar='N',
                        help='Number of worker processes for bin & calc.')
    parser_pipe.set_defaults(func=run_pipeline)



    # -----------------------
    # Subcommand: simbg
    # -----------------------
    parser_simbg = subparsers.add_parser(                          
        'simbg',
        help='Simulate a multi-scale local background for a binned BedGraph.',
        description='Build a pseudo background BedGraph by applying multi-scale centered rolling windows.',
        epilog=textwrap.dedent('''\
            Usage Examples
            --------------
            1) Infer chrom list from chrom.sizes and run in parallel:
               localfinder simbg \\
                 --input_bedgraph ./track.binSize200.bedgraph \\
                 --output_dir ./bg_out \\
                 --binNum_window 11 \\
                 --chrom_sizes hg19.chrom.sizes \\
                 --chroms all \\
                 --threads 4

            2) Process specific chromosomes:
               localfinder simbg \\
                 --input_bedgraph ./track.binSize200.bedgraph \\
                 --output_dir ./bg_out \\
                 --binNum_window 11 \\
                 --chroms chr1 chr2
        '''),
        formatter_class=_wide_formatter(),
    )
    g_io = parser_simbg.add_argument_group("I/O")                   
    g_io.add_argument('--input_bedgraph', required=True, metavar='BEDGRAPH',
                      help='Input BedGraph (already binned).')
    g_io.add_argument('--output_dir', required=True, metavar='DIR',
                      help='Directory for background outputs.')
    g_opts = parser_simbg.add_argument_group("Options")
    g_opts.add_argument('--binNum_window', type=int, default=11,
                        help='Odd number W; windows used are [W, (W-1)*10+1, (W-1)*100+1].')
    g_opts.add_argument('--bg_method', choices=['poisson', 'zinb'], default='poisson',
                        help='Window estimator (must be supported by utils).')
    g_genome = parser_simbg.add_argument_group("Genomics")
    g_genome.add_argument('--chroms', nargs='+', default=['all'], metavar='CHR',
                          help="Chromosomes to process (e.g. chr1 chr2) or 'all'.")
    g_genome.add_argument('--chrom_sizes', metavar='FILE',
                          help='Chrom sizes; used to expand "all" into explicit chrom list.')
    g_perf = parser_simbg.add_argument_group("Performance")
    g_perf.add_argument('--threads', '-t', type=int, default=1, metavar='N',
                        help='Number of worker processes (per chromosome).')
    parser_simbg.set_defaults(func=simbg_main)



    # Enable auto-completion
    argcomplete.autocomplete(parser)

    # Parse the arguments
    args = parser.parse_args()

    # Execute the appropriate function based on the subcommand
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
