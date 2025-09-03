import argparse
import importlib.metadata
import sys
import textwrap

from argparse import RawDescriptionHelpFormatter

from localfinder.commands.bin import main as bin_tracks_main
from localfinder.commands.calc import main as calc_corr_main
from localfinder.commands.findreg import main as find_regions_main
from localfinder.commands.viz import main as visualize_main
from localfinder.pipeline import run_pipeline


def main():
    # Retrieve package version
    try:
        version = importlib.metadata.version("localfinder")
    except importlib.metadata.PackageNotFoundError:
        version = "0.0.0"  # Fallback version

    parser = argparse.ArgumentParser(
        prog="localfinder",
        description=(
            "localfinder – calculate weighted local correlation (HMC) and enrichment "
            "significance (ES) between two genomic tracks, optionally discover "
            "significantly different regions, and visualise results. "
            "GitHub: https://github.com/astudentfromsustech/localfinder"
        ),
        formatter_class=RawDescriptionHelpFormatter,  ### <<< CHANGED
    )
    parser.add_argument(
        "--version",
        "-V",
        action="version",
        version=f"localfinder {version}",
        help="Show program's version number and exit.",
    )

    # Create subparsers for subcommands
    subparsers = parser.add_subparsers(dest="command", required=True, title="Sub-commands")

    # Subcommand: bin (alias: bin_tracks)
    parser_bin = subparsers.add_parser(
        "bin",
        help="Convert input files into bins with BedGraph format.",
        description="Bin genomic tracks into fixed-size bins and output BedGraph format.",
        epilog=textwrap.dedent("""\
            Usage Example 1:
                localfinder bin --input_files track1.bw track2.bw --output_dir ./binned_tracks --bin_size 200 --chrom_sizes mm10.chrom.sizes --chroms chr1 chr2

            Usage Example 2:
                localfinder bin --input_files track1.bigwig track2.bigwig --output_dir ./binned_tracks --bin_size 200 --chrom_sizes hg19.chrom.sizes --chroms all --threads 4
            """),
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser_bin.add_argument(
        "--input_files",
        nargs="+",
        required=True,
        help="Input files in BigWig/BedGraph/BAM format.",
    )
    parser_bin.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for binned data.",
    )
    parser_bin.add_argument(
        "--bin_size",
        type=int,
        default=200,
        help="Size of each bin (default: 200).",
    )
    parser_bin.add_argument(
        "--chrom_sizes",
        type=str,
        required=True,
        help="Path to the chromosome sizes file.",
    )
    parser_bin.add_argument(
        "--chroms",
        nargs="+",
        default=["all"],
        help="'all' or specific chromosomes (e.g. chr1 chr2).",
    )
    parser_bin.add_argument(
        "--threads",
        "-t",
        type=int,
        metavar="N",
        default=1,
        help="Number of worker processes to run in parallel (default: 1).",
    )
    parser_bin.set_defaults(func=bin_tracks_main)

    # Subcommand: calc (alias: calculate_localCorrelation_and_enrichmentSignificance)
    parser_calc = subparsers.add_parser(
        "calc",
        help="Compute HMC & ES tracks of two binned BedGraphs",
        description="Calculate weighted local correlation and enrichment significance between two BedGraph tracks.",
        epilog=textwrap.dedent("""\
            Usage Example 1:
                localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --output_dir ./results --method locP_and_ES --FDR --binNum_window 11 --step 1 --percentile 90 --percentile_mode all --binNum_peak 3 --FC_thresh 1.5 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2 --threads 4

            Usage Example 2:
                localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --output_dir ./results --percentile 99 --binNum_peak 2 --chrom_sizes hg19.chrom.sizes
            """),
        formatter_class=RawDescriptionHelpFormatter,  # Preserve formatting
    )
    parser_calc.add_argument("--track1", required=True, help="First input BedGraph file.")
    parser_calc.add_argument("--track2", required=True, help="Second input BedGraph file.")
    parser_calc.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for results.",
    )
    parser_calc.add_argument(
        "--method",
        choices=[
            "locP_and_ES",
            "locS_and_ES",
        ],
        default="locP_and_ES",
        help="Method for calculate weighted local correlation and enrichment significance. P and S denote Pearson and Spearman respectively (default: locP_and_ES).",
    )
    parser_calc.add_argument(
        "--FDR",
        action="store_true",
        help="Use Benjamini–Hochberg FDR-corrected q-values instead of raw P-values.",
    )
    parser_calc.add_argument(
        "--binNum_window",
        type=int,
        default=11,
        help="Number of bins in the sliding window (default: 11).",
    )
    parser_calc.add_argument(
        "--step",
        type=int,
        default=1,
        help="Step size for the sliding window (default: 1)",
    )
    parser_calc.add_argument(
        "--percentile",
        type=float,
        default=90,
        help="Percentile for floor correction of low-coverage bins (default: 90). High percentile such as 95 or 98 is recommended, when tracks mainly contains some high sharp peaks, while small percentile like 90 is recommended when tracks mainly contain broad and relatively low peaks.",
    )
    # NEW: small window size for the 3-bin ES calculation
    parser_calc.add_argument(
        "--percentile_mode",
        choices=["all", "nonzero"],
        default="all",
        help="Use all bins or only non-zero bins for percentile (default: all). When choosing 'all' mode, we will choose larger percentile value than that of 'nonzero' mode to get the same floor.",
    )
    parser_calc.add_argument(
        "--binNum_peak",
        type=int,
        default=3,
        help="Number of bins of the peak for ES (default: 3). When the peak is around 400bp and the bin_size=200bp, binNum_peak=2 is recommended.",
    )
    # NEW: FC threshold base for the log-fold enrichment
    parser_calc.add_argument(
        "--FC_thresh",
        type=float,
        default=1.5,
        help="Fold-change threshold used as log base in enrichment (default: 1.5). When FC_thresh=1.5, the null hypothesis is that log1.5(track1/track2)=0, which is quite similar to the FC_thresh in the vocalno plot. Wald, a statistical value following a normal distribution here, euqal to log1.5(track1/track2) / SE can be used to calculate the p value, whose -log10 represents for ES here.",
    )
    parser_calc.add_argument(
        "--norm_method",
        choices=["scale", "cpm", "rpkm"],
        default="rpkm",
        help="Normalisation: scale (match totals) or cpm (counts per million) or rpkm (reads per kilobase per million reads).",
    )
    parser_calc.add_argument(
        "--HMC_scale_pct",
        type=float,
        default=0.9995,
        help="Quantile used to linearly rescale HMC into [0,1] clip at this percentile then divide; default: 0.9995, i.e. top 0.05 percent clipped).",
    )
    parser_calc.add_argument(
        "--chrom_sizes",
        type=str,
        required=True,
        help="Path to the chromosome sizes file.",
    )
    parser_calc.add_argument(
        "--chroms",
        nargs="+",
        default=["all"],
        help="Chromosomes to process (e.g., chr1 chr2). Defaults to 'all'.",
    )
    parser_calc.add_argument(
        "--threads",
        "-t",
        type=int,
        metavar="N",
        default=1,
        help="Worker processes per-chromosome (default: 1).",
    )
    parser_calc.set_defaults(func=calc_corr_main)

    # Subcommand: findreg (alias: find_significantly_different_regions)
    parser_find = subparsers.add_parser(
        "findreg",
        help="Merge significant bins into regions and filter by ES & HMC.",
        description="Merge consecutive significant bins into regions. And find significantly different regions from ES & HMC tracks.",
        epilog=textwrap.dedent("""\
            Example 1:
              localfinder findreg --track_E track_ES.bedgraph --track_C track_HMC.bedgraph --output_dir ./findreg_out --p_thresh 0.05 --binNum_thresh 2 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

            Example 2:
              localfinder findreg --track_E track_ES.bedgraph --track_C track_HMC.bedgraph --output_dir ./findreg_out --chrom_sizes hg19.chrom.sizes
        """),
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser_find.add_argument("--track_E", required=True, help="track_ES.bedgraph")
    parser_find.add_argument("--track_C", required=True, help="track_HMC.bedgraph")
    parser_find.add_argument("--output_dir", required=True)
    parser_find.add_argument(
        "--p_thresh",
        type=float,
        default=0.05,
        help="P-value threshold (default: 0.05)",
    )
    parser_find.add_argument(
        "--binNum_thresh",
        type=int,
        default=2,
        help="Min consecutive significant bins per region (default: 2)",
    )
    parser_find.add_argument(
        "--max_gap_bins",
        type=int,
        default=0,  ### <<< NEW
        help="Allow merging two significant runs on the same chromosome if the empty gap between them is ≤ this number of bins (default: 0).",
    )
    parser_find.add_argument("--chroms", nargs="+", default=["all"])
    parser_find.add_argument("--chrom_sizes", required=True)
    parser_find.set_defaults(func=find_regions_main)

    # Subcommand: viz (alias: visualize_tracks_or_scatters)
    parser_visualize = subparsers.add_parser(
        "viz",
        help="Visualize genomic tracks.",
        description="Visualize genomic tracks.",
        epilog=textwrap.dedent("""\
            Usage Example 1:
                localfinder viz --input_files track1.bedgraph track2.bedgraph --output_file output.html --method plotly --region chr1 1000000 2000000 --colors blue red

            Usage Example 2:
                localfinder viz --input_files track1.bedgraph track2.bedgraph --output_file output.png --method pyGenomeTracks --region chr1 1000000 2000000 --colors
            """),
        formatter_class=RawDescriptionHelpFormatter,  # Preserve formatting
    )
    parser_visualize.add_argument(
        "--input_files",
        nargs="+",
        required=True,
        help="Input BedGraph files to visualize.",
    )
    parser_visualize.add_argument(
        "--output_file",
        required=True,
        help="Output visualization file (e.g., PNG, HTML).",
    )
    parser_visualize.add_argument(
        "--method",
        choices=["pyGenomeTracks", "plotly"],
        required=True,
        help="Visualization method to use.",
    )
    parser_visualize.add_argument(
        "--region",
        nargs=3,
        metavar=("CHROM", "START", "END"),
        help="Region to visualize in the format: CHROM START END (e.g., chr20 1000000 2000000).",
    )
    parser_visualize.add_argument(
        "--colors",
        nargs="+",
        help="Colors for the tracks (optional).",
    )
    parser_visualize.set_defaults(func=visualize_main)

    # Subcommand: pipeline
    parser_pipe = subparsers.add_parser(
        "pipeline",
        help="Run bin → calc → findreg in one command.",
        description="Run localfinder sequentially.",
        epilog=textwrap.dedent("""\
            Usage Example 1:
                localfinder pipeline --input_files track1.bedgraph track2.bedgraph --output_dir ./results --chrom_sizes hg19.chrom.sizes --bin_size 200 --method locP_and_ES --FDR --binNum_window 11 --binNum_peak 3 --step 1 --percentile 90 --percentile_mode all --FC_thresh 1.5 --norm_method rpkm --p_thresh 0.05 --binNum_thresh 2 --chroms chr1 chr2 --threads 4

            Usage Example 2:
                localfinder pipeline --input_files track1.bigwig track2.bigwig --output_dir ./results --chrom_sizes hg19.chrom.sizes --binNum_peak 3 --percentile 95 --binNum_thresh 2
            """),
        formatter_class=RawDescriptionHelpFormatter,  # Preserve formatting
    )

    parser_pipe.add_argument(
        "--input_files",
        nargs="+",
        required=True,
        help="Input files in BigWig/BedGraph format.",
    )
    parser_pipe.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for the pipeline results.",
    )
    parser_pipe.add_argument(
        "--chrom_sizes",
        required=True,
        help="Path to the chromosome sizes file.",
    )
    parser_pipe.add_argument(
        "--bin_size",
        type=int,
        default=200,
        help="Size of each bin (default: 200).",
    )
    # calc options forwarded
    parser_pipe.add_argument(
        "--method",
        choices=["locP_and_ES", "locS_and_ES"],
        default="locP_and_ES",
        help="Method for calculate weighted local correlation and enrichment significance. P and S denote Pearson and Spearman respectively (default: locP_and_ES).",
    )
    parser_pipe.add_argument(
        "--FDR",
        action="store_true",
        help="Use Benjamini–Hochberg FDR-corrected q-values instead of raw P-values",
    )
    parser_pipe.add_argument(
        "--binNum_window",
        type=int,
        default=11,
        help="Number of bins in the sliding window (default: 11).",
    )
    parser_pipe.add_argument(
        "--binNum_peak",
        type=int,
        default=3,
        help="Number of bins of the peak for ES (default: 3). When the peak is around 400bp and the bin_size=200bp, binNum_peak=2 is recommended.",
    )
    parser_pipe.add_argument(
        "--step",
        type=int,
        default=1,
        help="Step size for the sliding window (default: 1)",
    )
    parser_pipe.add_argument(
        "--percentile",
        type=float,
        default=90,
        help="Percentile for floor correction of low-coverage bins (default: 90). High percentile such as 95 or 98 is recommended, when tracks mainly contains some high sharp peaks, while small percentile like 90 is recommended when tracks mainly contain broad and relatively low peaks.",
    )
    parser_pipe.add_argument(
        "--percentile_mode",
        choices=["all", "nonzero"],
        default="all",
        help='Use all bins or only non-zero bins for percentile (default: all). When choosing "all" mode, we will choose larger percentile value than that of "nonzero" mode to get the same floor.',
    )
    parser_pipe.add_argument(
        "--FC_thresh",
        type=float,
        default=1.5,
        help="Fold-change threshold used as log base in enrichment (default: 1.5). When FC_thresh=1.5, the null hypothesis is that log1.5(track1/track2)=0, which is quite similar to the FC_thresh in the volcano plot. Wald, a statistical value following a normal distribution here, equal to log1.5(track1/track2) / SE can be used to calculate the p value, whose -log10 represents for ES here.",
    )
    parser_pipe.add_argument(
        "--norm_method",
        choices=["scale", "cpm", "rpkm"],
        default="rpkm",
        help="Normalisation: scale (match totals) or cpm (counts per million) or rpkm (reads per kilobase per million reads).",
    )
    parser_pipe.add_argument(
        "--HMC_scale_pct",
        type=float,
        default=0.9995,
        help="Quantile used to linearly rescale HMC into [0,1] clip at this percentile then divide; default: 0.9995, i.e. top 0.05 percent clipped).",
    )
    # findreg thresholds forwarded
    parser_pipe.add_argument(
        "--p_thresh",
        type=float,
        default=0.05,
        help="P-value threshold for merging significant bins into regions (default: 0.05)",
    )
    parser_pipe.add_argument(
        "--binNum_thresh",
        type=int,
        default=2,
        help="Min consecutive significant bins per region",
    )
    parser_pipe.add_argument(
        "--max_gap_bins",
        type=int,
        default=0,  ### <<< NEW
        help="Allow merging two significant runs on the same chromosome if the empty gap between them is ≤ this number of bins (default: 0).",
    )
    parser_pipe.add_argument(
        "--chroms",
        nargs="+",
        default=["all"],
        help='Chromosomes to process (e.g., chr1 chr2). Defaults to "all".',
    )
    parser_pipe.add_argument(
        "--threads",
        "-t",
        type=int,
        default=1,
        help="Number of worker processes for bin & calc (default: 1)",
    )
    parser_pipe.set_defaults(func=run_pipeline)

    # Parse the arguments
    args = parser.parse_args()

    # Execute the appropriate function based on the subcommand
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
