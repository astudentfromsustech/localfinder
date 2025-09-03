import os
import shutil
import sys

import pandas as pd

from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

from localfinder.utils import get_chromosomes_from_chrom_sizes, locCor_and_ES


# ----------------------------------------------------------------------
# helper for multiprocessing
# ----------------------------------------------------------------------
def _calc_one_chrom(
    df,
    chrom,
    bin_number_of_window,
    step,
    percentile,
    percentile_mode,
    FC_thresh,
    bin_number_of_peak,
    norm_method,
    corr_method,
    FDR,
    HMC_scale_pct,
    output_dir,
):
    """Run locCor_and_ES for a single chromosome and return its two file paths."""
    print(f"\n=== {chrom} ===")
    locCor_and_ES(
        df,
        bin_number_of_window=bin_number_of_window,
        step=step,
        percentile=percentile,
        percentile_mode=percentile_mode,
        FC_thresh=FC_thresh,
        bin_number_of_peak=bin_number_of_peak,
        norm_method=norm_method,
        corr_method=corr_method,
        FDR=FDR,
        HMC_scale_pct=HMC_scale_pct,
        output_dir=output_dir,
        chrom=chrom,
    )
    es_path = os.path.join(output_dir, f"track_ES.{chrom}.bedgraph")
    HMC_path = os.path.join(output_dir, f"track_HMC.{chrom}.bedgraph")
    return chrom, es_path, HMC_path


def main(args):
    track1_file = args.track1
    track2_file = args.track2
    output_dir = args.output_dir
    method = args.method
    FDR = args.FDR
    percentile = args.percentile
    percentile_mode = getattr(args, "percentile_mode", "all")  # --- NEW ---
    bin_number_of_window = args.binNum_window
    bin_number_of_peak = args.binNum_peak
    FC_thresh = args.FC_thresh
    step = args.step
    chroms = args.chroms
    chrom_sizes = args.chrom_sizes
    HMC_scale_pct = getattr(args, "HMC_scale_pct", 0.9995)
    norm_method = getattr(args, "norm_method", "rpkm")  # --- NEW ---
    n_threads = getattr(args, "threads", 1)  # ← NEW

    os.makedirs(output_dir, exist_ok=True)

    # **Modification Start**
    # If chroms is 'all', retrieve all chromosomes from chrom_sizes
    if chroms == ["all"] or chroms is None:
        chroms = get_chromosomes_from_chrom_sizes(chrom_sizes)
        print(f"'chroms' set to all chromosomes from chrom_sizes: {chroms}")
    else:
        print(f"'chroms' set to specified chromosomes: {chroms}")
    # **Modification End**

    # Read the binned tracks
    try:
        df1 = pd.read_csv(
            track1_file,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "readNum_1"],
        )
        df2 = pd.read_csv(
            track2_file,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "readNum_2"],
        )
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Merge the dataframes
    df = pd.merge(df1, df2, on=["chr", "start", "end"], how="inner")

    # decide which correlation to use
    if method == "locP_and_ES":
        corr_method = "pearson"
    elif method == "locS_and_ES":
        corr_method = "spearman"
    else:
        print(f"Unsupported method: {method}")
        sys.exit(1)

    # ---------------------------------------------------------------
    #  parallel execution per chromosome
    # ---------------------------------------------------------------
    worker = partial(
        _calc_one_chrom,
        df,
        bin_number_of_window=bin_number_of_window,
        step=step,
        percentile=percentile,
        percentile_mode=percentile_mode,
        FC_thresh=FC_thresh,
        bin_number_of_peak=bin_number_of_peak,
        norm_method=norm_method,
        corr_method=corr_method,
        FDR=FDR,
        output_dir=output_dir,
        HMC_scale_pct=HMC_scale_pct,
    )

    produced = {}  # --- NEW ---
    with ProcessPoolExecutor(max_workers=n_threads) as pool:
        futures = {pool.submit(worker, chrom): chrom for chrom in chroms}
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                chrom_ret, es_path, HMC_path = fut.result()
                produced[chrom_ret] = (
                    es_path,
                    HMC_path,
                )  # --- NEW ---
                print(f"[DONE] {chrom_ret}")
            except Exception as e:
                print(f"[ERROR] {chrom}: {e}")
                raise

    # ---------------------------------------------------------------
    #  concatenate per-chrom files → combined BedGraphs, then delete
    # ---------------------------------------------------------------
    combo_ES = os.path.join(output_dir, "track_ES.bedgraph")
    combo_HMC = os.path.join(output_dir, "track_HMC.bedgraph")

    print("[COMBINE] building combined BedGraphs")
    with (
        open(combo_ES, "wb") as es_out,
        open(combo_HMC, "wb") as HMC_out,
    ):
        for chrom in chroms:  # preserve order
            es_path, HMC_path = produced.get(chrom, (None, None))  ### <<< NEW guard
            for src, dst in [(es_path, es_out), (HMC_path, HMC_out)]:
                if src is None:  # --- NEW ---
                    # this chrom was skipped entirely (e.g. not in BigWig)
                    print(f"[SKIP] {chrom}: no {os.path.basename(dst.name)}")
                    continue
                if not os.path.exists(src):  # --- NEW ---
                    print(f"[SKIP] {src} missing – nothing to append")
                    continue
                with open(src, "rb") as fh:  # UNCHANGED
                    shutil.copyfileobj(fh, dst)
                os.remove(src)  # UNCHANGED

    print(f"[COMBINE] saved {combo_ES} and {combo_HMC}")
