# File: commands/simbg.py

import os, sys, shutil, pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from localfinder.utils import (
    build_pseudo_background_bedgraph,
    get_chromosomes_from_chrom_sizes,   
)

# ----------------------------------------------------------------------
# helper for multiprocessing: one chromosome
# ----------------------------------------------------------------------
def _build_bg_one_chrom(
    input_bedgraph,
    chrom,
    bin_number_of_window,
    bg_method,
    output_dir,
):
    """
    Build pseudo background for a single chromosome and return its file path.
    """
    basename = os.path.splitext(os.path.basename(input_bedgraph))[0]
    out_path = os.path.join(output_dir, f"{basename}.bg.{chrom}.bedgraph")

    # Call utils; some versions may return (bg_path, es_path). We only need bg.
    ret = build_pseudo_background_bedgraph(
        input_bedgraph=input_bedgraph,
        output_bedgraph=out_path,
        bin_number_of_window=bin_number_of_window,
        bg_method=bg_method,
        chrom=chrom,
    )
    # return to bg_path
    if isinstance(ret, tuple):
        bg_path = ret[0]
    else:
        bg_path = ret

    print(f"[DONE] {chrom}")
    return chrom, bg_path


def main(args):
    """
    Simulate background for an input BedGraph using a multi-scale local model.

    Args (Namespace):
      - input_bedgraph: path to the source BedGraph (binned signal)
      - output_dir:     where to write outputs
      - binNum_window:  odd int; base window W (we also use (W-1)*10+1 and (W-1)*100+1)
      - bg_method:      'poisson' (or 'zinb' if your utils supports it)
      - chroms:         list like ['chr1','chr2'] or ['all'] or None
      - chrom_sizes:    path to chrom.sizes (used only when chroms == ['all'] or None)
      - threads:        int, number of worker processes
    """
    input_bedgraph      = args.input_bedgraph
    output_dir          = args.output_dir
    bin_number_of_window= getattr(args, 'binNum_window', 11)
    bg_method           = getattr(args, 'bg_method', 'poisson')
    chroms              = getattr(args, 'chroms', None)
    chrom_sizes         = getattr(args, 'chrom_sizes', None)
    n_threads           = getattr(args, 'threads', 1)

    if not os.path.exists(input_bedgraph):
        print(f"Error: input BedGraph not found: {input_bedgraph}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    # -----------------------------
    # resolve chromosome list
    # -----------------------------
    if chroms == ['all'] or chroms is None:
        if chrom_sizes:
            chroms = get_chromosomes_from_chrom_sizes(chrom_sizes)
            print(f"'chroms' set to all chromosomes from chrom_sizes: {chroms}")
        else:
            # infer chroms from the BedGraph itself
            print("[INFO] inferring chromosomes from input BedGraph…")
            try:
                df_head = pd.read_csv(input_bedgraph, sep='\t', header=None,
                                      names=['chr','start','end','value'], usecols=[0])
                chroms = list(pd.unique(df_head['chr']))
                print(f"'chroms' inferred from input: {chroms}")
            except Exception as e:
                print(f"Error inferring chromosomes: {e}")
                sys.exit(1)
    else:
        print(f"'chroms' set to specified chromosomes: {chroms}")

    # ---------------------------------------------------------------
    # parallel execution per chromosome
    # ---------------------------------------------------------------
    worker = partial(
        _build_bg_one_chrom,
        input_bedgraph,
        bin_number_of_window=bin_number_of_window,
        bg_method=bg_method,
        output_dir=output_dir,
    )

    produced = {}
    with ProcessPoolExecutor(max_workers=n_threads) as pool:
        futures = {pool.submit(worker, chrom): chrom for chrom in chroms}
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                chrom_ret, bg_path = fut.result()
                produced[chrom_ret] = bg_path
            except Exception as e:
                print(f"[ERROR] {chrom}: {e}")
                raise

    # ---------------------------------------------------------------
    # concatenate per-chrom outputs -> combined BedGraph, then delete
    # ---------------------------------------------------------------
    base = os.path.splitext(os.path.basename(input_bedgraph))[0]
    combo_bg = os.path.join(output_dir, f"{base}.bg.bedgraph")

    print("[COMBINE] building combined background BedGraph")
    with open(combo_bg, "wb") as outfh:
        for chrom in chroms:  # preserve requested order
            src = produced.get(chrom)
            if not src:
                print(f"[SKIP] {chrom}: no background file recorded")
                continue
            if not os.path.exists(src):
                print(f"[SKIP] {src} missing – nothing to append")
                continue
            with open(src, "rb") as fh:
                shutil.copyfileobj(fh, outfh)
            os.remove(src)

    print(f"[COMBINE] saved {combo_bg}")
