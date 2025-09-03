import os
import shutil

from concurrent.futures import ProcessPoolExecutor, as_completed

from localfinder.utils import (
    check_external_tools,
    get_chromosomes_from_chrom_sizes,
    process_and_bin_file,
)


# ----------------------------------------------------------------------
# helper placed at module level so it is picklable by multiprocessing
# ----------------------------------------------------------------------
def _run_one_chrom(input_file, output_dir, stem, bin_size, chrom_sizes, chrom):
    """
    Worker that bins one chromosome of one input file.
    Written at top-level so ProcessPool can pickle it.
    """
    out_name = f"{stem}.binSize{bin_size}.{chrom}.bedgraph"
    output_file = os.path.join(output_dir, out_name)

    print(f"[WORKER] {chrom} → {out_name}")
    process_and_bin_file(
        input_file=input_file,
        output_file=output_file,
        bin_size=bin_size,
        chrom_sizes=chrom_sizes,
        chrom=chrom,
    )
    return chrom, output_file


def main(args):
    input_files = args.input_files
    output_dir = args.output_dir
    bin_size = args.bin_size
    chrom_sizes = args.chrom_sizes
    chroms = args.chroms  # list or ['all']
    n_threads = args.threads
    # Ensure required external tools are available
    check_external_tools()

    os.makedirs(output_dir, exist_ok=True)

    if chroms == ["all"] or chroms is None:
        chroms = get_chromosomes_from_chrom_sizes(chrom_sizes)
        print(f"[INFO] Using ALL chromosomes from {chrom_sizes}: {chroms}")
    else:
        print(f"[INFO] Processing specified chromosomes: {chroms}")

    tasks = []
    stems = {}
    for input_file in input_files:
        print(f"[START] Processing {input_file}")

        base = os.path.basename(input_file)
        stem, _ = os.path.splitext(base)
        stems[input_file] = stem
        for chrom in chroms:
            tasks.append((input_file, stem, chrom))

    # track outputs per file → per chromosome
    produced = {inf: {} for inf in input_files}  # --- NEW ---

    with ProcessPoolExecutor(max_workers=n_threads) as pool:  # --- CHANGED ---
        futures = {
            pool.submit(
                _run_one_chrom,
                in_file,
                output_dir,
                stem,
                bin_size,
                chrom_sizes,
                chrom,
            ): (in_file, chrom)
            for in_file, stem, chrom in tasks  # --- NEW ---
        }

        for fut in as_completed(futures):
            in_file, chrom = futures[fut]
            try:
                chrom_ret, out_path = fut.result()
                produced[in_file][chrom_ret] = out_path  # ← store path
                print(f"[DONE] {chrom_ret} ({os.path.basename(in_file)}) → {out_path}")
            except Exception as e:
                print(f"[ERROR] {chrom} ({in_file}): {e}")
                raise

    # ------------------------------------------------------------------
    # concatenate per-chromosome bedGraphs → combined bedGraph
    # ------------------------------------------------------------------
    for in_file in input_files:  # --- NEW ---
        stem = stems[in_file]
        combo_name = f"{stem}.binSize{bin_size}.bedgraph"
        combo_path = os.path.join(output_dir, combo_name)

        print(f"[COMBINE] creating {combo_name}")
        with open(combo_path, "wb") as out_fh:
            for chrom in chroms:  # keep requested order
                part = produced[in_file].get(chrom)
                if part and os.path.exists(part):
                    with open(part, "rb") as part_fh:
                        shutil.copyfileobj(part_fh, out_fh)
                    os.remove(part)  # ← delete per-chrom file
        print(f"[COMBINE] saved {combo_path}")

    print("[FINISHED] all inputs\n")
