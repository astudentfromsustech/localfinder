import os
import shutil
import subprocess
import sys
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from plotly.subplots import make_subplots
from scipy.stats import norm


def check_external_tools():
    required_tools = ["bedtools", "bigWigToBedGraph", "samtools"]
    missing_tools = []
    for tool in required_tools:
        if shutil.which(tool) is None:
            missing_tools.append(tool)

    if missing_tools:
        print(
            "Error: The following required external tools are not installed or not found in PATH:"
        )
        for tool in missing_tools:
            print(f" - {tool}")
        print("Please install them and ensure they are available in your system PATH.")
        sys.exit(1)


def get_chromosomes_from_chrom_sizes(chrom_sizes_file):
    """
    Retrieves all chromosome names from the chromosome sizes file.

    Parameters:
    - chrom_sizes_file (str): Path to the chromosome sizes file.

    Returns:
    - chromosomes (list): List of chromosome names.
    """
    try:
        df_chrom = pd.read_csv(
            chrom_sizes_file,
            sep="\t",
            header=None,
            names=["chr", "size"],
        )
        chromosomes = df_chrom["chr"].tolist()
        return chromosomes
    except Exception as e:
        print(f"Error reading chromosome sizes file {chrom_sizes_file}: {e}")
        sys.exit(1)


def process_and_bin_file(input_file, output_file, bin_size, chrom_sizes, chrom):
    """
    Processes an input file by detecting its format, converting it to BedGraph if necessary,
    and binning it into fixed-size bins for a single chromosome.

    Parameters
    ----------
    input_file : str
        Path to the input BigWig or BedGraph file.
    output_file : str
        Path for the final binned BedGraph file.
    bin_size : int
        Bin size in base pairs.
    chrom_sizes : str
        Path to the chromosome-sizes file (2-column TSV).
    chrom : str
        Chromosome to process (e.g. 'chr1').
    """

    print(f"[START] {chrom}: processing {os.path.basename(input_file)}")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_bedgraph = os.path.join(temp_dir, "temp.bedgraph")

        input_format = detect_format(input_file)
        if input_format is None:
            print(f"Could not determine the format of the file: {input_file}")
            sys.exit(1)

        try:
            # ------------------------------------------------------------------
            # 1. Convert to BedGraph (single chromosome)
            # ------------------------------------------------------------------
            if input_format == "bigwig":
                subprocess.run(
                    [
                        "bigWigToBedGraph",
                        f"-chrom={chrom}",
                        input_file,
                        temp_bedgraph,
                    ],
                    check=True,
                )
            elif input_format == "bedgraph":
                with (
                    open(input_file) as inp,
                    open(temp_bedgraph, "w") as outp,
                ):
                    for line in inp:
                        if line.startswith(f"{chrom}\t"):
                            outp.write(line)
                if os.path.getsize(temp_bedgraph) == 0:
                    print(f"No data for chromosome {chrom} in {input_file}")
                    sys.exit(1)
            elif input_format == "bam":
                with open(temp_bedgraph, "w") as out_bg:
                    # 1) stream just this chromosome
                    p_view = subprocess.Popen(
                        ["samtools", "view", "-b", input_file, chrom],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,  # keep stderr in case of error
                    )
                    # 2) genomecov on the stream
                    p_cov = subprocess.run(
                        [
                            "bedtools",
                            "genomecov",
                            "-bg",
                            "-ibam",
                            "-",
                        ],
                        stdin=p_view.stdout,
                        stdout=out_bg,
                        stderr=subprocess.PIPE,
                    )
                    p_view.stdout.close()  # allow p_view to receive SIGPIPE if genomecov fails
                    p_view.wait()

                    if p_view.returncode != 0:
                        raise RuntimeError(
                            f"samtools view failed ({p_view.returncode}):\n{p_view.stderr.read().decode()}"
                        )
                    if p_cov.returncode != 0:
                        raise RuntimeError(
                            f"bedtools genomecov failed ({p_cov.returncode}):\n{p_cov.stderr.decode()}"
                        )
            else:
                print(f"Unsupported format: {input_format}")
                sys.exit(1)

            # ------------------------------------------------------------------
            # 2. Bin the BedGraph data
            # ------------------------------------------------------------------
            print(f"[INFO] Binning {chrom} at {bin_size}-bp resolution")
            bin_bedgraph(
                temp_bedgraph,
                output_file,
                bin_size,
                chrom_sizes,
                chrom,
            )

        finally:
            # Temporary directory and its contents are deleted automatically
            print(f"[DONE] {chrom}: intermediate files cleaned up")


def detect_format(filename):
    """
    Infers the file format based on the file extension.

    Parameters:
    - filename (str): The name or path of the file.

    Returns:
    - format (str): The inferred format ('bam', 'bedgraph', 'bigwig'), or None if unknown.
    """
    extension = os.path.splitext(filename)[1].lower()
    if extension == ".bam":
        return "bam"
    elif extension in [".bedgraph", ".bdg"]:
        return "bedgraph"
    elif extension in [".bigwig", ".bw"]:
        return "bigwig"
    else:
        return None


def bin_bedgraph(input_bedgraph, output_bedgraph, bin_size, chrom_sizes, chrom):
    """
    Bins BedGraph data for a single chromosome into fixed-size bins,
    replacing missing values with zero.  Windows produced by
    `bedtools makewindows` are already coordinate-sorted, so no extra sort
    is run on them; we still sort the BedGraph itself for safety.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        # --------------------------------------------------------------
        # 1. Tiny genome file for the target chromosome
        # --------------------------------------------------------------
        temp_genome = os.path.join(temp_dir, "temp_genome.txt")
        df_sizes = pd.read_csv(
            chrom_sizes,
            sep="\t",
            header=None,
            names=["chrom", "size"],
        )
        row = df_sizes[df_sizes["chrom"] == chrom]
        if row.empty:
            print(f"Chromosome {chrom} not found in {chrom_sizes}")
            sys.exit(1)
        with open(temp_genome, "w") as g:
            g.write(f"{chrom}\t{int(row.iloc[0]['size'])}\n")

        # Temporary paths
        temp_bins = os.path.join(temp_dir, "temp_bins.bed")
        temp_bg_sorted = os.path.join(temp_dir, "temp_bedgraph_sorted.bedgraph")
        temp_binned = os.path.join(temp_dir, "temp_binned.bedgraph")

        # --------------------------------------------------------------
        # 2. Generate fixed-size windows (already sorted)
        # --------------------------------------------------------------
        with open(temp_bins, "w") as f_bins:
            subprocess.run(
                [
                    "bedtools",
                    "makewindows",
                    "-g",
                    temp_genome,
                    "-w",
                    str(bin_size),
                ],
                stdout=f_bins,
                check=True,
            )

        # --------------------------------------------------------------
        # 3. Coordinate-sort the BedGraph
        # --------------------------------------------------------------
        with open(temp_bg_sorted, "w") as f_bg_sorted:
            subprocess.run(
                [
                    "bedtools",
                    "sort",
                    "-faidx",
                    temp_genome,
                    "-i",
                    input_bedgraph,
                ],
                stdout=f_bg_sorted,
                check=True,
            )

        # --------------------------------------------------------------
        # 4. Map values onto bins (`-sorted` because both inputs are sorted)
        # --------------------------------------------------------------
        with open(temp_binned, "w") as f_binned:
            subprocess.run(
                [
                    "bedtools",
                    "map",
                    "-sorted",
                    "-a",
                    temp_bins,
                    "-b",
                    temp_bg_sorted,
                    "-c",
                    "4",
                    "-o",
                    "mean",
                ],
                stdout=f_binned,
                check=True,
            )

        # --------------------------------------------------------------
        # 5. Fill missing values with zero and save
        # --------------------------------------------------------------
        df = pd.read_csv(
            temp_binned,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "value"],
            na_values=".",
        )
        df["value"] = df["value"].fillna(0).astype(float)
        df.to_csv(output_bedgraph, sep="\t", header=False, index=False)
        print(f"[SUCCESS] ({chrom}) Saved {output_bedgraph}")


def locCor_and_ES(
    df,
    column1="readNum_1",
    column2="readNum_2",
    bin_number_of_window=11,
    step=1,
    percentile=5,
    percentile_mode="all",
    FC_thresh=1.5,
    bin_number_of_peak=11,
    norm_method="rpkm",
    corr_method="pearson",
    FDR=False,
    HMC_scale_pct=0.9995,
    output_dir="output",
    chrom=None,
):
    """
    Calculate weighted local correlation and enrichment significance.
    Writes two bedgraph files: track_HMC.bedgraph, track_ES.bedgraph.
    """
    if chrom is None:
        raise ValueError("Argument 'chrom' (e.g. 'chr1') must be provided")

    print(
        f"calculate weighted local correlation and enrichment significance "
        f"(corr_method={corr_method})"
    )
    print(
        f"parameters: percentile={percentile}, FC_thresh={FC_thresh}, "
        f"bin_number_of_window={bin_number_of_window}, "
        f"bin_number_of_peak={bin_number_of_peak}, "
        f"FDR={'ON' if FDR else 'OFF'}",
        f"norm_method={norm_method}, HMC_scale_pct={HMC_scale_pct}",
    )

    EPS = 1e-9
    os.makedirs(output_dir, exist_ok=True)

    # ---------- 1  global normalisation ---------------------------------
    print(f"step1: global normalisation ({norm_method})")
    cov1, cov2 = df[column1].sum(), df[column2].sum()
    if norm_method == "scale":
        # existing “match smaller total” behavior
        if cov1 and cov2:
            if cov1 > cov2:
                df[column1] *= cov2 / cov1
                print(f"Scaled {column1} down to match {column2}")
            else:
                df[column2] *= cov1 / cov2
                print(f"Scaled {column2} down to match {column1}")
    elif norm_method == "cpm":
        # counts per million: each column → counts / (total/1e6)
        if cov1 > 0:
            df[column1] = df[column1] * 1e6 / cov1
        if cov2 > 0:
            df[column2] = df[column2] * 1e6 / cov2
        print(f"Converted {column1} and {column2} to CPM")

    elif norm_method == "rpkm":  ### <<< NEW BLOCK
        bin_lengths_kb = (df["end"] - df["start"]) / 1_000.0
        if cov1 > 0:
            df[column1] = df[column1] / bin_lengths_kb * 1e6 / cov1
        if cov2 > 0:
            df[column2] = df[column2] / bin_lengths_kb * 1e6 / cov2
        print(f"Converted {column1} and {column2} to RPKM")
    else:
        raise ValueError(f"Unknown norm_method: {norm_method}")

    # ---------- output paths -------------------------------------------
    out_ES = os.path.join(output_dir, f"track_ES.{chrom}.bedgraph")
    out_HMC = os.path.join(output_dir, f"track_HMC.{chrom}.bedgraph")

    half_w = (bin_number_of_window - 1) // 2
    half_p = (bin_number_of_peak - 1) // 2

    df_final = pd.DataFrame()

    # ===================================================================
    # per-chromosome processing
    # ===================================================================
    print("step2: calculate weighted correlation and ES")
    df_raw = df[df["chr"] == chrom].reset_index(drop=True)  # --- CHANGED ---
    n = len(df_raw)
    if n == 0:
        print(f"[WARN] No data for {chrom}; nothing written")
        return

    print(f"\nProcessing {chrom}: {n} bins")
    print(f"  half_window = {half_w}, half_peak = {half_p}")

    # 2·1  percentile floor (only for ES)
    all_vals = pd.concat([df_raw[column1], df_raw[column2]])
    if percentile_mode == "nonzero":
        vals = all_vals[all_vals > 0]
    else:
        vals = all_vals
    if len(vals):
        pct = np.percentile(vals, percentile)
    else:
        pct = 0.0
    print(
        f"  {percentile}th percentile over "
        f"{'non-zero' if percentile_mode == 'nonzero' else 'all'} bins = {pct:.4f}"
    )

    df_floor = df_raw.copy()  ### floor
    mask = (df_floor[column1] <= pct) & (df_floor[column2] <= pct)
    df_floor.loc[mask, [column1, column2]] = pct
    df_floor[column1] = df_floor[column1].astype("float64")
    df_floor[column2] = df_floor[column2].astype("float64")

    # ---------- 2·2  sliding-window statistics (REFINED) -----------------
    print(
        f"Step2.2: calculate weighted local correlation ({corr_method}) "
        f"and ES (window={bin_number_of_window}, step={step})"
    )

    W = 2 * ((bin_number_of_window - 1) // 2) + 1  # full window length
    P = 2 * ((bin_number_of_peak - 1) // 2) + 1  # peak window length

    # Convert to NumPy arrays once
    x = df_raw[column1].to_numpy(float)
    y = df_raw[column2].to_numpy(float)
    xf = df_floor[column1].to_numpy(float)
    yf = df_floor[column2].to_numpy(float)

    # Prefix sums (pad with leading zero for easier slicing)
    sx = np.concatenate(([0.0], x.cumsum()))
    sy = np.concatenate(([0.0], y.cumsum()))
    sxy = np.concatenate(([0.0], (x * y).cumsum()))
    sx2 = np.concatenate(([0.0], (x * x).cumsum()))
    sy2 = np.concatenate(([0.0], (y * y).cumsum()))

    sxf = np.concatenate(([0.0], xf.cumsum()))
    syf = np.concatenate(([0.0], yf.cumsum()))
    sx2f = np.concatenate(([0.0], (xf * xf).cumsum()))
    sy2f = np.concatenate(([0.0], (yf * yf).cumsum()))

    # Allocate result arrays (same names as before)
    corr = np.zeros(n, float)
    m_corr = np.zeros(n, float)
    hmw = np.zeros(n, float)

    m1_raw = np.zeros(n, float)
    m2_raw = np.zeros(n, float)

    mES1_fl = np.zeros(n, float)
    mES2_fl = np.zeros(n, float)

    m1_fl = np.zeros(n, float)
    m2_fl = np.zeros(n, float)
    v1_fl = np.zeros(n, float)
    v2_fl = np.zeros(n, float)
    d1_fl = np.zeros(n, float)
    d2_fl = np.zeros(n, float)

    # Helper: prefix-sum slice mean / var
    def box_mean(pref, lo, hi, L):
        return (pref[hi] - pref[lo]) / L

    def box_var(pref, pref2, lo, hi, L, mean):
        return (pref2[hi] - pref2[lo]) / L - mean * mean

    half_w = (bin_number_of_window - 1) // 2
    half_p = (bin_number_of_peak - 1) // 2

    for i in range(half_w, n - half_w, step):
        lo_w = i - half_w
        hi_w = i + half_w + 1
        lo_p = i - half_p
        hi_p = i + half_p + 1

        # RAW means & variances for correlation window
        mx = box_mean(sx, lo_w, hi_w, W)
        my = box_mean(sy, lo_w, hi_w, W)
        vx = box_var(sx, sx2, lo_w, hi_w, W, mx)
        vy = box_var(sy, sy2, lo_w, hi_w, W, my)
        cxy = box_mean(sxy, lo_w, hi_w, W) - mx * my

        if vx > 0 and vy > 0:
            r = cxy / np.sqrt(vx * vy)
            corr[i] = max(r, 0.0)

        m1_raw[i], m2_raw[i] = mx, my
        hmw[i] = (mx * my) / (mx + my + EPS)

        # FLOORED means / vars in same window (for SE)
        m1f = box_mean(sxf, lo_w, hi_w, W)
        m2f = box_mean(syf, lo_w, hi_w, W)
        v1f = box_var(sxf, sx2f, lo_w, hi_w, W, m1f)
        v2f = box_var(syf, sy2f, lo_w, hi_w, W, m2f)

        m1_fl[i], m2_fl[i] = m1f, m2f
        v1_fl[i], v2_fl[i] = v1f, v2f
        d1_fl[i] = max((v1f - m1f) / (m1f**2 + EPS), 0.0)
        d2_fl[i] = max((v2f - m2f) / (m2f**2 + EPS), 0.0)

        # FLOORED means in peak region
        mES1_fl[i] = box_mean(sxf, lo_p, hi_p, P)
        mES2_fl[i] = box_mean(syf, lo_p, hi_p, P)

    # 2·3  calculate mean correlation
    for i in range(half_w, n - half_w, step):
        m_corr[i] = corr[i - half_w : i + half_w + 1].mean()

    # ---------- logFC in peak region -----------------------------
    safe1 = np.maximum(mES1_fl, pct)
    safe2 = np.maximum(mES2_fl, pct)
    logFC = np.log(safe2 / (safe1 + EPS)) / np.log(FC_thresh)
    logFC[:half_w] = 0.0
    logFC[-half_w:] = 0.0

    idx = np.arange(half_w, n - half_w, step)
    mu1 = np.maximum(m1_fl[idx], pct)
    mu2 = np.maximum(m2_fl[idx], pct)
    SE = np.sqrt(1 / mu1 + 1 / mu2 + d1_fl[idx] + d2_fl[idx])  ### <<< CHANGED
    Wald = logFC[idx] / SE
    p = 2 * (1 - norm.cdf(np.abs(Wald)))
    if FDR:  # ── new branch ──
        # q  = multipletests(p, alpha=0.05, method='fdr_bh')[1]
        # lP = -np.log10(np.where(q == 0, np.nan, q))    # log-q
        lP = -np.log10(np.where(p == 0, np.nan, p))  # log-p
    else:
        lP = -np.log10(np.where(p == 0, np.nan, p))  # log-p

    # ---------- assemble dataframe ----------
    dfc = df_raw.copy()
    dfc["corr"] = corr
    dfc["m_corr"] = m_corr
    dfc["hmw"] = hmw
    dfc["logFC"] = logFC
    dfc["SE"] = 0.0
    dfc["Wald"] = 0.0
    dfc["Wald_pValue"] = 0.0
    dfc["log_Wald_pValue"] = 0.0
    dfc.loc[idx, "SE"] = SE
    dfc.loc[idx, "Wald"] = Wald
    dfc.loc[idx, "Wald_pValue"] = p
    dfc.loc[idx, "log_Wald_pValue"] = np.nan_to_num(lP, nan=0.0)

    df_final = dfc

    # ---------- 3  stack & write  --------------------------------------
    print("step3: write track_HMC, track_ES")
    df_final["signed_log_Wald_pValue"] = np.sign(df_final["logFC"]) * df_final["log_Wald_pValue"]
    df_final["HMC"] = df_final["m_corr"] * df_final["hmw"]
    # # --- NEW: log10+min-max scale HMC into [0,1] while preserving order ---
    # log_HMC      = np.log10(df_final['HMC'] + 1)
    # max_log_HMC  = log_HMC.max()
    # df_final['HMC'] = log_HMC / max_log_HMC

    # --- NEW: linear scale to [0,1] using the 99th percentile ---
    p_thr = df_final["HMC"].quantile(HMC_scale_pct)
    # clip the top 1% to p99, then divide so that p99 → 1
    df_final["HMC"] = df_final["HMC"].clip(upper=p_thr) / p_thr

    df_final[["chr", "start", "end", "HMC"]].to_csv(out_HMC, sep="\t", header=False, index=False)
    df_final[["chr", "start", "end", "signed_log_Wald_pValue"]].to_csv(
        out_ES, sep="\t", header=False, index=False
    )

    print(f"Saved outputs to {out_HMC} and {out_ES}")


def visualize_tracks(
    input_files,
    output_file,
    method="pyGenomeTracks",
    region=None,
    colors=None,
):
    if method == "pyGenomeTracks":
        visualize_with_pygenometracks(
            input_files=input_files,
            output_file=output_file,
            region=region,
            colors=colors,
        )
    elif method == "plotly":
        visualize_with_plotly(
            input_files=input_files,
            output_file=output_file,
            region=region,
            colors=colors,
        )
    else:
        raise ValueError(f"Unsupported visualization method: {method}")


def visualize_with_pygenometracks(input_files, output_file, region=None, colors=None):
    with tempfile.TemporaryDirectory() as temp_dir:
        config_file = os.path.join(temp_dir, "tracks.ini")
        with open(config_file, "w") as config:
            for idx, bedgraph_file in enumerate(input_files):
                track_id = f"track{idx + 1}"
                config.write(f"[{track_id}]\n")
                config.write(f"file = {bedgraph_file}\n")
                config.write(f"title = {os.path.basename(bedgraph_file)}\n")
                config.write("height = 4\n")
                config.write("type = line\n")

                # Use the provided or default color for each track
                try:
                    track_color = colors[idx]
                except (TypeError, IndexError):
                    print("Error: Insufficient colors provided for the number of input files.")
                    track_color = get_plotly_default_colors(1)[0]  # Fallback to default color

                config.write(f"color = {track_color}\n")
                config.write("\n")

        cmd = [
            "pyGenomeTracks",
            "--tracks",
            config_file,
            "--outFileName",
            output_file,
        ]

        if region:
            chrom, start, end = region
            region_str = f"{chrom}:{start}-{end}"
            cmd.extend(["--region", region_str])

        try:
            subprocess.run(cmd, check=True)
            print(f"Visualization saved to {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error running pyGenomeTracks: {e}")
            raise


def visualize_with_plotly(input_files, output_file, region=None, colors=None):
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    num_tracks = len(input_files)
    fig = make_subplots(
        rows=num_tracks,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.02,
    )

    dataframes = []  # List to hold dataframes if needed for y-axis scaling

    for idx, bedgraph_file in enumerate(input_files):
        df = pd.read_csv(bedgraph_file, sep="\t", header=None)
        df.columns = ["chrom", "start", "end", "value"]

        # Filter by region if specified
        if region:
            chrom, start, end = region
            df = df[(df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)]

        if df.empty:
            print(f"No data to plot for track {bedgraph_file}")
            continue

        df["position"] = (df["start"] + df["end"]) / 2
        dataframes.append(df)

        try:
            track_color = colors[idx]
        except (TypeError, IndexError):
            print("Error: Insufficient colors provided for the number of input files.")
            track_color = get_plotly_default_colors(1)[0]  # Fallback to default color

        fig.add_trace(
            go.Scattergl(
                x=df["position"],
                y=df["value"],
                mode="lines",
                name=os.path.basename(bedgraph_file),
                line=dict(width=1, color=track_color),  # Apply color here
            ),
            row=idx + 1,
            col=1,
        )

        # Add y-axis title for each subplot
        fig.update_yaxes(
            title_text=os.path.basename(bedgraph_file),
            row=idx + 1,
            col=1,
        )

    fig.update_layout(
        title="Genomic Track Visualization",
        xaxis_title="Genomic Position",
        showlegend=False,
        height=300 * num_tracks,
        hovermode="x unified",
        plot_bgcolor="white",
    )

    # Update x-axis line color and thickness
    fig.update_xaxes(showline=True, linewidth=1, linecolor="black")
    # Update y-axis line color and thickness for all subplots
    for idx in range(num_tracks):
        fig.update_yaxes(
            showline=True,
            linewidth=1,
            linecolor="black",
            row=idx + 1,
            col=1,
        )

    fig.write_html(output_file)
    print(f"Visualization saved to {output_file}")


def get_plotly_default_colors(num_colors):
    # Use Plotly's default color sequence (which has at least 10 colors by default)
    plotly_colors = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
    # Extend the color list if needed
    if num_colors > len(plotly_colors):
        # Repeat the color sequence to match the number of input files
        times = num_colors // len(plotly_colors) + 1
        plotly_colors = (plotly_colors * times)[:num_colors]
    return plotly_colors[:num_colors]


def find_significantly_different_regions(
    track_ES_file: str,
    track_HMC_file: str,
    output_dir: str,
    p_thresh: float = 0.05,
    binNum_thresh: int = 2,
    max_gap_bins: int = 0,
    chroms=None,
    chrom_sizes=None,
):
    """
    1) Read track_ES.bedgraph and track_HMC.bedgraph
    2) Merge on chr/start/end
    3) Keep bins with abs(ES) >= -log10(p_thresh) AND same-sign runs of length >= binNum_thresh.
       Then merge runs on the same chromosome if the empty gap between them is
       <= max_gap_bins bins (default: 1). Gaps do not contribute HMC values.   ### <<< CHANGED
    4) Write those bins (chr, start, end, ES, HMC) to '<output_dir>/pos_signif_bins.tsv'
       and '<output_dir>/neg_signif_bins.tsv'
    5) Merge adjacent/nearby bins into regions, average HMC across significant bins only,
       sort descending by avg_HMC, and write to
       '<output_dir>/signif_pos_regions.tsv' and '<output_dir>/signif_neg_regions.tsv'
    6) Plot rank vs log(HMC+1) for those regions, save
       '<output_dir>/HMC_rank_in_pos_regions.png' and '_neg_regions.png'
    """
    max_gap_bins = max(0, int(max_gap_bins))
    os.makedirs(output_dir, exist_ok=True)

    # 1) read
    df_es = pd.read_csv(
        track_ES_file,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "ES"],
    )
    print(df_es.shape)
    # print(df_es.head())
    df_c = pd.read_csv(
        track_HMC_file,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "HMC"],
    )
    print(df_c.shape)
    # print(df_c.head())
    # 2) merge
    df = pd.merge(df_es, df_c, on=["chr", "start", "end"])
    print(df.shape)
    # print(df.head())

    # 3) filter by ES
    thresh = -np.log10(p_thresh)
    print(f"thresh: {thresh}")
    # print(f'max_gap_binNum: {max_gap_binNum}')
    # split positive and negative extremes
    df_pos = df[df["ES"] >= thresh]
    df_neg = df[df["ES"] <= -thresh]
    # print(df_pos.shape)
    # print(df_pos.head(50))
    # print(df_neg.shape)
    # print(df_neg.head(50))

    # helper to merge runs on same chromosome
    def _merge_runs_allow_gaps(df_bins, max_gap_bins=max_gap_bins):  ### <<< NEW
        regs = []
        for chrom, grp in df_bins.groupby("chr"):
            grp = grp.sort_values("start")
            # infer bin size per chromosome (robust to a few irregular bins)
            if len(grp) == 0:
                continue
            bin_sizes = grp["end"].values - grp["start"].values
            # fallback to 1 if something weird
            bin_size = int(np.median(bin_sizes)) if len(bin_sizes) else 1
            bin_size = max(bin_size, 1)

            curr = None
            for _, row in grp.iterrows():
                if curr is None:
                    curr = dict(
                        chr=chrom,
                        start=int(row.start),
                        end=int(row.end),
                        wvals=[float(row.HMC)],
                        count=1,
                    )
                    continue

                gap_bp = int(row.start) - int(curr["end"])
                # allow merge if contiguous (gap==0) or small gap in bins
                if gap_bp >= 0 and gap_bp <= max_gap_bins * bin_size:  ### <<< NEW
                    # extend region but do NOT add any HMC for the gap
                    curr["end"] = int(row.end)
                    curr["wvals"].append(float(row.HMC))
                    curr["count"] += 1
                else:
                    # close previous run if long enough
                    if curr["count"] >= binNum_thresh:
                        regs.append(curr)
                    # start a new run
                    curr = dict(
                        chr=chrom,
                        start=int(row.start),
                        end=int(row.end),
                        wvals=[float(row.HMC)],
                        count=1,
                    )
            # flush last
            if curr and curr["count"] >= binNum_thresh:
                regs.append(curr)

        # build DataFrame
        if not regs:
            return pd.DataFrame(columns=["chr", "start", "end", "avg_HMC"])
        df_regs = pd.DataFrame(regs)
        df_regs["avg_HMC"] = df_regs["wvals"].apply(np.mean)
        return df_regs[["chr", "start", "end", "avg_HMC"]]

    # 4) write bins
    df_pos_sig = df_pos.sort_values(["chr", "start"])
    pos_sig_n = df_pos_sig.shape[0]
    pos_bins_out = os.path.join(output_dir, "pos_signif_bins.tsv")
    df_pos_sig[["chr", "start", "end", "ES", "HMC"]].to_csv(pos_bins_out, sep="\t", index=False)
    print(f"{pos_sig_n} pos ignificant bins written to {pos_bins_out}")
    df_neg_sig = df_neg.sort_values(["chr", "start"])
    neg_sig_n = df_neg_sig.shape[0]
    neg_bins_out = os.path.join(output_dir, "neg_signif_bins.tsv")
    df_neg_sig[["chr", "start", "end", "ES", "HMC"]].to_csv(neg_bins_out, sep="\t", index=False)
    print(f"{neg_sig_n} neg ignificant bins written to {neg_bins_out}")

    # 5) merge regions & write
    df_pos_regs = _merge_runs_allow_gaps(df_pos).sort_values("avg_HMC", ascending=False)
    # df_pos_regs['avg_HMC'] = df_pos_regs['avg_HMC'].round(2)
    df_pos_regs["avg_HMC"] = df_pos_regs["avg_HMC"]
    pos_n = df_pos_regs.shape[0]
    print(pos_n)
    print(df_pos_regs.shape)
    # print(df_pos_regs)

    pos_regs_out = os.path.join(output_dir, "signif_pos_regions.tsv")
    df_pos_regs.to_csv(pos_regs_out, sep="\t", index=False)
    print(f"{pos_n} merged pos regions written to {pos_regs_out}")

    df_neg_regs = _merge_runs_allow_gaps(df_neg).sort_values("avg_HMC", ascending=False)
    # df_neg_regs['avg_HMC'] = df_neg_regs['avg_HMC'].round(2)
    df_neg_regs["avg_HMC"] = df_neg_regs["avg_HMC"]
    neg_n = df_neg_regs.shape[0]
    print(neg_n)
    print(df_neg_regs.shape)
    # print(df_neg_regs)

    neg_regs_out = os.path.join(output_dir, "signif_neg_regions.tsv")
    df_neg_regs.to_csv(neg_regs_out, sep="\t", index=False)
    print(f"{neg_n} merged neg regions written to {neg_regs_out}")

    # --- NEW: also save BED files (no header), sorted by chr/start/end ---
    pos_bed_out = os.path.join(output_dir, "signif_pos_regions.bed")
    neg_bed_out = os.path.join(output_dir, "signif_neg_regions.bed")

    # ensure integer coordinates, then sort and save
    if not df_pos_regs.empty:
        # df_pos_bed = df_pos_regs[['chr','start','end','avg_HMC']].copy()
        df_pos_bed = df_pos_regs[["chr", "start", "end"]].copy()
        df_pos_bed[["start", "end"]] = df_pos_bed[["start", "end"]].astype(int)
        df_pos_bed = df_pos_bed.sort_values(["chr", "start", "end"])
        df_pos_bed.to_csv(pos_bed_out, sep="\t", header=False, index=False)
        print(f"BED written: {pos_bed_out}")

    if not df_neg_regs.empty:
        # df_neg_bed = df_neg_regs[['chr','start','end','avg_HMC']].copy()
        df_neg_bed = df_neg_regs[["chr", "start", "end"]].copy()
        df_neg_bed[["start", "end"]] = df_neg_bed[["start", "end"]].astype(int)
        df_neg_bed = df_neg_bed.sort_values(["chr", "start", "end"])
        df_neg_bed.to_csv(neg_bed_out, sep="\t", header=False, index=False)
        print(f"BED written: {neg_bed_out}")
    # --- end NEW ---

    # 6) plot rank vs log(w+1)
    df_pos_regs["rank"] = np.arange(1, len(df_pos_regs) + 1)
    plt.figure(figsize=(6, 4))
    # plt.plot(df_pos_regs['rank'], np.log(df_pos_regs['avg_HMC']+1), marker='.', linestyle='none',markersize=1, color='k')
    plt.plot(
        df_pos_regs["rank"],
        df_pos_regs["avg_HMC"],
        marker=".",
        linestyle="none",
        markersize=1,
        color="k",
    )
    plt.xlabel("Region rank")
    # plt.ylabel('log(HMC+1)')
    plt.ylabel("HMC")
    plt.title("Ranked region HMC")
    plt.tight_layout()
    fig_out = os.path.join(output_dir, "HMC_rank_in_pos_regions.png")
    plt.savefig(fig_out, dpi=600)
    plt.close()
    print(f"Rank plot of pos regions saved to {fig_out}")

    df_neg_regs["rank"] = np.arange(1, len(df_neg_regs) + 1)
    plt.figure(figsize=(6, 4))
    # plt.plot(df_neg_regs['rank'], np.log(df_neg_regs['avg_HMC']+1), marker='.', linestyle='none',markersize=1, color='k')
    plt.plot(
        df_neg_regs["rank"],
        df_neg_regs["avg_HMC"],
        marker=".",
        linestyle="none",
        markersize=1,
        color="k",
    )
    plt.xlabel("Region rank")
    # plt.ylabel('log(HMC+1)')
    plt.ylabel("HMC")
    plt.title("Ranked region HMC")
    plt.tight_layout()
    fig_out = os.path.join(output_dir, "HMC_rank_in_neg_regions.png")
    plt.savefig(fig_out, dpi=600)
    plt.close()
    print(f"Rank plot of neg regions saved to {fig_out}")
