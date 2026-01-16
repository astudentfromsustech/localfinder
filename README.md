# localfinder

localfinder – calculate weighted local correlation (HMC) and enrichment significance (ES) between two genomic tracks, optionally discover significantly different regions, and visualize results. 

## Installation Requirements

Before installing and using `localfinder`, please ensure that the following external tools are installed on your `$PATH`:

- **bedtools**: Used for genomic interval operations.
  - Installation: [https://bedtools.readthedocs.io/en/latest/content/installation.html](https://bedtools.readthedocs.io/en/latest/content/installation.html)
  - conda install -c bioconda -c conda-forge bedtools 
  - mamba install -c bioconda -c conda-forge bedtools
- **bigWigToBedGraph (UCSC utility)**: Used for converting BigWig files to BedGraph format.
  - Download: [http://hgdownload.soe.ucsc.edu/admin/exe/](http://hgdownload.soe.ucsc.edu/admin/exe/)
  - conda install -c bioconda -c conda-forge ucsc-bigwigtobedgraph
  - mamba install -c bioconda -c conda-forge ucsc-bigwigtobedgraph
- **samtools**: Used for processing SAM/BAM files.
  - Installation: [http://www.htslib.org/download/](http://www.htslib.org/download/)
  - conda install -c bioconda -c conda-forge samtools
  - mamba install -c bioconda -c conda-forge samtools

These tools are required for processing genomic data and must be installed separately.

## Installation

Install `localfinder` using `pip`:

```bash
pip install localfinder
```

## Usage
There are 5 subcommands (bin, calc, findreg, viz, pipeline) in localfinder, and you can check it using:
```bash
localfinder -h
```

### bin
```bash
localfinder bin -h
```

```
usage: localfinder bin [-h] --input_files FILE [FILE ...] --output_dir DIR [--bin_size BP] --chrom_sizes FILE
                       [--chroms CHR [CHR ...]] [--threads N]

Bin genomic tracks into fixed-size bins and output BedGraph format.

options:
  -h, --help                    show this help message and exit

I/O:
  --input_files FILE [FILE ...]
                                Input files in BigWig/BedGraph/BAM format. (default: None)
  --output_dir DIR              Output directory for binned data. (default: None)

Genomics:
  --bin_size BP                 Size of each bin. (default: 200)
  --chrom_sizes FILE            Path to the chromosome sizes file. (default: None)
  --chroms CHR [CHR ...]        'all' or specific chromosomes (e.g. chr1 chr2). (default: ['all'])

Performance:
  --threads N, -t N             Number of worker processes to run in parallel. (default: 1)

Usage Examples
--------------
1) Example 1:
   localfinder bin \
     --input_files track1.bw track2.bw \
     --output_dir ./binned_tracks \
     --bin_size 200 \
     --chrom_sizes mm10.chrom.sizes \
     --chroms chr1 chr2

2) Example 2:
   localfinder bin \
     --input_files track1.bigwig track2.bigwig \
     --output_dir ./binned_tracks \
     --bin_size 200 \
     --chrom_sizes hg19.chrom.sizes \
     --chroms all \
     --threads 4
```

### calc
```bash
localfinder calc -h
```
```
usage: localfinder calc [-h] --track1 BEDGRAPH --track2 BEDGRAPH --output_dir DIR [--method {locP_and_ES,locS_and_ES}]
                        [--FDR] [--binNum_window BINNUM_WINDOW] [--step STEP] [--percentile PERCENTILE]
                        [--percentile_mode {all,nonzero}] [--binNum_peak BINNUM_PEAK] [--FC_thresh FC_THRESH]
                        [--norm_method {scale,cpm,rpkm}] [--HMC_scale_pct HMC_SCALE_PCT] --chrom_sizes FILE
                        [--chroms CHR [CHR ...]] [--threads N]

Calculate hormonic-mean weighted local correlation (HMC) and enrichment significance (ES) between two BedGraph tracks.

options:
  -h, --help                    show this help message and exit

I/O:
  --track1 BEDGRAPH             First input BedGraph file. (default: None)
  --track2 BEDGRAPH             Second input BedGraph file. (default: None)
  --output_dir DIR              Output directory for results. (default: None)

Method:
  --method {locP_and_ES,locS_and_ES}
                                P=Pearson, S=Spearman for local correlation. (default: locP_and_ES)
  --FDR                         Use Benjamini–Hochberg q-values instead of raw P-values. (default: False)
  --binNum_window BINNUM_WINDOW
                                Number of bins in the sliding window. (default: 11)
  --step STEP                   Step size for the sliding window. (default: 1)
  --percentile PERCENTILE       Percentile for floor correction of low-coverage bins. (default: 90)
  --percentile_mode {all,nonzero}
                                Use all bins or only non-zero bins for percentile. (default: all)
  --binNum_peak BINNUM_PEAK     Number of bins in the ES peak window. (default: 3)
  --FC_thresh FC_THRESH         Fold-change threshold used as log base in enrichment. (default: 1.5)
  --norm_method {scale,cpm,rpkm}
                                Normalisation strategy. (default: rpkm)
  --HMC_scale_pct HMC_SCALE_PCT
                                Quantile used to clip+rescale HMC into [0,1]. (default: 0.9995)

Genomics:
  --chrom_sizes FILE            Path to the chromosome sizes file. (default: None)
  --chroms CHR [CHR ...]        Chromosomes to process (e.g. chr1 chr2) or 'all'. (default: ['all'])

Performance:
  --threads N, -t N             Worker processes per chromosome. (default: 1)

Usage Examples
--------------
1) Example 1:
   localfinder calc \
     --track1 track1.bedgraph --track2 track2.bedgraph \
     --output_dir ./results \
     --method locP_and_ES --FDR \
     --binNum_window 11 --step 1 \
     --percentile 90 --percentile_mode all \
     --binNum_peak 3 --FC_thresh 1.5 \
     --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2 \
     --threads 4

2) Example 2:
   localfinder calc \
     --track1 track1.bedgraph --track2 track2.bedgraph \
     --output_dir ./results \
     --percentile 99 --binNum_peak 2 \
     --chrom_sizes hg19.chrom.sizes
```

### findreg
```bash
localfinder findreg -h
```
```
usage: localfinder findreg [-h] --track_E BEDGRAPH --track_C BEDGRAPH --output_dir DIR [--p_thresh P_THRESH]
                           [--binNum_thresh BINNUM_THRESH] [--max_gap_bins MAX_GAP_BINS] [--chroms CHR [CHR ...]]
                           --chrom_sizes FILE

Find SDRs by ES and further classify SDRs by their HMC

options:
  -h, --help                    show this help message and exit

I/O:
  --track_E BEDGRAPH            ES track (BedGraph). (default: None)
  --track_C BEDGRAPH            HMC track (BedGraph). (default: None)
  --output_dir DIR

Thresholds:
  --p_thresh P_THRESH           Bin-level P-value threshold. (default: 0.05)
  --binNum_thresh BINNUM_THRESH
                                Minimum consecutive significant bins per region. (default: 2)
  --max_gap_bins MAX_GAP_BINS   Merge two significant runs if the intervening gap ≤ this number of bins. (default: 0)

Genomics:
  --chroms CHR [CHR ...]
  --chrom_sizes FILE

Usage Examples
--------------
1) Example 1:
   localfinder findreg \
     --track_E track_ES.bedgraph --track_C track_HMC.bedgraph \
     --output_dir ./findreg_out \
     --p_thresh 0.05 --binNum_thresh 2 --max_gap_bins 0 \
     --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

2) Example 2:
   localfinder findreg \
     --track_E track_ES.bedgraph --track_C track_HMC.bedgraph \
     --output_dir ./findreg_out \
     --chrom_sizes hg19.chrom.sizes
```

### pipeline
```bash
localfinder pipeline -h
```
```
usage: localfinder pipeline [-h] --input_files FILE [FILE ...] --output_dir DIR --chrom_sizes FILE [--bin_size BP]
                            [--chroms CHR [CHR ...]] [--method {locP_and_ES,locS_and_ES}] [--FDR]
                            [--binNum_window BINNUM_WINDOW] [--binNum_peak BINNUM_PEAK] [--step STEP]
                            [--percentile PERCENTILE] [--percentile_mode {all,nonzero}] [--FC_thresh FC_THRESH]
                            [--norm_method {scale,cpm,rpkm}] [--HMC_scale_pct HMC_SCALE_PCT] [--p_thresh P_THRESH]
                            [--binNum_thresh BINNUM_THRESH] [--max_gap_bins MAX_GAP_BINS] [--threads N]

Run localfinder sequentially: bin → calc → findreg.

options:
  -h, --help                    show this help message and exit

I/O:
  --input_files FILE [FILE ...]
                                Input files in BigWig/BedGraph format. (default: None)
  --output_dir DIR              Output directory for the pipeline results. (default: None)

Genomics:
  --chrom_sizes FILE            Path to the chromosome sizes file. (default: None)
  --bin_size BP                 Size of each bin. (default: 200)
  --chroms CHR [CHR ...]        Chromosomes to process or "all". (default: ['all'])

calc options:
  --method {locP_and_ES,locS_and_ES}
                                P=Pearson, S=Spearman for local correlation. (default: locP_and_ES)
  --FDR                         Use Benjamini–Hochberg q-values. (default: False)
  --binNum_window BINNUM_WINDOW
                                Number of bins in the sliding window. (default: 11)
  --binNum_peak BINNUM_PEAK     Number of bins of the peak for ES. (default: 3)
  --step STEP                   Step size for the sliding window. (default: 1)
  --percentile PERCENTILE       Percentile for floor correction of low-coverage bins. (default: 90)
  --percentile_mode {all,nonzero}
                                Use all bins or only non-zero bins for percentile. (default: all)
  --FC_thresh FC_THRESH         Fold-change threshold used as log base in enrichment. (default: 1.5)
  --norm_method {scale,cpm,rpkm}
                                Normalisation strategy. (default: rpkm)
  --HMC_scale_pct HMC_SCALE_PCT
                                Quantile used to clip+rescale HMC into [0,1]. (default: 0.9995)

findreg options:
  --p_thresh P_THRESH           P-value threshold for merging significant bins into regions. (default: 0.05)
  --binNum_thresh BINNUM_THRESH
                                Minimum consecutive significant bins per region. (default: 2)
  --max_gap_bins MAX_GAP_BINS   Merge two significant runs if the intervening gap ≤ this number of bins. (default: 0)

Performance:
  --threads N, -t N             Number of worker processes for bin & calc. (default: 1)

Usage Examples
--------------
1) Example 1:
   localfinder pipeline \
     --input_files track1.bedgraph track2.bedgraph \
     --output_dir ./results \
     --chrom_sizes hg19.chrom.sizes \
     --bin_size 200 \
     --method locP_and_ES --FDR \
     --binNum_window 11 --binNum_peak 3 --step 1 \
     --percentile 90 --percentile_mode all \
     --FC_thresh 1.5 --norm_method rpkm --HMC_scale_pct 0.9995 \
     --p_thresh 0.05 --binNum_thresh 2 --max_gap_bins 0 \
     --chroms chr1 chr2 --threads 4

2) Example 2:
   localfinder pipeline \
     --input_files track1.bigwig track2.bigwig \
     --output_dir ./results \
     --chrom_sizes hg19.chrom.sizes
```

### viz
```bash
localfinder viz -h
```
```
usage: localfinder viz [-h] --input_files BEDGRAPH [BEDGRAPH ...] --output_file FILE --method {pyGenomeTracks,plotly}
                       [--region CHROM START END] [--colors COLORS [COLORS ...]]

Visualize genomic tracks.

options:
  -h, --help                    show this help message and exit

I/O:
  --input_files BEDGRAPH [BEDGRAPH ...]
                                Input BedGraph files to visualize. (default: None)
  --output_file FILE            Output visualization file (e.g. PNG, HTML). (default: None)

Options:
  --method {pyGenomeTracks,plotly}
                                Visualization method to use. (default: None)
  --region CHROM START END      Region to visualize, e.g. chr20 1000000 2000000. (default: None)
  --colors COLORS [COLORS ...]  Colors for the tracks (optional). (default: None)

Usage Examples
--------------
1) Example 1:
   localfinder viz \
     --input_files track1.bedgraph track2.bedgraph \
     --output_file output.html \
     --method plotly \
     --region chr1 1000000 2000000 \
     --colors blue red

2) Example 2:
   localfinder viz \
     --input_files track1.bedgraph track2.bedgraph \
     --output_file output.png \
     --method pyGenomeTracks \
     --region chr1 1000000 2000000
```

## Run an example step by step
Create a conda env called localfinder and enter this conda environment
```bash
conda create -n localfinder
conda activate  localfinder
```

Install external tools and localfinder
```bash
conda install -c conda-forge -c bioconda samtools bedtools ucsc-bigwigtobedgraph
pip install localfinder
```

Download the souce code of [localfinder](https://github.com/astudentfromsustech/localfinder)  
```bash
git clone https://github.com/astudentfromsustech/localfinder.git
```

Run the examples under localfinder/tests/ (scripts have been preprared in tests folder)  
