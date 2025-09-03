#!/bin/bash
export PYTHONUNBUFFERED=1

chrom_sizes="../annotations/mm9.chrom.sizes"
data_dir="../test_CTCF_ChIP-seq_data_from_Cell2017"
mkdir -p ${data_dir}
wget -t 0 -c ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2609nnn/GSM2609185/suppl/GSM2609185_CTCF_ChIP-seq_CTCF-AID_untreated_rep1_ENC124_tagDensity.bw -O ${data_dir}/CTCF_untreated.bw
wget -t 0 -c ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2609nnn/GSM2609186/suppl/GSM2609186_CTCF_ChIP-seq_CTCF-AID_auxin2days_rep1_ENC125_tagDensity.bw -O ${data_dir}/CTCF_treated.bw

result_dir="../result_test_CTCF_ChIP-seq_data_from_Cell2017_using_pipeline"
mkdir -p ${result_dir}

track_1="${data_dir}/CTCF_untreated.bw"
track_2="${data_dir}/CTCF_treated.bw"


# After visulaization of the two tracks, we can find that most peaks are sharp peaks and most bins are with low values, so we choose large percentile around 99 (after testing, we choose 99.4) and binNum_peak as 3
localfinder pipeline \
    --input_files ${track_1} ${track_2} \
    --output_dir  ${result_dir} \
    --chrom_sizes $chrom_sizes \
    --bin_size    200 \
    --method      locP_and_ES \
    --binNum_window 11 --binNum_peak 3 --step 1 \
    --percentile 99.4 --percentile_mode all --FC_thresh 1.5 \
    --p_thresh 0.05 --binNum_thresh 2 \
    --threads 4 --norm_method rpkm \
    --chroms chr1 chr2
