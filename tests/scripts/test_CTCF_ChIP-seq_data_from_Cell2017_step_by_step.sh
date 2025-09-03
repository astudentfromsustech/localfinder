#!/bin/bash
export PYTHONUNBUFFERED=1

chrom_size="../annotations/mm9.chrom.sizes"
data_dir="../test_CTCF_ChIP-seq_data_from_Cell2017"
mkdir -p ${data_dir}
wget -t 0 -c ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2609nnn/GSM2609185/suppl/GSM2609185_CTCF_ChIP-seq_CTCF-AID_untreated_rep1_ENC124_tagDensity.bw -O ${data_dir}/CTCF_untreated.bw
wget -t 0 -c ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2609nnn/GSM2609186/suppl/GSM2609186_CTCF_ChIP-seq_CTCF-AID_auxin2days_rep1_ENC125_tagDensity.bw -O ${data_dir}/CTCF_treated.bw

result_dir="../result_test_CTCF_ChIP-seq_data_from_Cell2017_step_by_step"
mkdir -p ${result_dir}
track_1="${data_dir}/CTCF_untreated.bw"
track_2="${data_dir}/CTCF_treated.bw"

# # After visulaization of the two tracks, we can find that most peaks are sharp peaks and most bins are with low values, so we choose large percentile around 99 (after testing, we choose 99.4) and binNum_peak as 3
bin_dir="${result_dir}/bin_tracks"
mkdir -p ${bin_dir}
localfinder bin --input_files ${track_1} ${track_2} --output_dir ${bin_dir} --bin_size 200 --chrom_sizes ${chrom_size} --chroms chr1 chr2

correlation_enrichment_dir="${result_dir}/correlation_enrichment"
# correlation_enrichment_dir="${result_dir}/correlation_enrichment"
mkdir -p ${correlation_enrichment_dir}
base1=$(basename "$track_1" .bw)
bin_track_1="${result_dir}/bin_tracks/${base1}.binSize200.bedgraph"
# echo ${bin_track_1}
base2=$(basename "$track_2" .bw)
bin_track_2="${result_dir}/bin_tracks/${base2}.binSize200.bedgraph"
# echo ${bin_track_2}
localfinder calc --track1 ${bin_track_1} --track2 ${bin_track_2} --output_dir ${correlation_enrichment_dir} \
--binNum_peak 3 --percentile 99.4 --threads 4 --norm_method rpkm --FDR --HMC_scale_pct 0.9995 --chrom_sizes ${chrom_size} \
--chroms chr1 chr2

significant_regions_dir="${result_dir}/significant_regions"
mkdir -p ${significant_regions_dir}
localfinder findreg --track_E ${correlation_enrichment_dir}/track_ES.bedgraph --track_C ${correlation_enrichment_dir}/track_HMC.bedgraph \
--output_dir ${significant_regions_dir} --binNum_thresh 2 --chrom_size ${chrom_size} --chroms chr1 chr2
