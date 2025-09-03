dest_dir=../annotations_test
mkdir -p ${dest_dir}
wget -t 0 -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -O ${dest_dir}/raw_hg19.chrom.sizes

# The chrom.sizes file downloaded from UCSC can be unsorted, so it’s convenient to sort the chromosomes here for downstream visualization.
grep -E '^chr([0-9]+|X|Y|M)\b' ${dest_dir}/raw_hg19.chrom.sizes | sort -V > ${dest_dir}/hg19.chrom.sizes
