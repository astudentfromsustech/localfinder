# conda create -n tools_igvtools
# conda activate tools_igvtools
# conda activate -c bioconda -c conda-forge igvtools=2.16.2
genome="hg19"
src_dir="../test_data_chr20"

# genome="mm10"
# src_dir="../test_ChIPseq_data_chr1_from_NatureGenetics2011"


for bg in "$src_dir"/*.bedgraph; do
  base="${bg%.bedgraph}"
  echo ${base}
  tdf=${base}.tdf

  igvtools toTDF ${bg} ${tdf} ${genome}
done
