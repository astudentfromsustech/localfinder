# conda create -n tools_igvtools
# conda activate tools_igvtools
# conda install -c bioconda -c conda-forge igvtools=2.16.2

# genome="hg19"
# src_dir="../test_data_chr20"

genome="mm9"
# src_dir="../result_test_CTCF_ChIP-seq_data_from_Cell2017_using_pipeline/correlation_enrichment"
src_dir="../result_test_CTCF_ChIP-seq_data_from_Cell2017_step_by_step/correlation_enrichment"


for bg in "$src_dir"/*.bedgraph; do
  base="${bg%.bedgraph}"
  echo ${base}
  tdf=${base}.tdf

  igvtools toTDF ${bg} ${tdf} ${genome}
done
