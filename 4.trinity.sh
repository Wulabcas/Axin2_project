trinity=path
for ((i=1;i<=6;i++));
do
$trinity/abundance_estimates_to_matrix.pl \
	 --est_method RSEM \
	 --cross_sample_norm TMM \
	 --gene_trans_map none \
	 --quant_files genes.quant_files.txt \
	 --out_prefix rsem_matrix_M${i}
done

