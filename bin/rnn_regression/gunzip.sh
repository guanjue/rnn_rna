for i in {0..16}
do
	echo $i
	gunzip '/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn/matrix_for_rnn/gene_list.gmp.'$i'.tab.gz'
done


#gunzip gene_list.gmp.0.tab.gz