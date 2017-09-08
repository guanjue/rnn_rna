while read -r celltype_id_num celltype; do
	echo $celltype
	### 
	for state_id_num in $(cat /Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/lists/label_list.txt)
	do
		echo $state_id_num
		time /Users/gzx103/Documents/apps/deepTools/bin/plotHeatmap -m 'gene_list.'$celltype'.'$state_id_num'.tab.gz' --sortUsing region_length --kmeans 1 --sortRegions no --colorMap binary --outFileName 'gene_list.'$celltype'.'$state_id_num'.tab.pdf' -max 1 -min 0
	done
done < /Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/lists/cell_type_label.txt
