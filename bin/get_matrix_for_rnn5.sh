cd /Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/matrix_for_rnn/
### get gene list
#cat /Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/rsem/LSK.rsem.1007.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4,$3-$2,$6}' | sort -k5,5rn > gene_list.bed

while read -r celltype_id_num celltype; do
	echo $celltype
	### 
	for state_id_num in $(cat /Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/lists/label_list.txt)
	do
		echo $state_id_num
		time /Users/gzx103/Documents/apps/deepTools/bin/computeMatrix scale-regions --scoreFile '/Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/ideas_labels/ideasVision.'$celltype'.'$state_id_num'.bw' --regionsFile gene_list.bed --beforeRegionStartLength 100000 --afterRegionStartLength 100000 --regionBodyLength 100000 --binSize 2000 -out 'gene_list.'$celltype'.'$state_id_num'.tab.gz'
		time /Users/gzx103/Documents/apps/deepTools/bin/plotHeatmap -m 'gene_list.'$celltype'.'$state_id_num'.tab.gz' --sortUsing region_length --kmeans 1 --sortRegions no --colorMap binary --outFileName 'gene_list.'$celltype'.'$state_id_num'.tab.pdf'
	done
done < /Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/lists/cell_type_label5.txt
