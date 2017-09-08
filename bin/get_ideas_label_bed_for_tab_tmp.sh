cd /Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/ideas_labels/

while read -r celltype_id_num celltype; do
	echo $celltype_id_num 
	echo $celltype
	### convert bigbed to bed file
	/Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/bin/bigBedToBed 'ideasVision'$celltype_id_num'.bb' 'ideasVision.'$celltype'.bed'
	### 
	rm $celltype'.totalregion.wg.txt'
	for state_id_num in $(cat /Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/lists/label_list.txt)
	do
		echo $state_id_num
		### extract ideas_label
		cat 'ideasVision.'$celltype'.bed' | awk -F '\t' -v OFS='\t' -v var="$state_id_num" '{if ($4==var) print $1,$2,$3,$4,1,$6; else print $1,$2,$3,$4,0,$6}' > 'ideasVision.'$celltype'.'$state_id_num'.bed'
		### get total number of region at whole genome level
		python get_total_region.py -i 'ideasVision.'$celltype'.'$state_id_num'.bed' -o 'ideasVision.'$celltype'.'$state_id_num'.totalnum.wg.bed' -f 'state_'$state_id_num
		cat 'ideasVision.'$celltype'.'$state_id_num'.totalnum.wg.bed' >> $celltype'.totalregion.wg.txt'
		### convert bed file to bedgraph then to bigwig
		cat 'ideasVision.'$celltype'.'$state_id_num'.bed' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$5}' | sort -k1,1 -k2,2n > 'ideasVision.'$celltype'.'$state_id_num'.bedgraph'
		/Users/gzx103/Documents/apps/ucsc/bedGraphToBigWig 'ideasVision.'$celltype'.'$state_id_num'.bedgraph' mm10.chrom.sizes 'ideasVision.'$celltype'.'$state_id_num'.bw'
	done
done < /Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/lists/cell_type_label_miss.txt

