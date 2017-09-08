
for celltype_id_num in $(cat celltype_list.txt)
do
	cat data_matrix_binarypattern_sorted.txt | sort -k3,3 -k4,4n | awk -F '\t' -v OFS='\t' '{print $3,$4,$5,$6,$2,$1}'> homerTable3.peaks.filtered.sort.210000.bed
	echo $celltype_id_num
	### convert bigbed to bed file
	/Volumes/MAC_Data/data/labs/zhang_lab/roadmap_rnaseq_prediction/bigBedToBed 'ideasVision'$celltype_id_num'.bb' ideasVision.bed
	time bash split_idea_labels.sh
	for id_num in $(cat label_list.txt)
	do
		echo $id_num
		### get total number at atac-seq peaks
		cat 'ideasVision.'$id_num'.bed' | sort -k1,1 -k2,2n > 'ideasVision.'$id_num'.sort.bed'
		bedtools intersect -b 'ideasVision.'$id_num'.sort.bed' -a homerTable3.peaks.filtered.sort.210000.bed -wa -u > 'ideasVision.'$id_num'.sort.interatac.bed'
		#-wa -u
		python get_total_region.py -i 'ideasVision.'$id_num'.sort.interatac.bed' -o 'ideasVision.'$id_num'.totalnum.bed' -f 'state_'$id_num
		cat 'ideasVision.'$id_num'.totalnum.bed' >> 'homerTable3.peaks.filtered.matrix.totalnum.txt'
	done
	mv homerTable3.peaks.filtered.matrix.txt 'homerTable3.peaks.filtered.matrix.'$celltype_id_num'.txt'
	mv 'homerTable3.peaks.filtered.matrix.totalnum.txt' 'homerTable3.peaks.filtered.matrix.'$celltype_id_num'.totalnum.atac.txt'
	mv 'homerTable3.peaks.filtered.matrix.totalnum.wg.txt' 'homerTable3.peaks.filtered.matrix.'$celltype_id_num'.totalnum.txt'
	rm *bed
done
