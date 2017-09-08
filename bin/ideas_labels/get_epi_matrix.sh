for celltype_id_num in $(cat celltype_list.txt)
do
	echo $celltype_id_num
	### convert bigbed to bed file
	/Volumes/MAC_Data/data/labs/zhang_lab/roadmap_rnaseq_prediction/bigBedToBed 'ideasVision'$celltype_id_num'.bb' ideasVision.bed
	### get the label_list.txt
	#cat ideasVision.bed | awk '{print $4}' | sort -u > label_list.txt
	### split idea label into seperate bed files
	time bash split_idea_labels.sh 
	### total number of region intersect with atac-seq peak
	time bash total_region_atacseq.sh
	### intersect
	time bash intersect_labels.sh
	mv homerTable3.peaks.filtered.matrix.txt 'homerTable3.peaks.filtered.matrix.'$celltype_id_num'.txt'
	mv 'homerTable3.peaks.filtered.matrix.totalnum.txt' 'homerTable3.peaks.filtered.matrix.'$celltype_id_num'.totalnum.atac.txt'
	mv 'homerTable3.peaks.filtered.matrix.totalnum.wg.txt' 'homerTable3.peaks.filtered.matrix.'$celltype_id_num'.totalnum.txt'
	rm *bed
done
