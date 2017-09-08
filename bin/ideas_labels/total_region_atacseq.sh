tail -n+2 /Volumes/MAC_Data/data/labs/zhang_lab/vision_project/clustering/data/supervised_clustering/homerTable3.peaks.filtered.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k4,4 > homerTable3.peaks.filtered.bed 
cat homerTable3.peaks.filtered.bed | sort -k1,1 -k2,2n > homerTable3.peaks.filtered.sort.bed
for id_num in $(cat label_list.txt)
do
	echo $id_num
	### get total number at atac-seq peaks
	cat 'ideasVision.'$id_num'.bed' | sort -k1,1 -k2,2n > 'ideasVision.'$id_num'.sort.bed'
	bedtools intersect -a 'ideasVision.'$id_num'.sort.bed' -b homerTable3.peaks.filtered.sort.bed -wa -u > 'ideasVision.'$id_num'.sort.interatac.bed'

	python get_total_region.py -i 'ideasVision.'$id_num'.sort.interatac.bed' -o 'ideasVision.'$id_num'.totalnum.bed' -f 'state_'$id_num
	cat 'ideasVision.'$id_num'.totalnum.bed' >> 'homerTable3.peaks.filtered.matrix.totalnum.txt'
done