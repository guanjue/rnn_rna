tail -n+2 /Volumes/MAC_Data/data/labs/zhang_lab/vision_project/clustering/data/supervised_clustering/homerTable3.peaks.filtered.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k4,4 > homerTable3.peaks.filtered.bed 
cp homerTable3.peaks.filtered.bed homerTable3.peaks.filtered.matrix.txt
for state_id_num in $(cat label_list.txt)
do
	echo $state_id_num
	bedtools window -a homerTable3.peaks.filtered.bed -b 'ideasVision.'$state_id_num'.bed' -w 0 -c > homerTable3.ideasVision.tmp.bed
	cat homerTable3.ideasVision.tmp.bed | awk -F '\t' -v OFS='\t' '{print $5}' > tmp.labelnum.txt
	paste homerTable3.peaks.filtered.matrix.txt tmp.labelnum.txt > homerTable3.peaks.filtered.matrix.txt.tmp && mv homerTable3.peaks.filtered.matrix.txt.tmp homerTable3.peaks.filtered.matrix.txt
	rm tmp.labelnum.txt
	rm homerTable3.ideasVision.tmp.bed
done

