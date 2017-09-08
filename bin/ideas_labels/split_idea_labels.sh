rm 'homerTable3.peaks.filtered.matrix.totalnum.wg.txt'

for state_id_num in $(cat label_list.txt)
do
	echo $state_id_num
	cat ideasVision.bed | awk -F '\t' -v OFS='\t' -v var="$state_id_num" '{if ($4==var) print $1,$2,$3,$4,$5,$6}' > 'ideasVision.'$state_id_num'.bed'
	### get total number of region at whole genome level
	python get_total_region.py -i 'ideasVision.'$state_id_num'.bed' -o 'ideasVision.'$state_id_num'.totalnum.wg.bed' -f 'state_'$state_id_num
	cat 'ideasVision.'$state_id_num'.totalnum.wg.bed' >> 'homerTable3.peaks.filtered.matrix.totalnum.wg.txt'
done

