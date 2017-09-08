#cp idea_label_matrix/data_matrix_rep2_cluster_epilabel_all.1.cdt data_matrix_rep2_cluster_epilabel_all.1.cdt

#python vlookup.py -t reads_count_matrix_sort.bed -m 4 -s data_matrix_rep2_cluster_epilabel_all.1.cdt -n 1 -o used_DNA_intervals.txt

#cat used_DNA_intervals.txt | sort -k1,1 -k2,2n | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' > used_DNA_intervals.bed

#python vlookup.py -t data_matrix_binarypattern.txt -m 1 -s used_DNA_intervals.bed -n 4 -o data_matrix_binarypattern_sorted.txt

#scp gzx103@biostar.psu.edu:/gpfs/home/gzx103/scratch/vision_clustering/used_DNA_intervals.bed ./
#scp gzx103@biostar.psu.edu:/gpfs/home/gzx103/scratch/vision_clustering/epi_cell_reoreder.txt ./
#scp gzx103@biostar.psu.edu:/gpfs/home/gzx103/scratch/vision_clustering/data_matrix_binarypattern_sorted.txt ./
#hsc,cmp,gmp,mep,cfu_e,ery,cfu_mk,meg,mono,neutrophil,b,nk,t_cd4,t_cd8,g1e,er4

for celltype_id_num in $(cat celltype_list.txt)
do
	echo $celltype_id_num
	### convert bigbed to bed file
	/Volumes/MAC_Data/data/labs/zhang_lab/roadmap_rnaseq_prediction/bigBedToBed 'ideasVision'$celltype_id_num'.bb' 'ideasVision'$celltype_id_num'.bed'
	cat 'ideasVision'$celltype_id_num'.bed' | sort -k1,1 -k2,2n > 'ideasVision'$celltype_id_num'.sort.bed'
	bedtools window -a used_DNA_intervals.bed -b 'ideasVision'$celltype_id_num'.sort.bed' -w 0 > 'used_DNA_intervals.'$celltype_id_num'.tmp.bed'
	### get input
	cat 'used_DNA_intervals.'$celltype_id_num'.tmp.bed' | awk -F '\t' -v OFS='\t' '{
		if ($6-$2>=0 && $7-$3<=0) print $1,$2,$3,$4,$7-$6,$8,$13,($7+$6-$3-$2)/2,$7-$6; 
		else if ($6-$2<0 && $7-$3>0) print $1,$2,$3,$4,$3-$2,$8,$13,($7+$6-$3-$2)/2,$7-$6;
		else if ($6-$2<0 && $7-$3<=0) print $1,$2,$3,$4,$7-$2,$8,$13,($7+$6-$3-$2)/2,$7-$6;
		else if ($6-$2>=0 && $7-$3>0) print $1,$2,$3,$4,$3-$6,$8,$13,($7+$6-$3-$2)/2,$7-$6
		}' > 'used_DNA_intervals.'$celltype_id_num'.bed'
	### get bed file with color
	python get_genomebrowser.py -a data_matrix_binarypattern_sorted.txt -i 'used_DNA_intervals.'$celltype_id_num'.bed' -o 'used_DNA_intervals.'$celltype_id_num'.colored.bed'
done

### reorder binary label
#cat data_matrix_binarypattern_sorted.txt | awk -F '\t' -v OFS='\t' '{print $2}' | awk -F '_' -v OFS='\t' '{print $1,$2,$4,$3,$8,$9,$10,$5,$6,$15,$16,$11,$12,$13,$14}' > data_matrix_binarypattern_sorted.colreordered.txt


paste used_DNA_intervals.9.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($10==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.lsk.9.colored.bed
paste used_DNA_intervals.3.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($11==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.cmp.3.colored.bed
paste used_DNA_intervals.11.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($12==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.mep.11.colored.bed
paste used_DNA_intervals.7.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($13==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.gmp.7.colored.bed
paste used_DNA_intervals.10.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($14==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.meg.10.colored.bed
paste used_DNA_intervals.12.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($15==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.mono.12.colored.bed
paste used_DNA_intervals.8.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($16==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.neutrophil.8.colored.bed
paste used_DNA_intervals.2.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($17==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.cfue.2.colored.bed
paste used_DNA_intervals.5.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($18==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.ery.5.colored.bed
paste used_DNA_intervals.6.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($19==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.g1e.6.colored.bed
paste used_DNA_intervals.4.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($20==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.er4.4.colored.bed
paste used_DNA_intervals.1.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($21==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.bcell.1.colored.bed
paste used_DNA_intervals.13.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($22==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.nkcell.13.colored.bed
paste used_DNA_intervals.14.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($23==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.tcd4.14.colored.bed
paste used_DNA_intervals.15.colored.bed data_matrix_binarypattern_sorted.colreordered.txt | awk -F '\t' -v OFS='\t' '{if ($24==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > used_DNA_intervals.tcd8.15.colored.bed

### add label names
time python add_label_name.py -i used_DNA_intervals.lsk.9.colored.bed -o used_DNA_intervals.lsk.9.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.cmp.3.colored.bed -o used_DNA_intervals.cmp.3.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.mep.11.colored.bed -o used_DNA_intervals.mep.11.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.gmp.7.colored.bed -o used_DNA_intervals.gmp.7.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.meg.10.colored.bed -o used_DNA_intervals.meg.10.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.mono.12.colored.bed -o used_DNA_intervals.mono.12.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.neutrophil.8.colored.bed -o used_DNA_intervals.neutrophil.8.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.cfue.2.colored.bed -o used_DNA_intervals.cfue.2.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.ery.5.colored.bed -o used_DNA_intervals.ery.5.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.g1e.6.colored.bed -o used_DNA_intervals.g1e.6.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.er4.4.colored.bed -o used_DNA_intervals.er4.4.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.bcell.1.colored.bed -o used_DNA_intervals.bcell.1.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.nkcell.13.colored.bed -o used_DNA_intervals.nkcell.13.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.tcd4.14.colored.bed -o used_DNA_intervals.tcd4.14.colored_named.bed
time python add_label_name.py -i used_DNA_intervals.tcd8.15.colored.bed -o used_DNA_intervals.tcd8.15.colored_named.bed

###format
#chr1	130318000	130318800	3	1000	.	130318000	130318800	218,218,221

### to bigbed
#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
#chmod +x ./bedToBigBed ./fetchChromSizes
#./fetchChromSizes mm10 > mm10.chrom.sizes

#ls used_DNA_intervals.*.*.colored.bed > bb_list.txt
#for bed in $(cat bb_list.txt)
#do
#	./bedToBigBed $bed mm10.chrom.sizes $bed'.bb'
#done


#track name="lsk" description="Item RGB demonstration" visibility=2 itemRgb="On"

#chr1:4776776-4827908
#cat data_matrix_binarypattern_sorted.txt | sort -k1,1 -k2,2n > homerTable3.peaks.filtered.sort.210000.bed

#for celltype_id_num in $(cat celltype_list.txt)
#do
#	echo $celltype_id_num
	### convert bigbed to bed file
	#/Volumes/MAC_Data/data/labs/zhang_lab/roadmap_rnaseq_prediction/bigBedToBed 'ideasVision'$celltype_id_num'.bb' ideasVision.bed
#	time bash split_idea_labels.sh
#	for id_num in $(cat label_list.txt)
#	do
#		echo $id_num
		### get total number at atac-seq peaks
#		cat 'ideasVision.'$id_num'.bed' | sort -k1,1 -k2,2n > 'ideasVision.'$id_num'.sort.bed'
#		bedtools intersect -a 'ideasVision.'$id_num'.sort.bed' -b homerTable3.peaks.filtered.sort.210000.bed -wa -u > 'ideasVision.'$id_num'.sort.interatac.bed'

#		python get_total_region.py -i 'ideasVision.'$id_num'.sort.interatac.bed' -o 'ideasVision.'$id_num'.totalnum.bed' -f 'state_'$id_num
#		cat 'ideasVision.'$id_num'.totalnum.bed' >> 'homerTable3.peaks.filtered.matrix.totalnum.txt'
#	done
#	mv homerTable3.peaks.filtered.matrix.txt 'homerTable3.peaks.filtered.matrix.'$celltype_id_num'.txt'
#	mv 'homerTable3.peaks.filtered.matrix.totalnum.txt' 'homerTable3.peaks.filtered.matrix.'$celltype_id_num'.totalnum.atac.txt'
#	mv 'homerTable3.peaks.filtered.matrix.totalnum.wg.txt' 'homerTable3.peaks.filtered.matrix.'$celltype_id_num'.totalnum.txt'
#	rm ideasVision.*.bed
#done





