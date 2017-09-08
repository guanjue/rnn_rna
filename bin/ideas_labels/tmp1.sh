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

