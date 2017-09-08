cd /Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/rsem/

paste CFU-E.rsem.1046.bed CFU-E.rsem.1087.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16,($4+$16)/2,$3-$2,$6}' | sort -k6,6rn > cfuE_ad.rsem.txt
paste CFU-Mk.rsem.1049.bed CFU-Mk.rsem.1050.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16,($4+$16)/2,$3-$2,$6}' | sort -k6,6rn > cfuMK_ad.rsem.txt
paste CMP.rsem.1009.bed CMP.rsem.1010.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16,($4+$16)/2,$3-$2,$6}' | sort -k6,6rn > cmp.rsem.txt
paste ER4.rsem.985.bed ER4.rsem.986.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16,($4+$16)/2,$3-$2,$6}' | sort -k6,6rn > er4.rsem.txt
paste G1E.rsem.983.bed G1E.rsem.984.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16,($4+$16)/2,$3-$2,$6}' | sort -k6,6rn > g1e.rsem.txt
paste GMP.rsem.1017.bed GMP.rsem.1018.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16,($4+$16)/2,$3-$2,$6}' | sort -k6,6rn > gmp.rsem.txt
paste LSK.rsem.1007.bed LSK.rsem.1063.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16,($4+$16)/2,$3-$2,$6}' | sort -k6,6rn > lsk_bm.rsem.txt
paste MEP.rsem.1019.bed MEP.rsem.1064.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16,($4+$16)/2,$3-$2,$6}' | sort -k6,6rn > mep.rsem.txt
paste Mono.rsem.1158.bed Mono.rsem.1208.bed Mono.rsem.1233.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16"_"$28,($4+$16+$28)/3,$3-$2,$6}' | sort -k6,6rn > mono_bm.rsem.txt
paste ery.rsem.1047.bed ery.rsem.1088.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16,($4+$16)/2,$3-$2,$6}' | sort -k6,6rn > ery.rsem.txt
paste megs.rsem.1051.bed megs.rsem.1052.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16,($4+$16)/2,$3-$2,$6}' | sort -k6,6rn > meg.rsem.txt
paste neutrophil.rsem.1156.bed neutrophil.rsem.1157.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=2000) print $1,$2,$3,$4"_"$16,($4+$16)/2,$3-$2,$6}' | sort -k6,6rn > gra_bm.rsem.txt

