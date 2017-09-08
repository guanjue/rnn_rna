cd /Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/atac_seq_pk

cat atac_21k_binary.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' > atac_21k_binary_split.txt

tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$9}' > atac_21k_binary_split.cfuE_ad.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$11}' > atac_21k_binary_split.cfuMK_ad.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$6}' > atac_21k_binary_split.cmp.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$20}' > atac_21k_binary_split.er4.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$10}' > atac_21k_binary_split.ery_ad.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$19}' > atac_21k_binary_split.g1e_ad.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$7}' > atac_21k_binary_split.gmp.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$14}' > atac_21k_binary_split.gra_bm.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5}' > atac_21k_binary_split.lsk_bm.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$12}' > atac_21k_binary_split.meg.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$8}' > atac_21k_binary_split.mep.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$13}' > atac_21k_binary_split.mono_bm.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$15}' > atac_21k_binary_split.b_spl.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$16}' > atac_21k_binary_split.nk_spl.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$17}' > atac_21k_binary_split.t_cd4_spl.txt
tail -n+2 atac_21k_binary_split.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$18}' > atac_21k_binary_split.t_cd8_spl.txt


