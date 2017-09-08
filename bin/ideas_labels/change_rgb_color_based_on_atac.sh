while read -r ideas_id celltype bam_a bam_b; do
	echo $ideas_id 
	echo $celltype
	echo $bam_a
	echo $bam_b
	time python /gpfs/home/gzx103/scratch/vision_clustering/ideas_with_atac_color.py -i 'ideas_bed_split.'$ideas_id'.bed' -a 'ideas_bed_split.'$bam_a'_result3.bam.sort.bam.bed' -b 'ideas_bed_split.'$bam_b'_result3.bam.sort.bam.bed' -o $celltype'_result3.ideas_bed_split.'$ideas_id'.bed' -s 3
done < cell_type_label_atacsample.txt


curl http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
curl http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
chmod +x ./bedToBigBed ./fetchChromSizes
./fetchChromSizes mm10 > mm10.chrom.sizes


/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb
/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb
/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb
/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb
/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb
/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb
/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb
/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb
/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb
/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb
/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb
/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed b_spl_result3.ideas_bed_split.int.1.bed /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes b_spl_result3.ideas_bed_split.int.1.bb


wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
chmod +x ./bedToBigBed ./fetchChromSizes
./fetchChromSizes mm10 > mm10.chrom.sizes


for filename in $(cat od_bed_list.txt)
do
	echo $filename
	cat $filename | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | awk -F ',' -v OFS='\t' '{print $1,$2,$3}' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,int($9)","int($10)","int($11)}' > $filename'.int.bed'
	/gpfs/home/gzx103/scratch/vision_clustering/bedToBigBed $filename'.int.bed' /gpfs/home/gzx103/scratch/vision_clustering/mm10.chrom.sizes $filename'.int.bb'
done
