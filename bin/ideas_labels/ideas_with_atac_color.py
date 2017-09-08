import colorsys
import numpy as np

def ideas_with_atac_color(inputbed, intput_atac1, intput_atac2, outputfile, std_num):
	data1a=open(intput_atac1,'r')
	data1b=open(intput_atac2,'r')

	### get atac-seq average reads count relative to max
	average_atac_sig=[]
	for atac_sigsa,atac_sigsb in zip(data1a,data1b):
		tmp_siga=[x.strip() for x in atac_sigsa.split('\t')]
		tmp_sigb=[x.strip() for x in atac_sigsb.split('\t')]
		average_atac_sig.append( (float(tmp_siga[9])+float(tmp_sigb[9]))/2 )
	data1a.close()
	data1b.close()
	### get atac-seq reads count relative to mean+1.5*std (in case of outliers)
	average_atac_sig=np.array(average_atac_sig,dtypes=float)
	average_atac_sig_mean=np.mean(average_atac_sig)
	average_atac_sig_std=np.std(average_atac_sig)
	average_atac_sig_threshold=average_atac_sig_mean+std_num*average_atac_sig_std
	average_atac_sig_percentage_withthreshold=np.clip(a, 0.0, average_atac_sig_threshold)/average_atac_sig_threshold

	### change color
	data0=open(inputbed,'r')
	result=open(outputfile,'w')
	data01=[]
	i=0
	for records in data0:
		tmp_pk=[x.strip() for x in records.split('\t')]

		### rgb color to hsv
		tmp_od_rgb_color=tmp_pk[8].split(',')
		tmp_od_hsv_color=colorsys.rgb_to_hsv(int(tmp_od_rgb_color[0]),int(tmp_od_rgb_color[1]),int(tmp_od_rgb_color[2]))

		### get Saturation atac_seq reads count based
		tmp_sat_percentage=average_atac_sig_percentage_withthreshold[i]
		tmp_atacsat_hsv_color_array=[tmp_od_hsv_color[0],tmp_od_hsv_color[0]*tmp_sat_percentage,tmp_od_hsv_color[0]]
		tmp_atacsat_rgb_color_array=colorsys.hsv_to_rgb(tmp_atacsat_hsv_color_array[0], tmp_atacsat_hsv_color_array[1], tmp_atacsat_hsv_color_array[2])

		### write data
		for records in range(0,len(tmp_pk)-1):
			result.write(records+'\t')

		result.write(tmp_atacsat_rgb_color_array[0]+','+tmp_atacsat_rgb_color_array[1]+','+tmp_atacsat_rgb_color_array[2]+'\n')

	data0.close()
	result.close()

############################################################################
### python ideas_with_atac_color.py -i tmp.bed -a tmp_bam_a.bed -b tmp_bam_b.bed -o tmp.atac_colored.bed -s 1.5
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h:i:a:b:o:s:")
	except getopt.GetoptError:
		print 'python ideas_with_atac_color.py -i tmp.bed -a tmp_bam1.bed -b tmp_bam2.bed -o tmp.atac_colored.bed -s 1.5'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python ideas_with_atac_color.py -i tmp.bed -a tmp_bam1.bed -b tmp_bam2.bed -o tmp.atac_colored.bed -s 1.5'
			sys.exit()
		elif opt=="-i":
			inputbed=str(arg.strip())
		elif opt=="-a":
			intput_atac1=str(arg.strip())
		elif opt=="-b":
			intput_atac2=str(arg.strip())
		elif opt=="-o":
			outputfile=str(arg.strip())
		elif opt=="-s":
			std_num=float(arg.strip())

	ideas_with_atac_color(inputbed, intput_atac1, intput_atac2, outputfile, std_num)
if __name__=="__main__":
	main(sys.argv[1:])

