def split_bed_window_size(inpufile,outputfile,win):
	data=open(inpufile,'r')
	result=open(outputfile,'w')

	### split bed based on window size
	for records in data:
		tmp=[x.strip() for x in records.split('\t')]
		num=( int(tmp[2])-int(tmp[1]) )/win
		for i in range(0,num):
			result.write(tmp[0]+'\t'+str(int(tmp[1])+win*i)+'\t'+str(int(tmp[1])+win*i+win)+'\n' )

	data.close()
	result.close()

############################################################################
### python split_bed_window_size.py -i ideas_label_bed/ideasVision1.sort.bed -o ideas_bed_split.bed -w 200
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h:i:o:w:")
	except getopt.GetoptError:
		print 'python split_bed_window_size.py -i input_bedfile.bed -o output_bedfile.bed -w window_size'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python split_bed_window_size.py -i input_bedfile.bed -o output_bedfile.bed -w window_size'
			sys.exit()
		elif opt=="-i":
			inpufile=str(arg.strip())
		elif opt=="-o":
			outputfile=str(arg.strip())
		elif opt=="-w":
			win=int(arg.strip())

	split_bed_window_size(inpufile, outputfile, win)
if __name__=="__main__":
	main(sys.argv[1:])