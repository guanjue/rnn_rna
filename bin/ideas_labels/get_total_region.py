def get_total_region(input_file,output_file,info):
	data=open(input_file,'r')
	data1=0
	for records in data:
		tmp=[x.strip() for x in records.split('\t')]
		if int(tmp[4]) == 1:
			data1=data1+int(tmp[2])-int(tmp[1])
	result=open(output_file,'w')
	result.write(str(data1)+'\t'+info+'\n')

############################################################################
### python get_total_region.py -i <input file name> -o <output file name> -f <info>
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h:i:o:f:")
	except getopt.GetoptError:
		print 'python get_total_region.py -i <input file name> -o <output file name> -f <info>'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python get_total_region.py -i <input file name> -o <output file name> -f <info>'
			sys.exit()
		elif opt=="-i":
			input_file=str(arg.strip())
		elif opt=="-o":
			output_file=str(arg.strip())
		elif opt=="-f":
			info=str(arg.strip())


	get_total_region(input_file,output_file,info)
if __name__=="__main__":
	main(sys.argv[1:])
