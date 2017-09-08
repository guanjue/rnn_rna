def bed2matrix(inputfile,outputfile):
	data0=open(inputfile,'r')
	data1=[]
	for records in data0:
		tmp=[x.strip() for x in records.split('\t')]
		info=[0]*17
		if tmp[3]=='NA':
			data1.append([tmp[0],tmp[1],tmp[2],tmp[9]]+info)
		else:
			state_num=int(tmp[3])
			info[state_num]=1
			data1.append([tmp[0],tmp[1],tmp[2],tmp[9]]+info)
	data0.close()

	result=open(outputfile,'w')
	for records in data1:
		for i in range(0,len(records)-1):
			result.write(str(records[i])+'\t')
		result.write(str(records[len(records)-1])+'\n')
	result.close()



############################################################################
### python bed2matrix.py -i used_DNA_intervals.1.forbed2matrix.bed -o homerTable3.peaks.filtered.bed2matrix.1.txt
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h:i:o:")
	except getopt.GetoptError:
		print 'python prepare_bed2matrix.py -i used_DNA_intervals.1.forbed2matrix.bed -o homerTable3.peaks.filtered.bed2matrix.1.txt'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python prepare_bed2matrix.py -i used_DNA_intervals.1.forbed2matrix.bed -o homerTable3.peaks.filtered.bed2matrix.1.txt'
			sys.exit()
		elif opt=="-i":
			inputfile=str(arg.strip())
		elif opt=="-o":
			outputfile=str(arg.strip())

	bed2matrix(inputfile,outputfile)
if __name__=="__main__":
	main(sys.argv[1:])

