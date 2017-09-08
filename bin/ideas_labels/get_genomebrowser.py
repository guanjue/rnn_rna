def get_genomebrowser(input_all,input_celltype_ideas,output):
	data_all0=open(input_all,'r')
	data_all1=[]
	for records in data_all0:
		tmp=[x.strip() for x in records.split('\t')]
		data_all1.append(tmp)

	data_ideas0=open(input_celltype_ideas,'r')
	data_ideas1={}
	data_ideas1_maxcover={}
	data_ideas1_middist={}
	data_ideas1_statelen={}
	k=0
	for records in data_ideas0:
		tmp=[x.strip() for x in records.split('\t')]
		if not (tmp[3] in data_ideas1):
			data_ideas1[tmp[3]] = tmp
			data_ideas1_maxcover[tmp[3]] = int(tmp[4])
			data_ideas1_middist[tmp[3]] = float(tmp[7])
			data_ideas1_statelen[tmp[3]] = int(tmp[8])
		elif int(tmp[4]) > data_ideas1_maxcover[tmp[3]]:
			### if multiple cover; select the highest covering state
			data_ideas1[tmp[3]] = tmp
			data_ideas1_maxcover[tmp[3]] = int(tmp[4])
			data_ideas1_middist[tmp[3]] = float(tmp[7])
			data_ideas1_statelen[tmp[3]] = int(tmp[8])
		elif tmp[4] == data_ideas1_maxcover[tmp[3]]:
			if int(tmp[7]) < data_ideas1_middist[tmp[3]]:
				### if cover the same; check mid point distance
				data_ideas1[tmp[3]] = tmp
				data_ideas1_maxcover[tmp[3]] = int(tmp[4])
				data_ideas1_middist[tmp[3]] = float(tmp[7])
				data_ideas1_statelen[tmp[3]] = int(tmp[8])		
			elif tmp[7] == data_ideas1_middist[tmp[3]]:
				if int(tmp[8]) < data_ideas1_statelen[tmp[3]]:
					### if cover same & mid point distance same; check state len 
					data_ideas1[tmp[3]] = tmp
					data_ideas1_maxcover[tmp[3]] = int(tmp[4])
					data_ideas1_middist[tmp[3]] = float(tmp[7])
					data_ideas1_statelen[tmp[3]] = int(tmp[8])				
				else:
					k=k+1
					print('problem!')
					print(k)

	result=open(output,'w')

	for records in data_all1:
		if records[0] in data_ideas1:
			tmp=data_ideas1[records[0]]
			result.write(tmp[0]+'\t'+tmp[1]+'\t'+tmp[2]+'\t'+tmp[5]+'\t'+'1000'+'\t'+'.'+'\t'+tmp[1]+'\t'+tmp[2]+'\t'+tmp[6]+'\n')
		else:
			tmp=records
			result.write(tmp[2]+'\t'+tmp[3]+'\t'+tmp[4]+'\t'+'NA'+'\t'+'1000'+'\t'+'.'+'\t'+tmp[3]+'\t'+tmp[4]+'\t'+'255,255,255'+'\n')
	result.close()
	data_all0.close()
	data_ideas0.close()


############################################################################
### python get_genomebrowser.py -a data_matrix_binarypattern_sorted.txt -i used_DNA_intervals.1.bed -o used_DNA_intervals.1.colored.bed
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h:a:i:o:")
	except getopt.GetoptError:
		print 'python get_genomebrowser.py -a data_matrix_binarypattern_sorted.txt -i used_DNA_intervals.1.bed -o used_DNA_intervals.1.colored.bed'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python get_genomebrowser.py -a data_matrix_binarypattern_sorted.txt -i used_DNA_intervals.1.bed -o used_DNA_intervals.1.colored.bed'
			sys.exit()
		elif opt=="-a":
			input_all=str(arg.strip())
		elif opt=="-i":
			input_celltype_ideas=str(arg.strip())
		elif opt=="-o":
			output=str(arg.strip())

	get_genomebrowser(input_all,input_celltype_ideas,output)
if __name__=="__main__":
	main(sys.argv[1:])






