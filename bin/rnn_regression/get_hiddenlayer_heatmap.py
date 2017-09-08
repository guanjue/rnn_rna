import os,sys
from tensorflow.python.ops import rnn
import tensorflow as tf
import numpy as np
import math
import random
from sklearn import preprocessing
from sklearn.decomposition import PCA, KernelPCA
from sklearn import svm

def read_tab2matrix(inputfile,weight,colnum):
	import numpy as np
	data0 = open(inputfile,'r')
	#data0.readline()### skip the 1st line
	data01 = []
	seq_matrix = []
	for records in data0:
		tmp = [x.strip() for x in records.split('\t')]
		string = []
		for bp in tmp[colnum:]:
			if bp != 'nan':
				string.append(float(bp))
			else:
				string.append(0.0)
		for i in range(weight):
			seq_matrix.append(string)
	seq_matrix = np.array(seq_matrix)
	data0.close()
	return seq_matrix

def read_strtab2matrix(inputfile,weight,colnum):
	import numpy as np
	data0 = open(inputfile,'r')
	data0.readline()### skip the 1st line
	data01 = []
	seq_matrix = []
	for records in data0:
		tmp = [x.strip() for x in records.split('\t')]
		string = []
		for bp in tmp[colnum:]:
			if bp != 'nan':
				string.append(bp)
			else:
				string.append(0.0)
		for i in range(weight):
			seq_matrix.append(string)
	seq_matrix = np.array(seq_matrix)
	data0.close()
	return seq_matrix

def read_rna2matrix(inputfile,weight,colnum):
	import numpy as np
	data0 = open(inputfile,'r')
	data01 = []
	seq_matrix = []
	seq_matrix_sig = []
	for records in data0:
		tmp = [x.strip() for x in records.split('\t')]
		for i in range(weight):
			seq_matrix_sig.append([np.log2(float(tmp[colnum])+0.01)])
			if np.log2(float(tmp[colnum])+0.01) > 0:
				seq_matrix.append([1.0,0.0])
			else:
				seq_matrix.append([0.0,1.0])
	seq_matrix = np.array(seq_matrix)
	seq_matrix_sig = np.array(seq_matrix_sig)
	data0.close()
	return seq_matrix, seq_matrix_sig

def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

### randomly split testing and training data
def get_random_split_ids(data_matrix, random_seed, testing_data_p):
	### split testing & training
	n_sample=data_matrix.shape[0]
	index_array_p=np.arange(n_sample)
	np.random.seed(seed=random_seed)
	np.random.shuffle(index_array_p)
	test_id=index_array_p[0:int(n_sample*testing_data_p)]
	train_id=index_array_p[int(n_sample*testing_data_p):n_sample]
	return test_id, train_id


### read pk_info
data_0 = read_strtab2matrix('/Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/matrix_for_rnn/gene_list.lsk_bm.0.tab',1,0)
test_id, train_id = get_random_split_ids(data_0, 2017, 0.2)
pk_info = data_0[test_id,0:6]

### read rna-seq signal
data_rna_sig_all = read_tab2matrix('noatac3/data_rna_sig.lsk_bm.txt',1,0)
cell_types = ['cmp','mep','gmp','meg','cfuE_ad','ery_ad','gra_bm','g1e','er4']
cell_types_num = len(cell_types)+1
for ct in cell_types:
	data_rna_sig = read_tab2matrix('noatac3/data_rna_sig.'+ct+'.txt',1,0)
	data_rna_sig_all = np.concatenate((data_rna_sig_all, data_rna_sig), axis=1)
print('data_rna_sig_all.shape')
print(data_rna_sig_all.shape)

data_rna_sig_all_mean = np.mean(data_rna_sig_all,axis=1)
data_rna_sig_all_var = np.var(data_rna_sig_all,axis=1)
print('data_rna_sig_all_var.shape')
print(data_rna_sig_all_var.shape)
data_rna_sig_all_mean_sort_id = data_rna_sig_all_mean.argsort()

pk_info_mean_sort = pk_info[data_rna_sig_all_mean_sort_id]

data_rna_sig_all_mean_mean_sort = data_rna_sig_all[data_rna_sig_all_mean_sort_id]
data_rna_sig_all_var_mean_sort = data_rna_sig_all_var[data_rna_sig_all_mean_sort_id]

folder = 'withatac3'
cell_types = ['lsk_bm','cmp','mep','gmp','meg','cfuE_ad','ery_ad','gra_bm']#,'g1e','er4']
for ct in cell_types:
	label = read_tab2matrix(folder+'/y_conv_obs_all.'+ct+'.txt',1,0)
	pred = read_tab2matrix(folder+'/y_conv_pred_all.'+ct+'.txt',1,0)
	data_r1 = read_tab2matrix(folder+'/rnn_hidden_layer_pred_1.'+ct+'.txt',1,0)
	data_r2 = read_tab2matrix(folder+'/rnn_hidden_layer_pred_2.'+ct+'.txt',1,0)
	data_r3 = read_tab2matrix(folder+'/rnn_hidden_layer_pred_3.'+ct+'.txt',1,0)

	label_mean_sort = label[data_rna_sig_all_mean_sort_id]
	pred_mean_sort = pred[data_rna_sig_all_mean_sort_id]
	data_r1_mean_sort = data_r1[data_rna_sig_all_mean_sort_id] - pred_mean_sort[:,0].reshape((pred_mean_sort.shape[0]), 1)
	data_r2_mean_sort = data_r2[data_rna_sig_all_mean_sort_id] - pred_mean_sort[:,0].reshape((pred_mean_sort.shape[0]), 1)
	data_r3_mean_sort = data_r3[data_rna_sig_all_mean_sort_id] - pred_mean_sort[:,0].reshape((pred_mean_sort.shape[0]), 1)
	print('data_r1_mean_sort.shape')
	print(data_r1_mean_sort.shape)
	data_r1_mean_sort_out = np.concatenate((pk_info_mean_sort, data_r1_mean_sort), axis=1)
	data_r2_mean_sort_out = np.concatenate((pk_info_mean_sort, data_r2_mean_sort), axis=1)
	data_r3_mean_sort_out = np.concatenate((pk_info_mean_sort, data_r3_mean_sort), axis=1)

	#write2d_array(data_r1_mean_sort_out,'data_r1_mean_sort_out.'+ct+'.tab')
	#write2d_array(data_r2_mean_sort_out,'data_r2_mean_sort_out.'+ct+'.tab')
	#write2d_array(data_r3_mean_sort_out,'data_r3_mean_sort_out.'+ct+'.tab')

	data_r123_mean_sort_out = np.concatenate((pk_info_mean_sort, data_r1_mean_sort, data_r3_mean_sort, data_r2_mean_sort), axis=1)
	write2d_array(data_r123_mean_sort_out,folder+'/data_r123_mean_sort_out.'+ct+'.tab')

	data_rna_sig_all_mean_mean_sort0 = np.concatenate((pk_info_mean_sort, data_rna_sig_all_mean_mean_sort), axis=1)
	write2d_array(data_rna_sig_all_mean_mean_sort0,folder+'/data_rna_sig_all_mean_mean_sort.'+ct+'.tab')

	data_rna_sig_all_var_mean_sort0 = np.concatenate((pk_info_mean_sort, data_rna_sig_all_var_mean_sort.reshape((data_rna_sig_all_var_mean_sort.shape[0]), 1)), axis=1)
	write2d_array(data_rna_sig_all_var_mean_sort0,folder+'/data_rna_sig_all_var_mean_sort.'+ct+'.tab')

	label_mean_sort = np.concatenate((pk_info_mean_sort, pred_mean_sort, label_mean_sort), axis=1)
	write2d_array(label_mean_sort,folder+'/label_mean_sort.'+ct+'.tab')


