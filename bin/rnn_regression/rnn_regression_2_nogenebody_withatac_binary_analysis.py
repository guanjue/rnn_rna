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
	data0.readline()### skip the 1st line
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


sec_d=47 ### upstream_expand_num+gene_len_num
sec_d_rev=47+56-1 ### downstream_expand_num-1 
genebody_len = 56
sec_d_genebody = sec_d+genebody_len

all_len=47+56+47
thr_d=1
for_d=17
filter1_size_1st=8
filter1_size_2nd=1
filter1_size_out=32
filter1_max_pool_size=4

filter2_size_1st=8
filter2_size_out=32
filter2_max_pool_size=4

filter1_size2=1
first_filter_out=32
full_cn1_out=16
full_cn2_out=512
full_cn3_out=512
full_cn4_out=512
full_cn5_out=128
full_cn6_out=64
full_cn7_out=64
full_cn8_out=32
full_cn9_out=1

tail_size=28

epsilon = 0.0005
iter_num=50000
batch_size=32
training_speed=0.0005

rnn_cell_num = 32

def weight_variable(shape):
	initial = tf.truncated_normal(shape, stddev=0.01)
	return tf.Variable(initial)

def bias_variable(shape):
	initial = tf.constant(0.01, shape=shape)
	return tf.Variable(initial)

### Convolution and Pooling
def conv2d(x, W):
	return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')

def conv2d_1(x, W):
	return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')

def max_pool_n(x, max_pool_size):
	return tf.nn.max_pool(x, ksize=[1, max_pool_size, max_pool_size, 1], strides=[1, max_pool_size, max_pool_size, 1], padding='SAME')

######### Tensorflow model
tf.set_random_seed(2017)
x_1 = tf.placeholder(tf.float32, shape=[None, sec_d, thr_d, for_d])

y_ = tf.placeholder(tf.float32, shape=[None, 2])
keep_prob1 = tf.placeholder(tf.float32)


x_transpose_1 = tf.transpose(x_1, [1, 0, 2, 3])
x_reshape_1 = tf.reshape(x_transpose_1, [-1, for_d])

h_conv1_out_1 = tf.split(x_reshape_1,axis=0, num_or_size_splits=sec_d)
rnn_cell_1 = tf.contrib.rnn.LSTMCell(rnn_cell_num, forget_bias=0.0)

outputs_1, states_1 = tf.contrib.rnn.static_rnn(rnn_cell_1, h_conv1_out_1, dtype=tf.float32)

W_atten1_1=weight_variable([rnn_cell_num, 2])
b_atten1_1 = bias_variable([rnn_cell_num])



x_2 = tf.placeholder(tf.float32, shape=[None, sec_d, thr_d, for_d])

x_transpose_2 = tf.transpose(x_2, [1, 0, 2, 3])
x_reshape_2 = tf.reshape(x_transpose_2, [-1, for_d])

h_conv1_out_2 = tf.split(x_reshape_2,axis=0, num_or_size_splits=sec_d)
rnn_cell2 = tf.contrib.rnn.LSTMCell(rnn_cell_num,reuse=True, forget_bias=0.0)

outputs_2, states_2 = tf.contrib.rnn.static_rnn(rnn_cell2, h_conv1_out_2, dtype=tf.float32)

W_atten1_2=weight_variable([rnn_cell_num, 2])
b_atten1_2 = bias_variable([rnn_cell_num])

x_3 = tf.placeholder(tf.float32, shape=[None, genebody_len, thr_d, for_d])

x_transpose_3 = tf.transpose(x_3, [1, 0, 2, 3])
x_reshape_3 = tf.reshape(x_transpose_3, [-1, for_d])

h_conv1_out_3 = tf.split(x_reshape_3,axis=0, num_or_size_splits=genebody_len)
rnn_cell3 = tf.contrib.rnn.LSTMCell(rnn_cell_num,reuse=True, forget_bias=0.0)

outputs_3, states_3 = tf.contrib.rnn.static_rnn(rnn_cell3, h_conv1_out_3, dtype=tf.float32 )

W_atten1_3=weight_variable([rnn_cell_num, 2])
b_atten1_3 = bias_variable([rnn_cell_num])



lstm_tail_1=(outputs_1[-1])
lstm_tail_2=(outputs_2[-1])
lstm_tail_3=(outputs_3[-1])


y_conv = (tf.matmul(lstm_tail_1, W_atten1_1)) + (tf.matmul(lstm_tail_2, W_atten1_2)) #+ (tf.matmul(lstm_tail_3, W_atten1_3))

y_conv_pred = tf.nn.softmax(y_conv)
### Evaluate the Model
#cross_entropy = -tf.reduce_sum(y_*tf.log(y_conv))
cross_entropy = tf.nn.softmax_cross_entropy_with_logits(logits=y_conv, labels=y_) 
train_step = tf.train.AdamOptimizer(training_speed).minimize(cross_entropy)
correct_prediction = tf.equal(tf.argmax(y_conv,1), tf.argmax(y_,1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

au_roc = tf.metrics.auc(y_, y_conv_pred, curve='ROC')
au_prc = tf.metrics.auc(y_, y_conv_pred, curve='PR')


sess = tf.InteractiveSession()
sess.run(tf.initialize_all_variables())
sess.run(tf.initialize_local_variables())




print('Start!!! RNN')
saver = tf.train.Saver()
saver.restore(sess, "rnn_regression_2_nogenebody_withatac_binary.ckpt")





### read atac_seq pk

cell_types = ['lsk_bm','cmp','mep','gmp','meg','cfuE_ad','ery_ad','gra_bm']#,'g1e','er4']
cell_types_num = len(cell_types)+1
for ct in cell_types:
	data_atac = read_tab2matrix('/Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/matrix_for_rnn_atacpk/atac_21k_binary_split.tab',1,6)
	data_atac_cellspe = read_tab2matrix('/Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/matrix_for_rnn_atacpk/atac_21k_binary_split.'+ct+'.tab',1,6)
	### read ideas labels
	data_ideas_matrix = read_tab2matrix('/Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/matrix_for_rnn/gene_list.'+ct+'.0.tab',1,6)
	data_ideas_matrix = data_atac * data_ideas_matrix
	data_ideas_matrix = data_ideas_matrix.reshape([data_ideas_matrix.shape[0],data_ideas_matrix.shape[1],1,1])
	for i in range(1,17):
		print(i)
		data_ideas_tmp = read_tab2matrix('/Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/matrix_for_rnn/gene_list.'+ct+'.'+str(i)+'.tab',1,6)
		#data_ideas_tmp = data_ideas_tmp * data_atac
		data_ideas_tmp = data_ideas_tmp.reshape([data_ideas_tmp.shape[0],data_ideas_tmp.shape[1],1,1])
		data_ideas_matrix = np.concatenate((data_ideas_matrix, data_ideas_tmp), axis=3)

	print('data_ideas_matrix.shape')
	print(data_ideas_matrix.shape)
	### read rna-seq
	data_rna, data_rna_sig = read_rna2matrix('/Users/gzx103/Documents/zhang_lab/projects/atacseq_vision/rna_seq_pred/rsem/'+ct+'.rsem.txt',1,4)

	pos_test_id, pos_train_id = get_random_split_ids(data_rna, 2017, 0.2)

	x_test = data_ideas_matrix[pos_test_id,:]
	x_train = data_ideas_matrix[pos_train_id,:]

	y_test = data_rna[pos_test_id,:]
	y_train = data_rna[pos_train_id,:]

	y_test_sig = data_rna_sig[pos_test_id,:]
	y_train_sig = data_rna_sig[pos_train_id,:]

	test_accuracy = accuracy.eval(feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0})
	print("step %d, testing accuracy same cell!!! R2: %g"%(i, test_accuracy) )

	test_roc = sess.run(au_roc, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0})
	test_prc = sess.run(au_prc, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0})

	test_roc = sess.run(au_roc, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0})
	test_prc = sess.run(au_prc, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0})

	print("step %d, testing accuracy same cell!!! R2: %g"%(i, test_accuracy) )
	print("step %d, testing au roc same cell!!! R2: %g"%(i, test_roc[0]) )
	print("step %d, testing au prc same cell!!! R2: %g"%(i, test_prc[0]) )

	y_conv_pred_test = sess.run(y_conv_pred, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0})
	y_conv_obs_test = sess.run(y_, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0})
	write2d_array(y_conv_pred_test,'withatac2_nogenebody/y_conv_pred_all.'+ct+'.txt')
	write2d_array(y_conv_obs_test,'withatac2_nogenebody/y_conv_obs_all.'+ct+'.txt')
	write2d_array(y_test_sig,'withatac2_nogenebody/data_rna_sig.'+ct+'.txt')




	rnn_hidden_layer_1=np.array(sess.run(outputs_1, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0}))
	print('rnn_hidden_layer_1.shape')
	print(rnn_hidden_layer_1.shape)

	rnn_hidden_layer_t_1=np.transpose(np.array(rnn_hidden_layer_1,dtype=float),(1,0,2))
	print('rnn_hidden_layer_t_1.shape')
	print(rnn_hidden_layer_t_1.shape)

	fc_layer_1=np.array(sess.run(W_atten1_1, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0}))
	print('fc_layer_1.shape')
	print(fc_layer_1.shape)

	rnn_hidden_layer_2=np.array(sess.run(outputs_2, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0}))
	print('rnn_hidden_layer_2.shape')
	print(rnn_hidden_layer_2.shape)

	rnn_hidden_layer_t_2=np.transpose(np.array(rnn_hidden_layer_2,dtype=float),(1,0,2))
	print('rnn_hidden_layer_t_2.shape')
	print(rnn_hidden_layer_t_2.shape)

	fc_layer_2=np.array(sess.run(W_atten1_2, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0}))
	print('fc_layer_2.shape')
	print(fc_layer_2.shape)

	rnn_hidden_layer_3=np.array(sess.run(outputs_3, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0}))
	print('rnn_hidden_layer_3.shape')
	print(rnn_hidden_layer_3.shape)

	rnn_hidden_layer_t_3=np.transpose(np.array(rnn_hidden_layer_3,dtype=float),(1,0,2))
	print('rnn_hidden_layer_t_3.shape')
	print(rnn_hidden_layer_t_3.shape)

	fc_layer_3=np.array(sess.run(W_atten1_3, feed_dict={x_1: x_test[:,0:sec_d,:,:], x_2: x_test[:,:sec_d_rev:-1,:,:], x_3: x_test[:,sec_d:sec_d_genebody,:,:], y_: y_test, keep_prob1: 1.0}))
	print('fc_layer_3.shape')
	print(fc_layer_3.shape)

	#rnn_hidden_layer_pred_dif_1 = (np.dot(rnn_hidden_layer_t_1, fc_layer_1)[:,:,0] )#- y_test)
	print(np.dot(rnn_hidden_layer_t_2, fc_layer_2)[:,:,0].shape)
	print(np.dot(rnn_hidden_layer_t_2, fc_layer_2)[:,:,0][:,-1].shape)
	print(np.repeat(np.dot(rnn_hidden_layer_t_2, fc_layer_2)[:,:,0][:,-1],50,axis=0).shape)
	print(np.repeat(np.reshape(np.dot(rnn_hidden_layer_t_2, fc_layer_2)[:,:,0][:,-1],(y_test.shape[0],1)),50,axis=1).shape)

	rnn_hidden_layer_pred_1 = (np.dot(rnn_hidden_layer_t_1, fc_layer_1)[:,:,0]) + (np.reshape(np.dot(rnn_hidden_layer_t_2, fc_layer_2)[:,:,0][:,-1],(y_test.shape[0],1))) + (np.reshape(np.dot(rnn_hidden_layer_t_3, fc_layer_3)[:,:,0][:,-1],(y_test.shape[0],1)))
	rnn_hidden_layer_pred_1_neg = (np.dot(rnn_hidden_layer_t_1, fc_layer_1)[:,:,1]) + (np.reshape(np.dot(rnn_hidden_layer_t_2, fc_layer_2)[:,:,1][:,-1],(y_test.shape[0],1))) + (np.reshape(np.dot(rnn_hidden_layer_t_3, fc_layer_3)[:,:,1][:,-1],(y_test.shape[0],1)))
	rnn_hidden_layer_pred_softmax_1 = np.exp(rnn_hidden_layer_pred_1) / (np.exp(rnn_hidden_layer_pred_1) + np.exp(rnn_hidden_layer_pred_1_neg))

	rnn_hidden_layer_pred_2 = (np.dot(rnn_hidden_layer_t_2, fc_layer_2)[:,:,0]) + (np.reshape(np.dot(rnn_hidden_layer_t_1, fc_layer_1)[:,:,0][:,-1],(y_test.shape[0],1))) + (np.reshape(np.dot(rnn_hidden_layer_t_3, fc_layer_3)[:,:,0][:,-1],(y_test.shape[0],1)))
	rnn_hidden_layer_pred_2_neg = (np.dot(rnn_hidden_layer_t_2, fc_layer_2)[:,:,1]) + (np.reshape(np.dot(rnn_hidden_layer_t_1, fc_layer_1)[:,:,1][:,-1],(y_test.shape[0],1))) + (np.reshape(np.dot(rnn_hidden_layer_t_3, fc_layer_3)[:,:,1][:,-1],(y_test.shape[0],1)))	
	rnn_hidden_layer_pred_softmax_2 = np.exp(rnn_hidden_layer_pred_2) / (np.exp(rnn_hidden_layer_pred_2) + np.exp(rnn_hidden_layer_pred_2_neg))

	rnn_hidden_layer_pred_3 = (np.dot(rnn_hidden_layer_t_3, fc_layer_3)[:,:,0]) + (np.reshape(np.dot(rnn_hidden_layer_t_2, fc_layer_2)[:,:,0][:,-1],(y_test.shape[0],1))) + (np.reshape(np.dot(rnn_hidden_layer_t_1, fc_layer_1)[:,:,0][:,-1],(y_test.shape[0],1)))
	rnn_hidden_layer_pred_3_neg = (np.dot(rnn_hidden_layer_t_3, fc_layer_3)[:,:,1]) + (np.reshape(np.dot(rnn_hidden_layer_t_2, fc_layer_2)[:,:,1][:,-1],(y_test.shape[0],1))) + (np.reshape(np.dot(rnn_hidden_layer_t_1, fc_layer_1)[:,:,1][:,-1],(y_test.shape[0],1)))
	rnn_hidden_layer_pred_softmax_3 = np.exp(rnn_hidden_layer_pred_3) / (np.exp(rnn_hidden_layer_pred_3) + np.exp(rnn_hidden_layer_pred_3_neg))

	#rnn_hidden_layer_pred_dif[rnn_hidden_layer_pred_dif>=6] = 6
	#rnn_hidden_layer_pred_dif[rnn_hidden_layer_pred_dif<=-6] = -6

	write2d_array(rnn_hidden_layer_pred_softmax_1,'withatac2_nogenebody/rnn_hidden_layer_pred_1.'+ct+'.txt')
	write2d_array(rnn_hidden_layer_pred_softmax_2,'withatac2_nogenebody/rnn_hidden_layer_pred_2.'+ct+'.txt')
	write2d_array(rnn_hidden_layer_pred_softmax_3,'withatac2_nogenebody/rnn_hidden_layer_pred_3.'+ct+'.txt')

	write2d_array(y_test_sig,'withatac2_nogenebody/y_test_sig.'+ct+'.txt')
	write2d_array(y_test,'withatac2_nogenebody/y_test.'+ct+'.txt')

'''
step 16, testing accuracy same cell!!! R2: 0.903461
step 16, testing accuracy same cell!!! R2: 0.903461
step 16, testing au roc same cell!!! R2: 0.963226
step 16, testing au prc same cell!!! R2: 0.96304
'''
