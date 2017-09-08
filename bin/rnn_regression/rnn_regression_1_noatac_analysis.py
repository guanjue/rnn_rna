import os,sys
from tensorflow.python.ops import rnn
import tensorflow as tf
import numpy as np
import math
import random
from sklearn import preprocessing
from sklearn.decomposition import PCA, KernelPCA
from sklearn import svm
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from scipy import stats

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
	for records in data0:
		tmp = [x.strip() for x in records.split('\t')]
		for i in range(weight):
			seq_matrix.append([np.log2(float(tmp[colnum])+0.01)])
	seq_matrix = np.array(seq_matrix)
	data0.close()
	return seq_matrix

def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

### read atac_seq pk
data_atac = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn_atacpk/atac_21k_binary_split.tab',1,6)
### read ideas labels
data_ideas_matrix = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn/matrix_for_rnn/gene_list.er4.0.tab',1,6)
#data_ideas_matrix = data_ideas_matrix * data_atac
data_ideas_matrix = data_ideas_matrix.reshape([data_ideas_matrix.shape[0],data_ideas_matrix.shape[1],1,1])
for i in range(1,17):
	print(i)
	data_ideas_tmp = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn/matrix_for_rnn/gene_list.er4.'+str(i)+'.tab',1,6)
	#data_ideas_tmp = data_ideas_tmp * data_atac
	data_ideas_tmp = data_ideas_tmp.reshape([data_ideas_tmp.shape[0],data_ideas_tmp.shape[1],1,1])
	data_ideas_matrix = np.concatenate((data_ideas_matrix, data_ideas_tmp), axis=3)

print('data_ideas_matrix.shape')
print(data_ideas_matrix.shape)
### read rna-seq
data_rna = read_rna2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/rsem/er4.rsem.txt',1,4)

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

pos_test_id, pos_train_id = get_random_split_ids(data_rna, 2017, 0.2)

x_test = data_ideas_matrix[pos_test_id,:]
x_train = data_ideas_matrix[pos_train_id,:]
x_all = data_ideas_matrix[:,:]

y_test = data_rna[pos_test_id,:]
y_train = data_rna[pos_train_id,:]
y_all = data_rna[:,:]

'''
cell_types = ['er4']
cell_types_num = len(cell_types)+1
for ct in cell_types:
	### read atac_seq pk
	data_atac = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn_atacpk/atac_21k_binary_split.tab',1,6)
	### read ideas labels
	data_ideas_matrix = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn/matrix_for_rnn/gene_list.'+ct+'.0.tab',1,6)
	#data_ideas_matrix = data_ideas_matrix * data_atac
	data_ideas_matrix = data_ideas_matrix.reshape([data_ideas_matrix.shape[0],data_ideas_matrix.shape[1],1,1])
	for i in range(1,17):
		print(i)
		data_ideas_tmp = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn/matrix_for_rnn/gene_list.'+ct+'.'+str(i)+'.tab',1,6)
		#data_ideas_tmp = data_ideas_tmp * data_atac
		data_ideas_tmp = data_ideas_tmp.reshape([data_ideas_tmp.shape[0],data_ideas_tmp.shape[1],1,1])
		data_ideas_matrix = np.concatenate((data_ideas_matrix, data_ideas_tmp), axis=3)

	print('data_ideas_matrix.shape')
	print(data_ideas_matrix.shape)
	### read rna-seq
	data_rna = read_rna2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/rsem/'+ct+'.rsem.txt',1,4)

	x_test = np.concatenate((x_test, data_ideas_matrix[pos_test_id,:]), axis = 0)
	x_train = np.concatenate((x_train, data_ideas_matrix[pos_train_id,:]), axis = 0)

	y_test = np.concatenate((y_test, data_rna[pos_test_id,:]), axis = 0)
	y_train = np.concatenate((y_train, data_rna[pos_train_id,:]), axis = 0)

	print('x_train.shape')
	print(x_train.shape)
	print('y_train.shape')
	print(y_train.shape)
'''


sec_d=50 ### upstream_expand_num+gene_len_num
sec_d_rev=50+50-1 ### downstream_expand_num-1 
genebody_len = 50
sec_d_genebody = sec_d+genebody_len

all_len=50+50+50
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
x = tf.placeholder(tf.float32, shape=[None, all_len, thr_d, for_d])
#x_bn=tf.contrib.layers.batch_norm(x,center=True, scale=True,scope='bn')

y_ = tf.placeholder(tf.float32, shape=[None, 1])
keep_prob1 = tf.placeholder(tf.float32)


x_transpose = tf.transpose(x, [1, 0, 2, 3])
x_reshape = tf.reshape(x_transpose, [-1, for_d])
#print(x.shape)
#print(x_reshape.shape)
h_conv1_out = tf.split(x_reshape,axis=0, num_or_size_splits=all_len)
#tf.split(x_reshape, num_or_size_splits=sec_d, axis=1)
#print(h_conv1_out)
rnn_cell_num = 32
#rnn_cell = tf.contrib.rnn.LSTMCell(rnn_cell_num, forget_bias=0,use_peepholes=False)
rnn_cell = tf.contrib.rnn.GRUCell(rnn_cell_num)

outputs, states = tf.contrib.rnn.static_rnn(rnn_cell, h_conv1_out, dtype=tf.float32)

W_atten1=weight_variable([rnn_cell_num, 1])
b_atten1 = bias_variable([rnn_cell_num])




lstm_tail=(outputs[-1])


y_conv = (tf.matmul(lstm_tail, W_atten1)) 

cost=tf.reduce_mean(tf.pow(y_conv-y_, 2))
train_step = tf.train.AdamOptimizer(training_speed).minimize(cost)


sess = tf.InteractiveSession()
sess.run(tf.initialize_all_variables())
sess.run(tf.initialize_local_variables())




print('Start!!! RNN')
saver = tf.train.Saver()
saver.restore(sess, "trained_rnn_1_noatac.ckpt")


y_p = sess.run(y_conv, feed_dict={x: x_train[:,:,:,:], keep_prob1: 1})
y_obs=y_train[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_test=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, TRAINING accuracy same cell!!! R2: %g"%(i, r2_test)  )
print(stats.spearmanr(y_obs,y_p))
print('stats.pearsonr(y_obs,y_p)')
print(stats.pearsonr(y_obs,y_p))


y_p = sess.run(y_conv, feed_dict={x: x_test[:,:,:,:], keep_prob1: 1})
y_obs=y_test[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_test=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, TESTING accuracy same cell!!! R2: %g"%(i, r2_test)  )
print(stats.spearmanr(y_obs,y_p))
print(stats.spearmanr(y_obs,y_p))
print('stats.pearsonr(y_obs,y_p)')
print(stats.pearsonr(y_obs,y_p))


y_p = sess.run(y_conv, feed_dict={x: x_all[:,:,:,:], keep_prob1: 1})
y_obs=y_all[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_test=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, ALL accuracy same cell!!! R2: %g"%(i, r2_test)  )
print(stats.spearmanr(y_obs,y_p))
print('stats.pearsonr(y_obs,y_p)')
print(stats.pearsonr(y_obs,y_p))


rnn_hidden_layer=np.array(sess.run(outputs, feed_dict={x: x_all[:,:,:,:], keep_prob1: 1}))
print('rnn_hidden_layer.shape')
print(rnn_hidden_layer.shape)

rnn_hidden_layer_t=np.transpose(np.array(rnn_hidden_layer,dtype=float),(1,0,2))
print('rnn_hidden_layer_t.shape')
print(rnn_hidden_layer_t.shape)

fc_layer=np.array(sess.run(W_atten1, feed_dict={x: x_all[:,:,:,:], keep_prob1: 1}))
print('fc_layer.shape')
print(fc_layer.shape)

rnn_hidden_layer_pred_dif = (np.dot(rnn_hidden_layer_t, fc_layer)[:,:,0] - y_all)
rnn_hidden_layer_pred = (np.dot(rnn_hidden_layer_t, fc_layer)[:,:,0])
#rnn_hidden_layer_pred_dif[rnn_hidden_layer_pred_dif>=6] = 6
#rnn_hidden_layer_pred_dif[rnn_hidden_layer_pred_dif<=-6] = -6
print('rnn_hidden_layer_pred_dif.shape')
print(rnn_hidden_layer_pred_dif.shape)




print('kmeans')
kmeans_k=3
kmeans = KMeans(n_clusters=kmeans_k, init='k-means++', max_iter=1000, n_init=1, verbose=0, random_state=2017).fit(rnn_hidden_layer_pred_dif)
kmeans_label=kmeans.labels_

rnn_hidden_layer_pred_dif_k5 = rnn_hidden_layer_pred_dif[kmeans_label.argsort()]
y_all_k5 = y_all[kmeans_label.argsort()]

write2d_array(rnn_hidden_layer_pred_dif_k5,'rnn_hidden_layer_pred_dif_k5.txt')
fig = plt.figure(figsize=(200,100))
plt.imshow(rnn_hidden_layer_pred_dif_k5, cmap='bwr', interpolation='nearest', aspect='auto')
plt.savefig('rnn_hidden_layer_pred_dif_k5.pdf', dpi=300)


write2d_array(y_all_k5,'rnn_hidden_layer_pred_dif_k5.txt')
fig = plt.figure(figsize=(200,100))
plt.imshow(y_all_k5, cmap='bwr', interpolation='nearest', aspect='auto')
plt.savefig('y_all_k5.pdf', dpi=300)


'''
step 49000, TESTING accuracy same cell!!! R2: 0.81622
-2.35682180703
step 49100, training accuracy same cell!!! R2: 0.833024
step 49200, training accuracy same cell!!! R2: 0.825792
step 49300, training accuracy same cell!!! R2: 0.881536
step 49400, training accuracy same cell!!! R2: 0.822821
step 49500, training accuracy same cell!!! R2: 0.816171
step 49600, training accuracy same cell!!! R2: 0.815814
step 49700, training accuracy same cell!!! R2: 0.884659
step 49800, training accuracy same cell!!! R2: 0.798853
step 49900, training accuracy same cell!!! R2: 0.809067
-2.8011434812
step 49999, TESTING accuracy same cell!!! R2: 0.818298

real	149m40.261s
user	159m30.246s
sys	31m59.367s
'''