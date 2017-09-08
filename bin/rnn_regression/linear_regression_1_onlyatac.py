import os,sys
from tensorflow.python.ops import rnn
import tensorflow as tf
import numpy as np
import math
import random
from sklearn import preprocessing
from sklearn.decomposition import PCA, KernelPCA
from sklearn import svm
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import learning_curve
from sklearn.kernel_ridge import KernelRidge

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
### read atac_seq pk
data_atac = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn_atacpk/atac_21k_binary_split.tab',1,6)
data_atac_cellspe = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn_atacpk/atac_21k_binary_split.gmp.tab',1,6)

### read ideas labels
#data_ideas_matrix = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn/matrix_for_rnn/gene_list.gmp.0.tab',1,6)
data_ideas_matrix = data_atac * data_atac_cellspe
#data_ideas_matrix = data_ideas_matrix.reshape([data_ideas_matrix.shape[0],data_ideas_matrix.shape[1],1,1])

print('data_ideas_matrix.shape')
print(data_ideas_matrix.shape)
### read rna-seq
data_rna = read_rna2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/rsem/gmp.rsem.txt',1,4)

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

y_test = data_rna[pos_test_id,:]
y_train = data_rna[pos_train_id,:]


cell_types = ['er4']
cell_types_num = len(cell_types)+1
for ct in cell_types:
	### read atac_seq pk
	data_atac = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn_atacpk/atac_21k_binary_split.tab',1,6)
	data_atac_cellspe = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn_atacpk/atac_21k_binary_split.er4.tab',1,6)

	### read ideas labels
	data_ideas_matrix = data_atac
	#data_ideas_matrix = data_ideas_matrix.reshape([data_ideas_matrix.shape[0],data_ideas_matrix.shape[1],1,1])

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


i=0
# Create linear regression object
regr = linear_model.LinearRegression()
print('training linear model')
# Train the model using the training sets
regr.fit(x_train, y_train)
print('training linear model Done!')
# Make predictions using the testing set
y_p = regr.predict(x_train)
y_obs=y_train[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_train=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, LM TRAINING accuracy same cell!!! R2: %g"%(i, r2_train)  )

y_p = regr.predict(x_test)
y_obs=y_test[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_test=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, LM TESTING accuracy same cell!!! R2: %g"%(i, r2_test)  )

'''
training model Done!
-2.78601956485
step 16, TRAINING accuracy same cell!!! R2: 0.785169
-2.8011434812
step 16, TESTING accuracy same cell!!! R2: 0.744683

real	9m36.230s
user	6m22.094s
sys	1m18.330s
'''

max_depth = 30
regr_multirf = MultiOutputRegressor(RandomForestRegressor(max_depth=max_depth,random_state=0))
print('training RF model')
regr_multirf.fit(x_train, y_train)
print('training RF model Done!')


y_p = regr_multirf.predict(x_train)
y_obs=y_train[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_train=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, RF TRAINING accuracy same cell!!! R2: %g"%(i, r2_train)  )

y_p = regr_multirf.predict(x_test)
y_obs=y_test[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_test=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, RF TESTING accuracy same cell!!! R2: %g"%(i, r2_test)  )

'''
-2.78601956485
step 16, RF TRAINING accuracy same cell!!! R2: 0.945249
-2.8011434812
step 16, RF TESTING accuracy same cell!!! R2: 0.787583

real	13m12.701s
user	10m18.736s
sys	0m45.913s
'''

svr = GridSearchCV(SVR(kernel='rbf', gamma=0.1), cv=5,param_grid={"C": [1e0, 1e1, 1e2, 1e3],"gamma": np.logspace(-2, 2, 5)})
print('training SVR model')
svr.fit(x_train, y_train)
print('training SVR model Done!')

y_p = svr.predict(x_train)
y_obs=y_train[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_train=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, SVR TRAINING accuracy same cell!!! R2: %g"%(i, r2_train)  )

y_p = svr.predict(x_test)
y_obs=y_test[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_test=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, SVR TESTING accuracy same cell!!! R2: %g"%(i, r2_test)  )



kr = GridSearchCV(KernelRidge(kernel='rbf', gamma=0.1), cv=5,param_grid={"alpha": [1e0, 0.1, 1e-2, 1e-3],"gamma": np.logspace(-2, 2, 5)})
print('training KR model')
kr.fit(x_train, y_train)
print('training KR model Done!')

y_p = kr.predict(x_train)
y_obs=y_train[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_train=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, KR TRAINING accuracy same cell!!! R2: %g"%(i, r2_train)  )

y_p = kr.predict(x_test)
y_obs=y_test[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_test=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, KR TESTING accuracy same cell!!! R2: %g"%(i, r2_test)  )

