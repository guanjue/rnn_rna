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


sec_d=50 ### upstream_expand_num+gene_len_num
sec_d_rev=50+50-1 ### downstream_expand_num-1 
genebody_len = 50
sec_d_genebody = sec_d+genebody_len

all_len=50+50+50
thr_d=1
for_d=17

exp_win = 1
### read atac_seq pk
data_atac = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn_atacpk/atac_21k_binary_split.tab',1,6)
### read ideas labels
data_ideas_matrix = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn/matrix_for_rnn/gene_list.gmp.0.tab',1,6)
#data_ideas_matrix = data_ideas_matrix * data_atac
for i in range(1,17):
	print(i)
	data_ideas_tmp = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn/matrix_for_rnn/gene_list.gmp.'+str(i)+'.tab',1,6)
	#data_ideas_tmp = data_ideas_tmp * data_atac
	data_ideas_tmp = np.concatenate((data_ideas_tmp[:,0:sec_d-exp_win], data_ideas_tmp[:,sec_d_genebody+exp_win:all_len]), axis=1)
	data_ideas_matrix = np.concatenate((data_ideas_matrix, data_ideas_tmp), axis=1)

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
	### read ideas labels
	data_ideas_matrix = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn/matrix_for_rnn/gene_list.'+ct+'.0.tab',1,6)
	#data_ideas_matrix = data_ideas_matrix * data_atac
	for i in range(1,17):
		print(i)
		data_ideas_tmp = read_tab2matrix('/Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/matrix_for_rnn/matrix_for_rnn/gene_list.'+ct+'.'+str(i)+'.tab',1,6)
		#data_ideas_tmp = data_ideas_tmp * data_atac
		data_ideas_tmp = np.concatenate((data_ideas_tmp[:,0:sec_d-exp_win], data_ideas_tmp[:,sec_d_genebody+exp_win:all_len]), axis=1)
		data_ideas_matrix = np.concatenate((data_ideas_matrix, data_ideas_tmp), axis=1)

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

print('LM stats.pearsonr(y_obs,y_p)')
print(stats.spearmanr(y_obs,y_p))
print(stats.pearsonr(y_obs,y_p))

'''
step 16, LM TRAINING accuracy same cell!!! R2: 0.489922
-2.8011434812
step 16, LM TESTING accuracy same cell!!! R2: 0.415026
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

print('RF stats.pearsonr(y_obs,y_p)')
print(stats.spearmanr(y_obs,y_p))
print(stats.pearsonr(y_obs,y_p))

'''
training RF model
training RF model Done!
-2.78601956485
step 16, RF TRAINING accuracy same cell!!! R2: 0.87104
-2.8011434812
step 16, RF TESTING accuracy same cell!!! R2: 0.461647
'''
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
print("step %d, KR TRAINING accuracy same cell!!! R2: %g"%(i, r2_train))


y_p = kr.predict(x_test)
y_obs=y_test[:,:]
y_mean=np.mean(y_obs)
print(y_mean)
r2_test=1-np.sum(np.square(y_obs-y_p))/np.sum(np.square(y_obs-y_mean))
print("step %d, KR TESTING accuracy same cell!!! R2: %g"%(i, r2_test)  )
'''
