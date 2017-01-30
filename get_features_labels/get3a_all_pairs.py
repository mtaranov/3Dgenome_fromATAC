#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from bed_manipulate import get_contacts_at_distance, get_pairs_distance_matched
from utils_atac import binarize
from sklearn.model_selection import train_test_split

atac_1kb_contacts_w_labels_file = sys.argv[1]
label_dir = sys.argv[2]
thres=10.0
#output dir

data=np.loadtxt(atac_1kb_contacts_w_labels_file, dtype=str)
X=data[:,:6]
y=data[:,6]
y=y.reshape((y.shape[0],1))
indx=data[:,7:9]

neg_indxs = np.where(y.astype(float)==0)[0]
pos_indxs = np.where(y.astype(float)>=thres)[0]
X_pos=X[pos_indxs]
X_neg=X[neg_indxs]
y_pos=y[pos_indxs]
y_neg=y[neg_indxs]
indx_pos=indx[pos_indxs]
indx_neg=indx[neg_indxs]

y_new=np.concatenate((y_pos, y_neg))
X_new=np.concatenate((X_pos, X_neg))
indx_new=np.concatenate((indx_pos, indx_neg))

X_train, X_test, indx_train, indx_test, y_train, y_test = train_test_split(X_new, indx_new, y_new,  test_size=0.1, random_state=42)

print 'train: ', X_train.shape, y_train.shape, indx_train.shape
print 'test: ', X_test.shape, y_test.shape, indx_test.shape


np.savetxt(label_dir+'/site1_train_all_pairs_thres10.bed', X_train[:,:3], delimiter="\t", fmt="%s")
np.savetxt(label_dir+'/site1_test_all_pairs_thres10.bed', X_test[:,:3], delimiter="\t", fmt="%s")
np.savetxt(label_dir+'/site2_train_all_pairs_thres10.bed', X_train[:,3:], delimiter="\t", fmt="%s")
np.savetxt(label_dir+'/site2_test_all_pairs_thres10.bed', X_test[:,3:], delimiter="\t", fmt="%s")
np.savetxt(label_dir+'/indx_train_all_pairs_thres10.bed', indx_train, delimiter="\t", fmt="%s")
np.savetxt(label_dir+'/indx_test_all_pairs_thres10.bed', indx_test, delimiter="\t", fmt="%s")

np.save(label_dir+'/labels_train_all_pairs_thres10.npy', binarize(y_train.astype(float)).astype(int))
np.save(label_dir+'/labels_test_all_pairs_thres10.npy', binarize(y_test.astype(float)).astype(int))
