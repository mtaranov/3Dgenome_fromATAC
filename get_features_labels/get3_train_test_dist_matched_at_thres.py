#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from bed_manipulate import get_contacts_at_distance, get_pairs_distance_matched
from utils_atac import binarize
from sklearn.model_selection import train_test_split

atac_1kb_contacts_w_labels_file = sys.argv[1]
min_dist = int(sys.argv[2])
max_dist = int(sys.argv[3])
dist_step = int(sys.argv[4])
thres=int(sys.argv[5])
#output dir
label_dir = sys.argv[6]

imbalance_ratio = 1

#### Get Contacts at Distance w/o Distance matching ########
#thres_min, thres_max =10000, 2000000
#get_contacts_at_distance(atac_1kb_contacts_w_labels_D0_file, thres_min, thres_max, atac_1kb_contacts_w_labels_btw_10kb_2Mb_D0_file)

### Load data/labels/indx ##############################

data=np.loadtxt(atac_1kb_contacts_w_labels_file, dtype=str)
X=data[:,:6]
y=data[:,6]
y=y.reshape((y.shape[0],1))
indx=data[:,7:9]

#### Split Train/Test ##################################

X_train, X_test, indx_train, indx_test, y_train, y_test = train_test_split(X, indx, y,  test_size=0.3, random_state=42)

#### Get Distance Matched Contacts at threshold ########
print "Getting Distance Matched Contacts in Train ..."
X_train_new, y_train_new, indx_train_new = get_pairs_distance_matched(X_train, y_train, indx_train, min_dist, max_dist, dist_step, imbalance_ratio, thres)
print "Getting Distance Matched Contacts in Test ..."
X_test_new, y_test_new, indx_test_new = get_pairs_distance_matched(X_test, y_test, indx_test, min_dist, max_dist, dist_step, imbalance_ratio, thres)

print 'train: ', X_train_new.shape, y_train_new.shape, indx_train_new.shape
print 'test: ', X_test_new.shape, y_test_new.shape, indx_test_new.shape


np.savetxt(label_dir+'/site1_train_dist_matched_thres10.bed', X_train_new[:,:3], delimiter="\t", fmt="%s")
np.savetxt(label_dir+'/site1_test_dist_matched_thres10.bed', X_test_new[:,:3], delimiter="\t", fmt="%s")
np.savetxt(label_dir+'/site2_train_dist_matched_thres10.bed', X_train_new[:,3:], delimiter="\t", fmt="%s")
np.savetxt(label_dir+'/site2_test_dist_matched_thres10.bed', X_test_new[:,3:], delimiter="\t", fmt="%s")
np.savetxt(label_dir+'/indx_train_dist_matched_thres10.bed', indx_train_new, delimiter="\t", fmt="%s")
np.savetxt(label_dir+'/indx_test_dist_matched_thres10.bed', indx_test_new, delimiter="\t", fmt="%s")

np.save(label_dir+'/labels_train_dist_matched_thres10.npy', binarize(y_train_new.astype(float)).astype(int))
np.save(label_dir+'/labels_test_dist_matched_thres10.npy', binarize(y_test_new.astype(float)).astype(int))
