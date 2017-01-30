#!/bin/bash


PROJDIR='/users/mtaranov/3D_fromATAC_by_chr'
DATADIR=$PROJDIR/'data'
LABELDIR=$PROJDIR/'labels'
JOBSDIR=$PROJDIR/'jobs'


min_dist=10000
max_dist=2000000
dist_step=10000
#imbalance_ratio=1
thres=10

#rm $LABELDIR/chr*_atac_1kb_bins.bed
atac_1kb_contacts_w_labels=$LABELDIR/all_atac_1kb_contacts_w_labels.txt
#cat $LABELDIR/chr*_atac_1kb_contacts_w_labels.bed > $atac_1kb_contacts_w_labels


python $PROJDIR/get3_train_test_dist_matched_at_thres.py $atac_1kb_contacts_w_labels $min_dist $max_dist $dist_step  $thres $LABELDIR

#gzip $LABELDIR/*w_labels*
