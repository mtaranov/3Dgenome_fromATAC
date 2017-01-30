#!/bin/bash


PROJDIR='/users/mtaranov/3D_fromATAC_by_chr'
DATADIR=$PROJDIR/'data'
LABELDIR=$PROJDIR/'labels'
JOBSDIR=$PROJDIR/'jobs'


#thres=10

#rm $LABELDIR/chr*_atac_1kb_bins.bed
atac_1kb_contacts_w_labels=$LABELDIR/all_atac_1kb_contacts_w_labels.txt
#cat $LABELDIR/chr*_atac_1kb_contacts_w_labels.bed > $atac_1kb_contacts_w_labels


python $PROJDIR/get3a_all_pairs.py $atac_1kb_contacts_w_labels $LABELDIR

#gzip $LABELDIR/*w_labels*
