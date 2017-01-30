#!/bin/bash


PROJDIR='/users/mtaranov/3D_fromATAC_by_chr'
DATADIR=$PROJDIR/'data'
LABELDIR=$PROJDIR/'labels'
JOBSDIR=$PROJDIR/'jobs'

mkdir $LABELDIR

PromoterCapture_file=$DATADIR/'PromoterCapture_Digest_Human_HindIII_baits_ID.bed'
# original was sorted with sort -k 1,1 -k2,2n
atacPeak_D0_file=$DATADIR/'SORTED_primary_keratinocyte-d00.GGR.Stanford_Greenleaf.ATAC-seq.b1.trim.PE2SE.nodup.tn5_pooled.pf.pval0.1.500000.naive_overlap.narrowPeak.gz'
CaptureC_bait_bait_D0_file=$DATADIR/'D0_D2D8_merge_BaitToBait_intra.bed.gz'
atac_1kb_bins_D0_file=$DATADIR/'atac_1kb_bins_D0.bed'

   
/users/mtaranov/local/anaconda2/bin/python $PROJDIR/get1_1kb_atac_to_Hind3_mapping.py $PromoterCapture_file $atacPeak_D0_file $CaptureC_bait_bait_D0_file $atac_1kb_bins_D0_file_chr 
