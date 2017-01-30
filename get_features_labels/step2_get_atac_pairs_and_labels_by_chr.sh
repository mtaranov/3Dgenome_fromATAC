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


mkdir $JOBSDIR
chrs=`cat /users/mtaranov/3D_fromATAC_by_chr/data/atac_1kb_bins_D0.bed | awk '{print $1}' | uniq`
for chr in $chrs
#for chr in "chrY"
do
        atac_1kb_bins_D0_file_chr=$LABELDIR/${chr}_atac_1kb_bins.bed
        atac_1kb_contacts_w_labels_D0_file_chr=$LABELDIR/${chr}_atac_1kb_contacts_w_labels.bed
        cat  $atac_1kb_bins_D0_file | awk '{if ($1==chr) print $0}' chr=$chr > $atac_1kb_bins_D0_file_chr
   
        command1="/users/mtaranov/local/anaconda2/bin/python $PROJDIR/get2_atac_pairs_and_labels_by_chr.py $PromoterCapture_file $atacPeak_D0_file $CaptureC_bait_bait_D0_file $atac_1kb_bins_D0_file_chr $atac_1kb_contacts_w_labels_D0_file_chr"
        jobfile=$JOBSDIR/${chr}.job
        echo "$command1" >> $jobfile
        chmod 777 ${jobfile}
        qsub -q q@nandi  -o $JOBSDIR/o.${chr} -e $JOBSDIR/e.${chr} ${jobfile}

done

