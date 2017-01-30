#!/usr/bin/env python

import sys
#from bed_manipulate import get_bait_bait_contacts, get_atac_bins_at_1kb, get_HindIII_ContactPairs, get_2D_atac_labels, get_contacts_at_distance
#from bed_manipulate_exp import get_bait_bait_contacts, get_atac_bins_at_1kb, get_HindIII_ContactPairs, get_2D_atac_labels, get_contacts_at_distance
from bed_manipulate import  get_atac_bins_at_1kb

PromoterCapture_file = sys.argv[1]
atacPeak_D0_file = sys.argv[2]
CaptureC_bait_bait_D0_file = sys.argv[3]
atac_1kb_bins_D0_file = sys.argv[4]

# don't need for Bait-Bait contacts since they are labeled in Rds file
#CaptureContacts_D0_file='/mnt/lab_data/kundaje/mtaranov/ChicagoCalls/superConfident_bait_bait_fromAdam/humKer_Sub_D2D8_merge.bed'
#site1_D0_file=DATADIR+'site1_D0.tmp'
#site2_D0_file=DATADIR+'site2_D0.tmp'
#CaptureC_bait_bait_D0_file=DATADIR+'CaptureC_bait-bait_D0.bed'

# STEP1
# intersects atac-peak file with promoter-capture file
# assigns each peak summit to a HindIII
# returns file with 1kb regions aroud peak summit and corresponding HindIII number - atac_1kb_bins_D0_file
get_atac_bins_at_1kb(atacPeak_D0_file, PromoterCapture_file, atac_1kb_bins_D0_file, bed_out=False, bed_out_name='',  f=0.5, F=0.5, e=True, wo=True)

## STEP 2
## don't need for Bait-Bait contacts since they are labeled in Rds file
#
## pulls bait-bait contacts from all-all capture contacts
## these are positive lables
## outputs CaptureC_bait_bait_D0_file
#get_bait_bait_contacts(CaptureContacts_D0_file, PromoterCapture_file, CaptureC_bait_bait_D0_file, site1_D0_file, site2_D0_file)
