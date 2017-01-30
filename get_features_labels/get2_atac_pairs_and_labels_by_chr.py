#!/usr/bin/env python

import sys
from bed_manipulate import get_HindIII_ContactPairs, get_2D_atac_labels

PromoterCapture_file = sys.argv[1]
atacPeak_D0_file = sys.argv[2]
CaptureC_bait_bait_D0_file = sys.argv[3]
atac_1kb_bins_D0_file = sys.argv[4]
atac_1kb_contacts_w_labels_D0_file = sys.argv[5]

# STEP 1
# returns dict with key=chr and val= dict with key = list of HindIII contact pairs for each chr and val=score
HindIII_pairs_D0 =get_HindIII_ContactPairs(CaptureC_bait_bait_D0_file)

# STEP 2
# assigns 1/0 labels to atac_bin interactions
# outputs atac_1kb_contacts_w_labels_D0_file for each chr
get_2D_atac_labels(atac_1kb_bins_D0_file, HindIII_pairs_D0, atac_1kb_contacts_w_labels_D0_file)
###label_list=get_2D_atac_labels(atac_1kb_bins_D0_file, HindIII_pairs_D0, atac_1kb_contacts_w_labels_D0_file)
###print "pos = ", len(filter(lambda x: x is not '0', label_list)), " neg = ", len(filter(lambda x: x is '0', label_list))

