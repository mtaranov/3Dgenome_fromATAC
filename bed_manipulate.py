import numpy as np
import pybedtools
from pybedtools import BedTool
import csv

def bed_intersection(bed1_in, bed2_in, bed_out=False, bed_out_name='',  f=0, F=0, e=True):
    """
    intersects bed1 with bed2 
    """
  
    bed1 = BedTool(bed1_in)
    bed2 = BedTool(bed2_in)
    
    if bed_out == True:
        bed1.intersect(bed2_in, c=True, f=f, F=F, e=e).saveas(bed_out_name)
    return bed1.intersect(bed2_in, c=True, f=f, F=F, e=e)  

def bed(bed1_in, bed2_in, bed_out=False, bed_out_name='',  f=0, F=0, e=True, wo=True):
    """
    intersects bed1 with bed2 
    """
  
    bed1 = BedTool(bed1_in)
    bed2 = BedTool(bed2_in)
    
    if bed_out == True:
        bed1.intersect(bed2_in, f=f, F=F, e=e, wo=wo).saveas(bed_out_name)
    return bed1.intersect(bed2_in, f=f, F=F, e=e, wo=wo)   

def featuretype_filter(feature, featuretype):
    """
    returns lines which match featuretype 
    """
    if int(feature[5]) == featuretype:
        return True
    return False

def bed_w_NoIntersection(bed1_in, bed2_in,  bed_out_NoInters, bed_out=False, bed_out_name='',featuretype=0, f=0, F=0, e=True):
    """
    filters based on the featuretype_filter
    """
    a=bed_intersection(bed1_in, bed2_in, bed_out=bed_out, f=f, F=F, e=e)
    result=a.filter(featuretype_filter, featuretype).saveas(bed_out_NoInters)
    return pybedtools.BedTool(result.fn)


def bed_closest(bed1_in, bed2_in, bed_out=False, bed_out_name='', io=True):
    """
    finds closest distance
    """

    bed1 = BedTool(bed1_in)
    bed2 = BedTool(bed2_in)

    if bed_out == True:
        bed1.closest(bed2, io=io).saveas(bed_out_name)
    return bed1.closest(bed2, io=io)

# pulls bait-bait contacts from all-all capture  contacts
def get_bait_bait_contacts(InteractionsFile_D0, PromoterCaptureFile, outFile, site1, site2):
    data = BedTool(InteractionsFile_D0)
    BedTool([i[0], i[1], i[2], i[6]] for i in data).saveas(site1)
    BedTool([i[3], i[4], i[5], i[6]] for i in data).saveas(site2)
    bed1_site1=BedTool(site1)
    bed1_site2=BedTool(site2)
    bed2 = BedTool(PromoterCaptureFile)

    with open(outFile, 'w') as outcsv:
        #configure writer to write standart csv file
            writer = csv.writer(outcsv, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
            for interval in zip(bed1_site1.intersect(bed2, wao=True), bed1_site2.intersect(bed2, wao=True)):
                if interval[0][0] == interval[1][0]: # intra-chr contacts only 
                    if (int(interval[0][9]) != 0 and int(interval[1][9]) != 0):
                        writer.writerow([interval[0][0], interval[0][1], interval[0][2], interval[1][0], interval[1][1], interval[1][2], interval[1][3], interval[0][7], interval[1][7]])

# intersects atac-peak file with promoter-capture file
# assigns each peak summit to a HindIII
# outputs 1kb regions aroud peak summit and corresponding HindIII
def get_atac_bins_at_1kb(bed1_in, bed2_in, outFile, bed_out=False, bed_out_name='',  f=0.0, F=0.0, e=True, wo=True):

    bed1 = BedTool(bed1_in)
    bed2 = BedTool(bed2_in)

    if bed_out == True:
        bed1.intersect(bed2, f=f, F=F, e=e, wo=wo).saveas(bed_out_name)
    
    atac_HindIII={}
    # if one ATAC-fragment/peak overlaps multiple HindIII fragments, keep ATAC-fragment w the largest overlap
    for interval in bed1.intersect(bed2, f=f, F=F, e=e, wo=wo):
        if interval[0] not in atac_HindIII:
            atac_HindIII[interval[0]]={interval[1] : {interval[9] : [interval[13], interval[15]]}}
        else:
            if interval[1] not in atac_HindIII[interval[0]]:
                atac_HindIII[interval[0]][interval[1]]={interval[9] : [interval[13], interval[15]]}
            else:
                if interval[9] not in atac_HindIII[interval[0]][interval[1]]:
                     atac_HindIII[interval[0]][interval[1]][interval[9]]=[interval[13], interval[15]]
                else:
                     if interval[15] > atac_HindIII[interval[0]][interval[1]][interval[9]][1]:
                        atac_HindIII[interval[0]][interval[1]].update({interval[9] : [interval[13], interval[15]]})          
                    
    with open(outFile, 'w') as outcsv:
    #configure writer to write standart csv file
        writer = csv.writer(outcsv, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
        for chr in atac_HindIII.keys():
            for start_pos in atac_HindIII[chr].keys():
                for peaksummit in atac_HindIII[chr][start_pos].keys():
                    
                    writer.writerow([str(chr), str(int(start_pos) + int(peaksummit) - 500),  str(int(start_pos) + int(peaksummit) + 500),  str(atac_HindIII[chr][start_pos][peaksummit][0])])
                
# get the list of HindIII contact pairs for each chr
def get_HindIII_ContactPairs(bait_bait_file):
    contact_dict={}
    for line in open(bait_bait_file,'r'):
        words=line.rstrip().split()
        chr1=words[0]
        chr2=words[3]
        q_value=words[6]
        HindIII1=words[7]
        HindIII2=words[8]
        if chr1 != chr2:
            raise ValueError("Inter-chr contact!")
            exit
        else:
            if float(q_value) > 0: # only confident contacts
                if chr1 not in contact_dict:
                    if int(HindIII1) < int(HindIII2):
                        contact_dict[chr1]=[(HindIII1, HindIII2)]
                    else:
                        contact_dict[chr1]=[(HindIII2, HindIII1)]
                else:
                    if int(HindIII1) < int(HindIII2):
                        contact_dict[chr1].append((HindIII1, HindIII2))
                    else:
                        contact_dict[chr1].append((HindIII2, HindIII1))
    return contact_dict

# assigns 1/0 labels to atac_bin interactions
def get_2D_atac_labels(atac_bins_file, pos_HindIII_pairs, outFile):
    atac_dict={}
    label_list=[]
    
    for line in open(atac_bins_file,'r'):
        words=line.rstrip().split()
        chr=words[0]
        start=words[1]
        end=words[2]
        HindIII=words[3]
        if chr not in atac_dict:
            atac_dict[chr]=[(start, end, HindIII)]
        else:
            atac_dict[chr].append((start, end, HindIII))

    with open(outFile, 'w') as outcsv:
    #configure writer to write standart csv file
        writer = csv.writer(outcsv, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
        label_list=[]
        for chr in atac_dict:
            #if chr=='chrX':
            for item1 in atac_dict[chr]:
                for item2 in atac_dict[chr]:   
                    if int(item1[0]) < int(item2[0]):
                        if ((item1[2],item2[2]) in pos_HindIII_pairs[chr]):
                            label=1
                        else:
                            label=0
                        label_list.append(label)
                        writer.writerow([chr, item1[0],item1[1], chr, item2[0],item2[1], label, item1[2],item2[2]])
    return label_list

def get_contacts_at_distance(atac_w_labels_file, thres_min, thres_max, outFile):
    atac_dict={}
    for line in open(atac_w_labels_file,'r'):
        words=line.rstrip().split()
        chr=words[0]
        start1=words[1]
        end1=words[2]
        start2=words[4]
        end2=words[5]
        label=words[6]
        HindIII1=words[7]
        HindIII2=words[8]
        if (abs(int(start2)-int(start1)) > thres_min and abs(int(start2)-int(start1)) < thres_max):
            if chr not in atac_dict:
                atac_dict[chr]=[(start1, end1, start2, end2, label, HindIII1, HindIII2)]
            else:
                atac_dict[chr].append((start1, end1, start2, end2, label, HindIII1, HindIII2))

    with open(outFile, 'w') as outcsv:
    #configure writer to write standart csv file
        writer = csv.writer(outcsv, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
        for chr in atac_dict:
                for item in atac_dict[chr]:
                    writer.writerow([chr, item[0], item[1], chr, item[2], item[3], item[4], item[5], item[6]])

