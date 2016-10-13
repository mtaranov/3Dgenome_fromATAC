import numpy as np
import pybedtools
from pybedtools import BedTool

def bed_intersection(bed1_in, bed2_in, bed_out=False, bed_out_name='',  f=0, F=0, e=True):
    """
    intersects bed1 with bed2 
    """
  
    bed1 = BedTool(bed1_in)
    bed2 = BedTool(bed2_in)
    
    if bed_out == True:
        bed1.intersect(bed2_in, c=True, f=f, F=F, e=e).saveas(bed_out_name)
    return bed1.intersect(bed2_in, c=True, f=f, F=F, e=e)  

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


