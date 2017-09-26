# Use CNV data to annotate the transcriptional states

import pandas as pd
import os
import numpy as np
from scipy import stats
from math import *

def common_elements(list1, list2):

    """
    Calculates common elements between two list
    :param list1: list1
    :param list2: list2
    :return common elements
    """
    return [element for element in list1 if element in list2]

def create_newhmatrix(hmatrix, cnv_file):
	
	"""
    Calculate new matrix based on common elements
    :param hmatrix: matrix decomposition
    :param cnv_file: dataframe of cnv file
    :return new_common elemen matrix
    """
    
    # make list of column names
    hcol = list(hmatrix)
    Rcol = list(cnv_file)
    
    #Identify common elements
    common_cols = list(common_elements(hcol,Rcol))
    
    # extract common columns 
    new_hmatrix = hmatrix[common_cols]
    
    return new_hmatrix

def cnv_correlation (hmatrix, cnv_file):

	"""
    Calculate correlation of cnv data with hmatrix
    :param hmatrix: matrix decomposition
    :param cnv_file: dataframe of cnv file
    :return output file with top 10 associations
    """

	  # length of hmatrix : for every component	
	  
      for i in range(len(hmatrix)):
            dict = {}
    
            # for every gene in cnv file
            for y in range (len(cnv_file)):
                
                 # calculate spearman correlation and p value
                 corr, p = stats.spearmanr(hmatrix.iloc[i], cnv_file.iloc[y])     
                 corround = "%.4f" % corr

                 # Generate a dictionary
                 dict[cnv_file.index[y]]=  corround, p
    
            #print (list(dict.values()))
            
            #sort data and get top 10 associations
            c = sorted(dict.items(), key=lambda value: value[1][1])
            test = str(i) + "\t" + str(c[1:10]), "\n"
     
            # write output file
            outputfile.writelines(test)
            print (i)
            
hmatrix = pd.read_table("hmatrix_9",sep='\s+')
cnv_file = pd.read_table("cnv_gistic.txt",sep='\s+')

new_hmatrix = create_newhmatrix(hmatrix, cnv_file )
new_cnvmatrix = create_newhmatrix(cnv_file, new_hmatrix )
outputfile = open('cnv_9.out', 'w')
cnv_correlation(new_hmatrix, new_cnvmatrix)

