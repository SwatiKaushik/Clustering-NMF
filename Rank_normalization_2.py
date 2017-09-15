__author__ = "swati"


#Performs rank normalization of the given gene expression matrix
# Sep 14 2017

import pandas as pd
import operator
import argparse

def column(matrix, i):
	
	"""
	returns the required column
	:param matrix: input matrix, pandas dataframe
	:param i: column number
	"""
    return [row[i] for row in matrix]
  
def ranknormalization(table, gcount):
	
	"""
	Performs sample wise/row wise rank normalization 
	:param table: input matrix, pandas dataframe
	:param gcount: initialize as zero
	"""
	
   for index,row in table.iterrows():
        count = len(row)
        adict = {}
   
        for i in range(count):
           adict[colnames[i]] = row[i]

        sorted_x = sorted(adict.items(), key=lambda t: t[1])
        #print (sorted_x) 
    
        count = 0
        newdict ={}
    
        for i in range(len(sorted_x)):
          count += 1
          #print (sorted_x[i][0], count)
          newdict[sorted_x[i][0]] = count/lenc
    
        index_map = {v: i for i, v in enumerate(colnames)}
        colone = sorted(newdict.items(), key=lambda pair: index_map[pair[0]])
        sorted_col = column(colone, 1)
    
        sortedvalues = ('\t'.join([str(x) for x in sorted_col]))
        output_row = rownames[gcount] +"\t" + sortedvalues
        #print(rownames[gcount], sorted_col)
        #print (output_row)
        test = output_row + "\n"
        outputfile.writelines(test)
        gcount += 1
  
parser = argparse.ArgumentParser()                                               

parser.add_argument("--file", "-f", type=str, required=True)
args = parser.parse_args()
    
table = pd.read_table(args.file, sep='\s+')
#raw_data_file

colnames= list(table.columns)
lenc = len(colnames)
rownames = list(table.index)
gcount =0

#write output file
outputfile = open('raw_data_file_rank_normalized.out', 'w')
outputfile.writelines(colnames)
outputfile.writelines("\n")

ranknormalization(table, gcount)
    
outputfile.close()

