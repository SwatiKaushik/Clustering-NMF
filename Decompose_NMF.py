__author__ = "swati"

# script to generate n number of NMF components from the given input matrix and generate a heatmap
# Sep 12 2017


import argparse
import os

import pandas as pd
import seaborn as sns
from sklearn.decomposition import NMF
from scipy.stats import zscore

def dfConcatenate(table, annotation):
    """
    This function concatenates two dataframes
    :param table: input gene expression matrix, pandas dataframe, (n_features, m_samples)
    :param annotation: sample annotation file, (n_annotation, m_samples)
    :return concatenated dataframe
    """

    table = table.sort_index(axis=1)
    annotation = annotation.sort_index(axis=1)
    concat = pd.concat([table,annotation])
    return concat

def Calczscore(filename):
    """
    This function calculates the Z score of the dataframe
    :param filename: pandas dataframe
    :return z scored dataframe
    """
    
    filenamet= filename.transpose()
    filenamezscore = filenamet.apply(zscore)
    out = filenamezscore.transpose()
    return out
    
def plotheatmap(filename, k, annotation):

	"""
	This function generates heatmap of h matrix 
	:param filename: h matrix, pandas dataframe, (n_features, m_samples)
	:param k: # of NMF components 
	:param annotation: sample annotation file, (n_annotation, m_samples)
	:output: save heatmap of h matrix
	"""

	table = pd.read_table(filename)
	
	#sort, merge dataframe and add annotation
	merged = dfConcatenate(table, annotation)
	#remove the sample annotation
	merged.drop(merged.columns[[-1,]], axis=1, inplace=True)
	#remove last row
	names = merged.iloc[-1]
	
	#generate heatmap
	annotation = merged.drop(['annotation'])
	annotation = annotation.astype(float)
	out = Calczscore(annotation)
	lut = dict(zip(names.unique(), "rbg"))
	col_colors = pd.Series(names, name='clusters').map(lut)
	sns.plot = sns.clustermap(out, cmap="coolwarm", col_colors= col_colors, robust=True)

	#save heatmap
	outplot = filename + ".png"
	sns.plot.savefig(outplot)
    
def decomposeNMF(matrix,clusterannotationfile, components):
	
	"""
	This function decomposes given gene expression matrix into multiple components 
	:param matrix: input gene expression matrix, pandas dataframe, (n_features, m_samples)
	:param clusterannotationfile: sample annotation file, (n_annotation, m_samples)
    :output: generate h and w matrix
	"""
	
	for k in range(2,components):
		
		model = NMF(n_components=k, init='random', random_state=0)
		w, h = model.fit_transform(matrix), model.components_
		w = pd.DataFrame(w, index= matrix.index)
		h = pd.DataFrame(h, columns= matrix.columns)

		# get weight matrix
		outfile_w = "wmatrix" + "_"+ str(k)
		w.to_csv(outfile_w, sep='\t', encoding='utf-8')

		# get feature matrix
		outfile_h = "hmatrix" + "_"+ str(k)
		h.to_csv(outfile_h, sep='\t', encoding='utf-8')
		
		plotheatmap(outfile_h, k, clusterannotationfile)


parser = argparse.ArgumentParser()                                               

parser.add_argument("--file", "-f", type=str, required=True)
args = parser.parse_args()

matrix = pd.read_table(args.file)
clusterannotationfile = pd.read_table("test_annotation.txt")
components = len(matrix.index)-1
decomposeNMF(matrix, clusterannotationfile, components)


