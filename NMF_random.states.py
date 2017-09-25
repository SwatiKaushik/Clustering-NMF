import os
from numpy import asarray, zeros, argmax
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from seaborn import heatmap
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, cophenet
from sklearn.decomposition import NMF
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr
from sklearn.cluster import AgglomerativeClustering


def nmf(matrix, k, random_seed):
    """
    Performs nmf on the given matrix
    :param matrix: input gene expression matrix, pandas dataframe, (n_features, m_samples)
    :param k: number of decompositions
    :param random_seed: Random seed number
    :return: Dictionary of NMF results 
    """
    
    nmf_results = {}

    #Generate NMF components
    model = NMF(n_components=k, init='random', random_state=random_seed + i)
    w, h = model.fit_transform(matrix), model.components_
   
   	#Save w/h matrices
    w = pd.DataFrame(w, index=matrix.index)
    h = pd.DataFrame(h, columns=matrix.columns)
    nmf_results[k] = {'w': w, 'h': h}  
    
    return nmf_results

def plotnmf (nmf_result):
    """
    Plot results of nmf (dictionary as input)
    :param nmf_result: generates heatmap pf h/w matrix
    """
    w_matrix = nmf_result
    ax = plt.axes()
    
    #select the w/h matrix
    heatmap(w_matrix, cmap="RdBu_r",ax=ax)
    ax.set_ylabel('Feature')
    ax.set_xlabel('Component')
    ax.set_title(k)
    #plt.show()

def dfConcatenate(table, annotation): 
    """
    Concatenates two dataframes
    :param table: input gene expression matrix, pandas dataframe, (n_features, m_samples)
    :param annotation: sample annotation file, (n_annotation, m_samples)
    :return concatenated dataframe
    """
    table = table.sort_index(axis=1)
    annotation = annotation.sort_index(axis=1)
    
    #concatenate datasets
    concat = pd.concat([table,annotation])
    
    return concat

def plot_nmf(table, annotation,  outfile_h):
    """
    Generates heatmap of h matrix 
    :param table: h matrix, pandas dataframe, (n_features, m_samples)
	:param annotation: sample annotation file, (n_annotation, m_samples)
	:param outputfile_h: Name of output file
	:output: save heatmap of h matrix

    """
    merged = dfConcatenate(table, annotation)
    merged.drop(merged.columns[[-1,]], axis=1, inplace=True)
    names = merged.iloc[-1]
    annotation = merged.drop(['annotation'])
    annotation = annotation.astype(float)
  
    lut = dict(zip(names.unique(), "rbg"))
    #print(lut)
    col_colors = pd.Series(names, name='clusters').map(lut)
  
    sns.plot = sns.clustermap(annotation, cmap="coolwarm", col_colors= col_colors, robust=True, z_score =0)
    outplot = outfile_h + ".png"
    sns.plot.savefig(outplot)
        
def plot_coph(cophen):
    """
    Plot cophenetic correlation coefficient of NMF consensus cluster
    :param cophen: Dictionary of cophenetic correlation coefficients
    """
    fig = plt.figure()
    
    # Get cophenetic values per k
    
    A=(list(cophen.keys()),list(cophen.values()))

    fig.suptitle('Cophenetic correlation coefficient vs. K', fontsize=20)
    plt.scatter(A[0], A[1])
    plt.plot(A[0], A[1])
    plt.xlabel('k', fontsize=18)
        
    plt.ylabel('Cophenetic correlation coefficient', fontsize=16)
    fig.savefig('CCC_K.jpg')
        
def _get_consensus(sample_x_clustering):
    """
    Count number of co-clusterings.
    :param sample_x_clustering: DataFrame; (n_samples, n_clusterings)
    :return: DataFrame; (n_samples, n_samples)
    """

    sample_x_clustering_array = asarray(sample_x_clustering)

    n_samples, n_clusterings = sample_x_clustering_array.shape

    # Make sample x sample matrix
    coclusterings = zeros((n_samples, n_samples))

    # Count the number of co-clusterings
    for i in range(n_samples):
        for j in range(n_samples):
            for c_i in range(n_clusterings):
                v1 = sample_x_clustering_array[i, c_i]
                v2 = sample_x_clustering_array[j, c_i]
                if v1 and v2 and (v1 == v2):
                    coclusterings[i, j] += 1

    # Normalize by the number of clusterings and return
    coclusterings /= n_clusterings
    return coclusterings    

def _hierarchical_cluster_consensus_matrix(consensus_matrix, force_diagonal=True, method='ward'):
    """
    Hierarchical cluster consensus_matrix and compute cophenetic correlation coefficient.
    Convert consensus_matrix into distance matrix. Hierarchical cluster the distance matrix. And compute the
    cophenetic correlation coefficient.
    :param consensus_matrix: DataFrame;
    :param force_diagonal: bool;
    :param method: str; method parameter for scipy.cluster.hierarchy.linkage
    :return: ndarray float; linkage (Z) and cophenetic correlation coefficient
    """ 

    # Convert consensus matrix into distance matrix
    distance_matrix = 1 - consensus_matrix
    if force_diagonal:
        for i in range(distance_matrix.shape[0]):
            distance_matrix.iloc[i, i] = 0

    # Cluster consensus matrix to assign the final label
    hierarchical_clustering = linkage(consensus_matrix, method=method)

    # Compute cophenetic correlation coefficient
    cophenetic_correlation_coefficient = pearsonr(pdist(distance_matrix), cophenet(hierarchical_clustering))[0]

    return hierarchical_clustering, cophenetic_correlation_coefficient

        
matrix = pd.read_table("raw_data_file_EGFRmutant-mlike.out", sep='\s+')

random_seed = 12381
n_clusterings= 100
ks = 25


nmf_results = {}
cophenetic_correlation_coefficients = {}
annotation = pd.read_table("test_annotation_EGFRmutant-mlike.txt")

for k in range(2,ks):
    
    sample_x_clustering = pd.DataFrame(index=matrix.columns, columns=range(n_clusterings), dtype=int)
    print ("Running for:", k, "\n")
    for i in range(n_clusterings):
        nmf_result = nmf (matrix, k, random_seed +i)[k]
        #plotnmf(nmf_result['h'])     
        #print (nmf_result['h'])    
        sample_x_clustering.iloc[:, i] = argmax(asarray(nmf_result['h']), axis=0)
        #print (sample_x_clustering)
    
        if i == 0:
            nmf_results[k] = nmf_result
            #print (nmf_results[k]['w'])
            w = pd.DataFrame(nmf_results[k]['w'], index= matrix.index)
            outfile_w = "wmatrix" + "_"+ str(k)
            w.to_csv(outfile_w, sep='\t', encoding='utf-8')
            
            h = pd.DataFrame(nmf_results[k]['h'], columns= matrix.columns)
            outfile_h = "hmatrix" + "_"+ str(k)
            h.to_csv(outfile_h, sep='\t', encoding='utf-8')
            plot_nmf(h, annotation, outfile_h)

    #print (sample_x_clustering)
    consensus_matrix = _get_consensus(sample_x_clustering)
    consensus_matrix = pd.DataFrame(consensus_matrix)
    #print (consensus_matrix)

    hierarchical_clustering, cophenetic_correlation_coefficient = _hierarchical_cluster_consensus_matrix(consensus_matrix)
    #print (hierarchical_clustering, cophenetic_correlation_coefficient)
    cophenetic_correlation_coefficients[k] = cophenetic_correlation_coefficient

#print (nmf_results)
print (cophenetic_correlation_coefficients)
plot_coph(cophenetic_correlation_coefficients)