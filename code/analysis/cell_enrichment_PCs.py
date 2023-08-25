import pandas as pd
import numpy as np
import scipy as sc
from scipy.stats import hypergeom

import os, glob

#import matplotlib.pyplot as plt
#import seaborn as sns

print(sc.__version__)
print(pd.__version__)
print(np.__version__)

SAVEDIR='./data/processed/cell_enrichment'

if not os.path.exists(SAVEDIR):
    os.mkdir(SAVEDIR)

# Some functions for gene enrichment analyses
def calculate_enrichment(hit_list, top_genes, full_gene_list):
    x = sum(pd.DataFrame(top_genes).isin(hit_list).values) # how many top genes in cell list
    n = sum(pd.DataFrame(hit_list).isin(full_gene_list).values)[0] # how many cell genes in full list
    N = len(top_genes)  # number of samples
    M = len(full_gene_list)  # total number in population

    enrichment = safe_div( (x/N) , ((n-x) / (M-N)) )
    p = hypergeom.sf(x-1, M, n, N)

    return enrichment, p

def safe_div(x,y):
    if y == 0:
        return np.array([0])
    return x / y

def run_enrichment(classes, gene_lists, positive_genes, negative_genes, background_genes):

    enrichment_results = []
    num_genes = []
    # for each cell class/type
    for i in np.arange(len(classes)):
        # calculate enrichment in the postively and negatively correlated lists
        for g in [positive_genes, negative_genes]:
            enrichment_results.append(calculate_enrichment(list(gene_lists[i]), list(g), list(background_genes)))
            num_genes.append(len(gene_lists[i]))

    # collate into dataframe
    results = pd.DataFrame(np.hstack(enrichment_results).T)
    results.columns=['enrichment', 'p']
    results['class'] = np.repeat(classes, 2)
    results['loading'] = ['positive', 'negative']*(len(classes))
    results['num_genes'] = np.squeeze(num_genes)
    results = results.loc[:,['class','loading','num_genes','enrichment','p']]

    return results

def run_enrichment_percentile(classes, gene_lists, test_genes, background_genes):

    enrichment_results = []
    num_genes = []
    # for each cell class/type
    for i in np.arange(len(classes)):
        # calculate enrichment in the postively and negatively correlated lists
        for g in test_genes:
            enrichment_results.append(calculate_enrichment(list(gene_lists[i]), list(g), list(background_genes)))
            num_genes.append(len(gene_lists[i]))

    # collate into dataframe
    results = pd.DataFrame(np.hstack(enrichment_results).T)
    results.columns=['enrichment', 'p']
    results['class'] = np.repeat(classes, 10)
    results['percentile'] = ['10', '20', '30', '40', '50', '60', '70', '80', '90', '100']*(len(classes))
    results['num_genes'] = np.squeeze(num_genes)
    results = results.loc[:,['class','percentile','num_genes','enrichment','p']]

    return results

# load Ensembl human-mouse genes
all_homologs = pd.read_csv('./data/gene_data/gene_lists/mart_export.txt', delimiter='\t')

# load genes expressed in human thalamus and filter homologs
thalamus_proteins = pd.read_csv('./data/gene_data/gene_lists/human_thalamus_protein.tsv', delimiter='\t', low_memory=False)
thalamus_proteins = thalamus_proteins['Gene']

all_homologs = all_homologs[all_homologs['Gene name'].isin(list(thalamus_proteins))]

# mouse genes with human homologs expressed in thalamus
mouse_homologs = all_homologs['Mouse gene name'].dropna()

# load genes assayed in Saunders et al (DropViz)
# this is a list of 'metacells' where each gene has a total count of transcripts for each cell cluster
# as many genes are assayed, but not all detected in the thalamus, we remove those with low UMI counts (<25%ile)
thalamic_cells = pd.read_csv('data/gene_data/gene_lists/DropViz.metacells.BrainCellAtlas_Saunders_version_2018.04.01.csv')
aggregated_UMI = thalamic_cells.loc[:,thalamic_cells.columns.str.contains('TH')].sum(axis=1)
thalamic_cells = thalamic_cells.loc[aggregated_UMI > np.percentile(aggregated_UMI, 25)]
# get list of genes
thalamic_cell_genes = thalamic_cells.iloc[:,0]

# filter to include only those with homologs expressed in human thalamus
thalamic_cell_genes = list(thalamic_cell_genes[thalamic_cell_genes.isin(mouse_homologs)])
print('The total number of background genes (measured in mouse thalamus using Drop-seq with human homologs also expressed in human thalamus):\n{:}'.format(len(thalamic_cell_genes)))

# get cell type data
cell_data = pd.read_csv('./data/gene_data/gene_lists/thalamus_cell_cluster_types.csv', delimiter='\t')
cell_data

# load differentially expressed gene lists from DropViz for each cell type
gene_lists = [list(pd.read_csv('data/gene_data/gene_lists/TH_' + cluster + '.csv')['gene']) for cluster in cell_data['Cluster']]

# filter out genes not included in the background set
filtered_gene_lists = [list(set(cluster_list) & set(thalamic_cell_genes)) for cluster_list in gene_lists]

# add to dataframe
cell_data['genes'] = filtered_gene_lists

#drop neuronal and oligodendrocyte subclusters for now
celltype_data = cell_data[~cell_data['Cluster'].isin(['Neuron_Habenula_Tac2', 'Neuron_Rora', 'Neuron_Gad2-Ahi1', 'Oligodendrocyte_Mbp', 'Oligodendrocyte_Tfr' ])]

# concatenate gene lists for each cell classes
cell_class = list(celltype_data.Class.unique())
cell_class_lists = [list(np.unique(np.concatenate(celltype_data[celltype_data['Class'] == cell]['genes'].values))) for cell in cell_class]

cell_class = [cell_class_.replace('Microglia_Macrophage','Microglia/Macrophage') for cell_class_ in cell_class]
print(cell_class)

print('Total number of DE genes per class after filtering\n')
for n,i in enumerate(cell_class):
    print('{:}: {:} genes'.format(i, len(cell_class_lists[n])))

# gene lists for each neuronal subcluster
cluster_lists =  sorted(glob.glob('./data/gene_data/gene_lists/subclusters/*csv'))
#print(cluster_lists)
# name for each cluster
cluster_names = [('Cluster '+ s.split('\\')[-1].split('.')[0]) for s in cluster_lists]

cluster_names = [cluster_name.replace('Cluster 1','Habenula').replace('Cluster 2','Rora').replace('Cluster 3','Gad2/Ahi1') for cluster_name in cluster_names]
print(cluster_names)
# get cluster genes
cluster_genes = [list(pd.read_csv(s)['gene']) for s in cluster_lists]

# filter out genes not included in thalamic background
filtered_cluster_genes = [list(set(s) & set(thalamic_cell_genes)) for s in cluster_genes]

# all genes expressed by neuron class
neuron_genes = pd.read_csv('./data/gene_data/gene_lists/all-neuronal-genes.csv')['gene']
# filter out any not in thalamic background
neuron_genes = neuron_genes[neuron_genes.isin(thalamic_cell_genes)]

PCs=['PC1','PC2','PC3']

#OUTNAME='PC3'

for OUTNAME in PCs:

	# PC+ and PC- genes - change for the updated list of 100
	# load top 100
	positive_genes = pd.read_csv('./data/processed/'+OUTNAME+'_HumanMostPositiveSpinTested.csv', header=None)[0].values
	# identify homologs
	positive_homologs = all_homologs[all_homologs['Gene name'].isin(list(positive_genes))]['Mouse gene name'].dropna()
	# remove any not in background
	filtered_positive_homologs = list(set(positive_homologs) & set(thalamic_cell_genes))

	print('From top {:} positive genes'.format(len(positive_genes)))
	print('{:} corresponding mouse genes'.format(len(positive_homologs)))
	print('with {:} in background set'.format(len(filtered_positive_homologs)))
	print('')

	# as above but for bottom 100
	negative_genes = pd.read_csv('./data/processed/'+OUTNAME+'_HumanMostNegativeSpinTested.csv', header=None)[0].values
	negative_homologs = all_homologs[all_homologs['Gene name'].isin(list(negative_genes))]['Mouse gene name'].dropna()
	filtered_negative_homologs = list(set(negative_homologs) & set(thalamic_cell_genes))

	print('From top {:} negative genes'.format(len(negative_genes)))
	print('{:} corresponding mouse genes'.format(len(negative_homologs)))
	print('with {:} in background set'.format(len(filtered_negative_homologs)))

	enrichment = run_enrichment(cell_class, cell_class_lists, filtered_positive_homologs, filtered_negative_homologs, thalamic_cell_genes)
	enrichment['p<0.01'] = enrichment['p']<0.01
	enrichment['p<0.001'] = enrichment['p']<0.001
	enrichment['p<0.0001'] = enrichment['p']<0.0001

	#enrichment

	enrichment.to_csv(SAVEDIR+'cell_enrichment_' + OUTNAME +'.csv')

	# repeat the enrichment analysis focusing on neuronal subtypes
	neuron_cell_data = cell_data[cell_data['Class'] == 'Neuron']
	# drop whole neuron group
	neuron_cell_data = neuron_cell_data.loc[~(neuron_cell_data['Cluster']=='Neuron_ALL')]
	#neuron_cell_data

	# subtype enrichment

	#neuron_clusters_orig = list(neuron_cell_data['Cluster'])

	#neuron_cluster_names = ['Cluster 1', 'Cluster 2', 'Cluster 3']
	neuron_cluster_names = ['Habenula', 'Rora', 'Gad2/Ahi1']

	neuron_enrichment = run_enrichment(neuron_cluster_names, list(neuron_cell_data['genes']), filtered_positive_homologs, filtered_negative_homologs, thalamic_cell_genes)
	neuron_enrichment['p<0.01'] = neuron_enrichment['p']<0.01
	neuron_enrichment['p<0.001'] = neuron_enrichment['p']<0.001
	neuron_enrichment['p<0.0001'] = neuron_enrichment['p']<0.0001
	#neuron_enrichment

	neuron_enrichment.to_csv(SAVEDIR+'/neuron_enrichment_' + OUTNAME +'.csv')

	for n, i in enumerate(list(neuron_cell_data['Cluster'])):
		print('')
		print('genes shared between postive top 100 and {:}'.format(i))
		print(list(set(filtered_positive_homologs) & set(list(neuron_cell_data['genes'])[n])))
		save_name=format(i)
		genelist=pd.DataFrame({'Gene' : list(set(filtered_positive_homologs) & set(list(neuron_cell_data['genes'])[n]))})
		genelist.to_csv(SAVEDIR+'/Positive'+save_name+'Genes_' + OUTNAME +'.csv')

	for n, i in enumerate(list(neuron_cell_data['Cluster'])):
		print('')
		print('genes shared between negative bottom 100 and {:}'.format(i))
		print(list(set(filtered_negative_homologs) & set(list(neuron_cell_data['genes'])[n])))
		save_name=format(i)
		genelist=pd.DataFrame({'Gene' : list(set(filtered_negative_homologs) & set(list(neuron_cell_data['genes'])[n]))})
		genelist.to_csv(SAVEDIR+'/Negative'+save_name+'Genes_' + OUTNAME +'.csv')

	subcluster_enrichment = run_enrichment(cluster_names, cluster_genes, filtered_positive_homologs, filtered_negative_homologs, neuron_genes)
	subcluster_enrichment.rename(columns={"class": "subcluster"}, inplace=True)
	subcluster_enrichment

	subcluster_enrichment.to_csv(SAVEDIR+'/subcluster_enrichment_' + OUTNAME +'.csv')

	# split into positive and negative lists
	positive_subcluster_enrichment = subcluster_enrichment.loc[subcluster_enrichment.loading=='positive'].copy().reset_index(drop=True)
	negative_subcluster_enrichment = subcluster_enrichment.loc[subcluster_enrichment.loading=='negative'].copy().reset_index(drop=True)
	#print(negative_subcluster_enrichment)
