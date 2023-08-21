import pandas as pd
import numpy as np
import os

SAVEDIR='./data/processed/gene_trajectories'

if not os.path.isdir(SAVEDIR)
    os.mkdir(SAVEDIR)

# load temporal DEX genes from Li et al
all_temporal_dex = pd.read_excel('data/gene_data/gene_lists/mRNA-seq.Temporal.DEX.xlsx')

# keep only the thalamic data
md_temporal = all_temporal_dex.loc[:,all_temporal_dex.columns[all_temporal_dex.columns.str.contains('MD')]]

# genes with pre-vs-post enrichment
# defined as:
# a minimum of 3 pairwise sig. difference in gene expression between a prenatal and a postnatal window
md_temporal = md_temporal.loc[:,[
                                 'W1-vs-W6.MD',
                                 'W1-vs-W7.MD',
                                 'W1-vs-W8.MD',
                                 'W1-vs-W9.MD',
                                 'W2-vs-W6.MD',
                                 'W2-vs-W7.MD',
                                 'W2-vs-W8.MD',
                                 'W2-vs-W9.MD',
                                 'W3-vs-W6.MD',
                                 'W3-vs-W7.MD',
                                 'W3-vs-W8.MD',
                                 'W3-vs-W9.MD',
                                 'W4-vs-W6.MD',
                                 'W4-vs-W7.MD',
                                 'W4-vs-W8.MD',
                                 'W4-vs-W9.MD',
                                 'W5-vs-W6.MD',
                                 'W5-vs-W7.MD',
                                 'W5-vs-W8.MD',
                                 'W5-vs-W9.MD']
                             ]

# add gene symbols columns back on
md_temporal = pd.concat((all_temporal_dex.iloc[:,:2], md_temporal), axis=1)
# keep only protein-coding genes 
md_temporal = md_temporal[md_temporal.geneType == 'protein_coding']

# at least 3 significant pre-post comparisons
md_temporal_minp = np.isnan(md_temporal.iloc[:,2:]).values.sum(axis=1)<=17
# exclude any with no significant temporal DEX
md_temporal_dex = md_temporal[md_temporal_minp]

# get gene symbols of MD DEX genes
md_temporal_dex.insert(1, column='symbol', value=md_temporal_dex['Geneid'].str.split('|', expand=True).iloc[:,1])
md_temporal_dex_genes = list(md_temporal_dex['symbol'].unique())

# collect gene trajectories (calculated across all PsychENCODE samples)
trajectories = pd.read_csv('data/gene_data/gene_trajectories/data-gene-data-modelled-no-age.csv')
trajectories = trajectories[trajectories.region=='MD']

PCs=['PC1','PC2','PC3']

#OUTNAME='PC1'

for OUTNAME in PCs:

    # collect top 100 positive and top 100 negative
    positive = pd.read_csv('data/processed/HumanMostPositiveSpinTested_'+ OUTNAME + '.csv', header=None)[0].values
    negative = pd.read_csv('data/processed/HumanMostNegativeSpinTested_'+ OUTNAME + '.csv', header=None)[0].values
    # PC
    all_genes = pd.read_csv('data/gene_data/gene_lists/AllHumanSpinTested.csv', header=None)
    all_genes.columns = ['gene', 'pc','p']

    # how many are also temporal DEX genes?
    positive_dex = list(set(positive) & set(md_temporal_dex_genes))
    negative_dex = list(set(negative) & set(md_temporal_dex_genes))

    # get trajectories for each
    positive_trajectories =  trajectories[trajectories.symbol.isin(positive_dex)]
    negative_trajectories =  trajectories[trajectories.symbol.isin(negative_dex)]

    # add in PC coordinates
    pc_dict = dict(all_genes.iloc[:,:2].values)
    positive_trajectories.insert(1, 'PC', positive_trajectories.loc[:,'symbol'].map(pc_dict).values)
    negative_trajectories.insert(1, 'PC', negative_trajectories.loc[:,'symbol'].map(pc_dict).values)

    # add in windows
    windows = np.array(['prenatal']*len(positive_trajectories))
    windows[positive_trajectories.age.values>280] = 'postnatal'
    positive_trajectories.insert(1, 'window', windows)

    windows = np.array(['prenatal']*len(negative_trajectories))
    windows[negative_trajectories.age.values>280] = 'postnatal'
    negative_trajectories.insert(1, 'window', windows)

    positive_trajectories.to_csv(SAVEDIR+'/positive_trajectories_'+OUTNAME+'.csv')
    negative_trajectories.to_csv(SAVEDIR+'/negative_trajectories_'+OUTNAME+'.csv')

    print('out of {:} positive genes: {:} are also differentially expressed over time'.format(len(positive), len(positive_dex)))
    print('out of {:} negative genes: {:} are also differentially expressed over time'.format(len(negative), len(negative_dex)))

    # load mouse cell genes list
    # load Ensembl human-mouse genes
    all_homologs = pd.read_csv('data/gene_data/gene_lists/mart_export.txt', delimiter='\t')
    neuron_lists = pd.read_csv('data/gene_data/gene_lists/TH_Neuron_ALL.csv')['gene']
    glial_lists = pd.concat((pd.read_csv('data/gene_data/gene_lists/TH_Astrocyte_Gja1.csv'), 
                             pd.read_csv('data/gene_data/gene_lists/TH_Oligodendrocyte_ALL.csv'), 
                             pd.read_csv('data/gene_data/gene_lists/TH_Ependyma.csv')))['gene'].unique()

    # get human homologs
    neuron_homologs = all_homologs[all_homologs['Mouse gene name'].isin(list(neuron_lists))]['Gene name'].dropna()
    glial_homologs = all_homologs[all_homologs['Mouse gene name'].isin(list(glial_lists))]['Gene name'].dropna()

    # select trajectories
    positive_neuron_traj = positive_trajectories[positive_trajectories.symbol.isin(neuron_homologs.values)]
    negative_neuron_traj = negative_trajectories[negative_trajectories.symbol.isin(neuron_homologs.values)]
    positive_glial_traj = positive_trajectories[positive_trajectories.symbol.isin(glial_homologs.values)]
    negative_glial_traj = negative_trajectories[negative_trajectories.symbol.isin(glial_homologs.values)]

    # Gene lists for supplemental
    # DOWN over time and medially enriched
    a = negative_trajectories.groupby(by=['symbol', 'window'])[['fit']].mean().reset_index()
    b = a.groupby('symbol')['fit']
    # pre vs post - negative values is higher in prenatal than postnatal
    prepost = b.first().sub(b.last())
    down_neg = list(set(negative_trajectories.symbol.unique()) & set(prepost[prepost<0].index))

    # UP over time and medially enriched
    up_neg = list(set(negative_trajectories.symbol.unique()) & set(prepost[prepost>0].index))

    # DOWN over time and medially enriched
    a = positive_trajectories.groupby(by=['symbol', 'window'])[['fit']].mean().reset_index()
    b = a.groupby('symbol')['fit']
    # pre vs post - negative values is higher in prenatal than postnatal
    prepost = b.first().sub(b.last())
    down_pos = list(set(positive_trajectories.symbol.unique()) & set(prepost[prepost<0].index))

    # UP over time and medially enriched
    up_pos = list(set(positive_trajectories.symbol.unique()) & set(prepost[prepost>0].index))

    genelist=pd.DataFrame({'Negative_prenatal' : pd.Series(down_neg), 'Negative_postnatal' : pd.Series(up_neg), 'Positive_prenatal' : pd.Series(down_pos), 'Positive_postnatal' : pd.Series(up_pos)})
    genelist.to_csv(SAVEDIR+'/SigDevelopmentalTrajectoryGenes_' + OUTNAME +'.csv', index=False)