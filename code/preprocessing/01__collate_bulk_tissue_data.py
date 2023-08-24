#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os

def run():
    ##################################################################
    # Load PsychENCODE data and save out expression data and metadata
    ##################################################################
    print("")
    print("preprocessing PsychENCODE bulk tissue data..")

    data_dir = './data/gene_data/'
    rna_data_dir = data_dir + 'gene_trajectories/'


    # LOAD
    # load metadata
    metadata = pd.read_excel(rna_data_dir + 'PsychENCODE/mRNA-seq_Sample metadata.xlsx', skiprows=[0,1,2])
    metadata.drop(labels=0, axis='rows', inplace=True)
    metadata.index=np.arange(len(metadata))

    # load expression data
    rpkm = pd.read_csv(rna_data_dir + 'PsychENCODE/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt', delimiter='\t')
    ## adjust column titles
    rpkm.insert(1, column='symbol', value=rpkm['Geneid'].str.split('|', expand=True).iloc[:,1])
    rpkm.insert(1, column='ensembl', value=rpkm['Geneid'].str.split('|', expand=True).iloc[:,0])
    rpkm.drop(labels='Geneid', axis='columns', inplace=True)


    # load list of protein-coding genes
    protein_coding_genes = pd.read_csv(data_dir + 'gene_lists/human_thalamus_protein.tsv', delimiter='\t', low_memory=False)
    protein_coding_gene_list = np.unique(protein_coding_genes['Gene'])


    # PROCESS
    # keep all samples for trajectory calculations
    all_metadata = metadata.loc[(metadata['Days']>0),:]

    # filter out non-coding genes
    rpkm = rpkm.loc[rpkm['symbol'].isin(protein_coding_gene_list),:]

    # get RPKM for each sample of interest
    long_rpkm = rpkm.melt(id_vars=['symbol'],
                          value_vars=rpkm.columns[pd.Series(rpkm.columns).str.contains('HSB')],
                          var_name='sampleid',
                          value_name='rpkm')
    long_rpkm = long_rpkm.groupby(by=['symbol', 'sampleid']).mean().reset_index().copy()

    long_rpkm['region'] = long_rpkm['sampleid'].str.split('.', expand=True).iloc[:,1]
    long_rpkm['sample'] = long_rpkm['sampleid'].str.split('.', expand=True).iloc[:,0]
    long_rpkm = long_rpkm.loc[:,['symbol', 'sample', 'region', 'rpkm']]

    all_data = long_rpkm.loc[long_rpkm.loc[:,'sample'].isin(all_metadata['Braincode'].values),:].copy()

    # add in metadata variables
    for s in np.unique(metadata['Braincode']):
        for n, (new, orig) in enumerate(zip(['age','sex','eth','RIN','seq_site','PMI'],
                                            ['Days', 'Sex', 'Ethnicity', 'RIN', 'Sequencing Site', 'PMI'])):
            all_data.loc[all_data.loc[:,'sample']==s,new] = metadata.loc[metadata.loc[:,'Braincode']==s, orig].values[0]

    # dorsal thalamic anlage in embryonic samples -> medial dorsal to match others
    region_dict = {'DTH':'MD'}
    all_data.loc[:,'region'] = all_data['region'].map(region_dict).fillna(all_data['region']).values

    # save out
    os.makedirs(rna_data_dir+'generated_data/', exist_ok=True)
    all_data.to_csv(rna_data_dir+'generated_data/PsychENCODE-ALL-bulk-RPKM-data.csv', index=False)
    all_metadata.to_csv(rna_data_dir+'generated_data/PsychENCODE-ALL-bulk-metadata.csv', index=None)
    print("see generated_data/PsychENCODE-ALL-bulk-RPKM-data.csv")
    print("see generated_data/PsychENCODE-ALL-bulk-metadata.csv")

if __name__ == '__main__':
    run()
