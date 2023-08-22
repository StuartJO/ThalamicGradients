# ThalamicGradients

This is code for "A phylogenetically-conserved axis of thalamocortical connectivity in the human brain", a preprint which can be found [here](https://www.biorxiv.org/content/10.1101/2022.11.15.516574v1)

Note that almost all of this code assumes relative paths based on running code from the main directory (i.e., the one this file is originally located in)

# Rerunning the code on the preprocessed data

The data I provide is at a minimum all preprocessed (basically it includes everything except the raw tractography data).

Note that to get the demographic data for the HCP subjects, you'll need to (apply for and then) download the restricted data from the HCP website

## Remaking all figures

In MATLAB run:




This will take around 20 minutes to run (with 75% of that just being the time it takes to save a generated plot).

The plots were then combined seperately in Inkscape

## Running WebGestalt

To run the data through [WebGestalt](https://www.webgestalt.org/), set it up with the following parameters:

Organism of Interest = Homo sapiens
Method of Interest = Over-Representation Analysis (ORA)

PCX_HumanMostPositiveSpinTested.txt
PCX_HumanMostNegativeSpinTested.txt

## Running from the start

If you realllllly want to rerun everything, I detail how to do so below (including how to get all the data from their original source). The one exception is I don't detail how to process the HCP data 

### Making seeds, and getting gene-expression and connectivity

The code doesn't include the raw diffusion or gene-expression data because that would use up far too much space.

The basic steps that were applied to the HCP data are as follows

5ttgen -tempdir ${WORKDIR} fsl ${T1} ${WORKDIR}/ACT.nii -premasked

If you did want to run things from scratch you need to have a parcellation registered to each individual (as a nifti volume) and for it to be aligned with the individuals diffusion space (which should be easy for HCP data).

First we need all of the gene-expression data so run

```
DownloadGeneData.sh
```

We also need a mask from the gene data. Using any of the 'X_mRNA.nii' files (I used gene number 12; but they should all give the exact same result) do the following:

fslmaths 12_mRNA.nii -bin GeneMask.nii.gz 

To generate the seeds run the following in bash:

MakeSeeds.sh

To obtain the final list of seeds to use, in MATLAB run

```
CheckDownloadedGenes()
SeedQC()
```

Once this is done you can then generate tractograms from each seed by running the following:
```
SUBJECT_LIST="/projects/kg98/stuarto/ThalamicGradients/VALIDSEED_UnrelatedSubs.txt"
nsubs=$(wc -l ${SUBJECT_LIST} | awk '{ print $1 }')
for ID in $(seq 1 $nsubs); do SUB=$(sed -n "${ID}p" ${SUBJECT_LIST}); sbatch ./MakeThalamicTracts.sh $SUB; done
```
## Making ancillary and other preprocessed data

I also provide ancillary data (which I define as data largely used to help plotting or assist other functions) already formatted, but if you wanted to create some of this from scratch, please run in MATLAB

MakeAncillaryData.m

This will make and format the mouse parcellation and flat map, get coordinates of mouse thalamic nuclei

### Getting mouse ancillary data

To run MakeAncillaryData.m you'll need to go and download the original flatmap [here](http://download.alleninstitute.org/publications/allen_mouse_brain_common_coordinate_framework/cortical_surface_views/ccf/annotation/flatmap_dorsal.nrrd)

Then to get a nifti file of the Allen Mouse Brain Atlas, go to the [Scalable Brain Atlas](https://scalablebrainatlas.incf.org/mouse/ABA_v3#downloads) and download the [P56_Atlas.nii.gz](https://scalablebrainatlas.incf.org/templates/ABA_v3/source/P56_Atlas.nii.gz) and [P56_Annotation.nii.gz](https://scalablebrainatlas.incf.org/templates/ABA_v3/source/P56_Annotation.nii.gz) files.


### Getting mouse cortical features

To get these, you need to download [this](https://github.com/benfulcher/mouseGradients) GitHub repository, download the required data (found [here](https://figshare.com/articles/dataset/Mouse_cortical_gradients/7775684/1)) and follow the instructions to run the code. To get the final output ('Data_AMBAcortex.mat') run the following in MATLAB

load('mouseGradients/DataOutputs/BigDamnMatrix.mat')
save('Data_AMBAcortex.mat','allProperties','dataMatrix','structInfo')

Then put 'Data_AMBAcortex.mat'  in ./data/preprocessed 

### Getting neuromaps data

In python, run

```
GetNeuroMaps.py
```
Note you need the [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) to be installed for the above to work.

I manually a seperate Excel file which gives the metadata about each brain map, using [this](https://docs.google.com/spreadsheets/d/1oZecOsvtQEh5pQkIf8cB6CyhPKVrQuko/edit#gid=1162991686) as a reference.

### Getting mouse connectivity and gene expression

Go to the following [link](https://doi.org/10.5281/zenodo.4609603) and download the "AllenGeneDataset_19419.mat" and "Mouse_Connectivity_Data.mat" files. Or just uses these links for the [gene](https://zenodo.org/record/4609603/files/AllenGeneDataset_19419.mat?download=1) and [connectivity](https://zenodo.org/record/4609603/files/Mouse_Connectivity_Data.mat?download=1) data.


### DropViz data

Use this [link](http://dropviz.org/?_state_id_=061c47a2583ad172) to go to the DropViz data and use all the parameters we used.

Under the "Target cluster" dropdown menu, load in a cluster. The table in the bottom of the webpage (called "Differentially Over-Expressed") will then update. There is a download button next to the header, use that to download the data. Repeat for all the clusters. This gives all the cell classes (we renamed each file according to the name of each cluster e.g.,"TH Fibroblast-Like_Dcn [#6]" would be named "TH_Fibroblast-Like_Dcn.csv". Note we made additional clusters by combining all neuron clusters into a single file, and combining all the oligodendrocyte clusters into a single file (with the "_ALL" suffix added.)

Next, up the top click the "Subclusters" header. Under the "Target Subcluster" dropdown menu, load in a subclusters. Repeat as above for only neuron subclusters (we renamed each file according to the _number_ of each cluster e.g., "TH Neuron.Gad1Gad2.Six3-Adcy1 [#3-1]" would be called "3-1.csv"). 

To get the TSNE coordinates, make sure you have selected the "Subclusters" tab and click the "tSNE" button underneath. Then click the "Display" tab in the Parameters box. Under "t-SNE Plot Settings", set the "Downsample Cell" to "Show all". Then click the download button up the top (or just click this link

To convert this to an Excel spreadsheet, unzip it and the run the following R commands:
```
library("writexl")
load("tsne_sub.Rdata")
write_xlsx(xy.data,"tsne_subdata.xlsx")
```

You'll also need "metacells.BrainCellAtlas_Saunders_version_2018.04.01.csv" and "annotation.BrainCellAtlas_Saunders_version_2018.04.01.csv" which you can find under the "Data" tab (found at the top of the page; Note we added "DropViz." as a prefix to the file name for the metacells file. No idea why we did this but we did).

### Gene trajectory data

Go to the [PsychENCODE website](http://development.psychencode.org/)

Under "Raw data", click "mRNA-seq" and download the [Sample Metadata](http://development.psychencode.org/files/raw_data/mRNA-seq_Sample%20metadata.xlsx) file

Under "Processed Data", click "mRNA-seq" and download the [Gene expression in RPKM](http://development.psychencode.org/files/processed_data/RNA-seq/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt), and [Regionally differential expressed genes](http://development.psychencode.org/files/processed_data/RNA-seq/mRNA-seq.Temporal.DEX.xlsx) files.
<!-- (Also under "Processed Data", click "Single cell/nucleus RNA-seq", and download all the Prenatal Single cell RNA-seq files (Sample QC; Gene expression CPM, counts and cell type cluster; Cell type signatures)) -->

### Human Protein Atlas

To get the list of genes with protein expression in the human thalamus you can go here https://v21.proteinatlas.org/humanproteome/brain/thalamus, click on the "14704" link and it will take you to the download page for all the genes.

You can also download the .tsv file directly from the following [link](https://v21.proteinatlas.org/api/search_download.php?search=NOT%20brain_category_rna%3Athalamus%3BNot%20detected&columns=g,gs,eg,gd,up,chr,chrp,pc,upbp,up_mf,di,pe,evih,eviu,evin,rnats,rnatd,rnatss,rnatsm,rnascs,rnascd,rnascss,rnascsm,rnacas,rnacad,rnacass,rnacasm,rnabrs,rnabrd,rnabrss,rnabrsm,rnabcs,rnabcd,rnabcss,rnabcsm,rnabls,rnabld,rnablss,rnablsm,rnacls,rnacld,rnaclss,rnaclsm,rnambrs,rnambrd,rnambrss,rnambrsm,rnapbrs,rnapbrd,rnapbrss,rnapbrsm,ab,relih,relmb,relce,scl,secl,ccdp,ccdt,ecblood,ectissue,eccellline,ecsinglecell,scml,scal,abrr,blconcia,blconcms,prognostic_breast_cancer,prognostic_cervical_cancer,prognostic_colorectal_cancer,prognostic_endometrial_cancer,prognostic_glioma,prognostic_head_and_neck_cancer,prognostic_liver_cancer,prognostic_lung_cancer,prognostic_melanoma,prognostic_ovarian_cancer,prognostic_pancreatic_cancer,prognostic_prostate_cancer,prognostic_renal_cancer,prognostic_stomach_cancer,prognostic_testis_cancer,prognostic_thyroid_cancer,prognostic_urothelial_cancer,t_RNA_adipose_tissue,t_RNA_adrenal_gland,t_RNA_amygdala,t_RNA_appendix,t_RNA_basal_ganglia,t_RNA_bone_marrow,t_RNA_breast,t_RNA_cerebellum,t_RNA_cerebral_cortex,t_RNA_cervix,t_RNA_choroid_plexus,t_RNA_colon,t_RNA_duodenum,t_RNA_endometrium_1,t_RNA_epididymis,t_RNA_esophagus,t_RNA_fallopian_tube,t_RNA_gallbladder,t_RNA_heart_muscle,t_RNA_hippocampal_formation,t_RNA_hypothalamus,t_RNA_kidney,t_RNA_liver,t_RNA_lung,t_RNA_lymph_node,t_RNA_medulla_oblongata,t_RNA_midbrain,t_RNA_ovary,t_RNA_pancreas,t_RNA_parathyroid_gland,t_RNA_pituitary_gland,t_RNA_placenta,t_RNA_pons,t_RNA_prostate,t_RNA_rectum,t_RNA_retina,t_RNA_salivary_gland,t_RNA_seminal_vesicle,t_RNA_skeletal_muscle,t_RNA_skin_1,t_RNA_small_intestine,t_RNA_smooth_muscle,t_RNA_spinal_cord,t_RNA_spleen,t_RNA_stomach_1,t_RNA_testis,t_RNA_thalamus,t_RNA_thymus,t_RNA_thyroid_gland,t_RNA_tongue,t_RNA_tonsil,t_RNA_urinary_bladder,t_RNA_vagina,t_RNA_white_matter,cell_RNA_A-431,cell_RNA_A549,cell_RNA_AF22,cell_RNA_AN3-CA,cell_RNA_ASC_diff,cell_RNA_ASC_TERT1,cell_RNA_BEWO,cell_RNA_BJ,cell_RNA_BJ_hTERT%2B,cell_RNA_BJ_hTERT%2B_SV40_Large_T%2B,cell_RNA_BJ_hTERT%2B_SV40_Large_T%2B_RasG12V,cell_RNA_CACO-2,cell_RNA_CAPAN-2,cell_RNA_Daudi,cell_RNA_EFO-21,cell_RNA_fHDF%2FTERT166,cell_RNA_GAMG,cell_RNA_HaCaT,cell_RNA_HAP1,cell_RNA_HBEC3-KT,cell_RNA_HBF_TERT88,cell_RNA_HDLM-2,cell_RNA_HEK_293,cell_RNA_HEL,cell_RNA_HeLa,cell_RNA_Hep_G2,cell_RNA_HHSteC,cell_RNA_HL-60,cell_RNA_HMC-1,cell_RNA_HSkMC,cell_RNA_hTCEpi,cell_RNA_hTEC%2FSVTERT24-B,cell_RNA_hTERT-HME1,cell_RNA_hTERT-RPE1,cell_RNA_HUVEC_TERT2,cell_RNA_JURKAT,cell_RNA_K-562,cell_RNA_Karpas-707,cell_RNA_LHCN-M2,cell_RNA_MCF7,cell_RNA_MOLT-4,cell_RNA_NB-4,cell_RNA_NTERA-2,cell_RNA_OE19,cell_RNA_PC-3,cell_RNA_REH,cell_RNA_RH-30,cell_RNA_RPMI-8226,cell_RNA_RPTEC_TERT1,cell_RNA_RT4,cell_RNA_SCLC-21H,cell_RNA_SH-SY5Y,cell_RNA_SiHa,cell_RNA_SK-BR-3,cell_RNA_SK-MEL-30,cell_RNA_SuSa,cell_RNA_T-47d,cell_RNA_THP-1,cell_RNA_TIME,cell_RNA_U-138_MG,cell_RNA_U-2_OS,cell_RNA_U-2197,cell_RNA_U-251_MG,cell_RNA_U-266%2F70,cell_RNA_U-266%2F84,cell_RNA_U-698,cell_RNA_U-87_MG,cell_RNA_U-937,cell_RNA_WM-115,blood_RNA_basophil,blood_RNA_classical_monocyte,blood_RNA_eosinophil,blood_RNA_gdT-cell,blood_RNA_intermediate_monocyte,blood_RNA_MAIT_T-cell,blood_RNA_memory_B-cell,blood_RNA_memory_CD4_T-cell,blood_RNA_memory_CD8_T-cell,blood_RNA_myeloid_DC,blood_RNA_naive_B-cell,blood_RNA_naive_CD4_T-cell,blood_RNA_naive_CD8_T-cell,blood_RNA_neutrophil,blood_RNA_NK-cell,blood_RNA_non-classical_monocyte,blood_RNA_plasmacytoid_DC,blood_RNA_T-reg,blood_RNA_total_PBMC,brain_RNA_amygdala,brain_RNA_basal_ganglia,brain_RNA_cerebellum,brain_RNA_cerebral_cortex,brain_RNA_hippocampal_formation,brain_RNA_hypothalamus,brain_RNA_medulla_oblongata,brain_RNA_midbrain,brain_RNA_pons,brain_RNA_spinal_cord,brain_RNA_thalamus,brain_RNA_white_matter,sc_RNA_Adipocytes,sc_RNA_Alveolar_cells_type_1,sc_RNA_Alveolar_cells_type_2,sc_RNA_Astrocytes,sc_RNA_B-cells,sc_RNA_Basal_keratinocytes,sc_RNA_Basal_prostatic_cells,sc_RNA_Basal_respiratory_cells,sc_RNA_Basal_squamous_epithelial_cells,sc_RNA_Bipolar_cells,sc_RNA_Breast_glandular_cells,sc_RNA_Breast_myoepithelial_cells,sc_RNA_Cardiomyocytes,sc_RNA_Cholangiocytes,sc_RNA_Club_cells,sc_RNA_Collecting_duct_cells,sc_RNA_Cone_photoreceptor_cells,sc_RNA_Cytotrophoblasts,sc_RNA_dendritic_cells,sc_RNA_Distal_enterocytes,sc_RNA_Distal_tubular_cells,sc_RNA_Ductal_cells,sc_RNA_Early_spermatids,sc_RNA_Endometrial_ciliated_cells,sc_RNA_Endometrial_stromal_cells,sc_RNA_Endothelial_cells,sc_RNA_Enteroendocrine_cells,sc_RNA_Erythroid_cells,sc_RNA_Excitatory_neurons,sc_RNA_Exocrine_glandular_cells,sc_RNA_Extravillous_trophoblasts,sc_RNA_Fibroblasts,sc_RNA_Gastric_mucus-secreting_cells,sc_RNA_Glandular_and_luminal_cells,sc_RNA_granulocytes,sc_RNA_Granulosa_cells,sc_RNA_Hepatic_stellate_cells,sc_RNA_Hepatocytes,sc_RNA_Hofbauer_cells,sc_RNA_Horizontal_cells,sc_RNA_Inhibitory_neurons,sc_RNA_Intestinal_goblet_cells,sc_RNA_Ionocytes,sc_RNA_Kupffer_cells,sc_RNA_Langerhans_cells,sc_RNA_Late_spermatids,sc_RNA_Leydig_cells,sc_RNA_Macrophages,sc_RNA_Melanocytes,sc_RNA_Microglial_cells,sc_RNA_monocytes,sc_RNA_Muller_glia_cells,sc_RNA_NK-cells,sc_RNA_Oligodendrocyte_precursor_cells,sc_RNA_Oligodendrocytes,sc_RNA_Pancreatic_endocrine_cells,sc_RNA_Paneth_cells,sc_RNA_Peritubular_cells,sc_RNA_Plasma_cells,sc_RNA_Prostatic_glandular_cells,sc_RNA_Proximal_enterocytes,sc_RNA_Proximal_tubular_cells,sc_RNA_Respiratory_ciliated_cells,sc_RNA_Rod_photoreceptor_cells,sc_RNA_Sertoli_cells,sc_RNA_Skeletal_myocytes,sc_RNA_Smooth_muscle_cells,sc_RNA_Spermatocytes,sc_RNA_Spermatogonia,sc_RNA_Squamous_epithelial_cells,sc_RNA_Suprabasal_keratinocytes,sc_RNA_Syncytiotrophoblasts,sc_RNA_T-cells,sc_RNA_Theca_cells,sc_RNA_Undifferentiated_cells,sc_RNA_Urothelial_cells&compress=no&format=tsv)

Note that we used v21.1 of the Human Protein Atlas. Now the atlas is up to version 23.0 (as of 22 August 2023) and has found an extra ~2000 genes expressed in the thalamus! As an exercise for the reader they could download this new list and see how things replicate!

## Homologue mouse-human genbes from Ensembl Biomart

If you click this [link](http://asia.ensembl.org/biomart/martview/4a48fd04cdba815796782ab3e8bc620b?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.homologs.ensembl_gene_id|hsapiens_gene_ensembl.default.homologs.mmusculus_homolog_ensembl_gene|hsapiens_gene_ensembl.default.homologs.mmusculus_homolog_associated_gene_name|hsapiens_gene_ensembl.default.homologs.external_gene_name&FILTERS=&VISIBLEPANEL=attributepanel), it will set up all the parameters (in the correct order as well) for you to download (download by setting "Export all results to" to "File" and "CSV", then click go and it will download a txt file) 

## FreeSurfer and FSL data

Just have FreeSurfer and FSL installed. The code should pick up things for FSL. For FreeSurfer you'll need to dig through and extract the surface files (e.g., inflated, white, pial, sphere .etc) from whereever the fsaverage 164k surface is.

## Melbourne subcortical atlas