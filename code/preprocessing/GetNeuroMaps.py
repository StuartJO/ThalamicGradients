from neuromaps import datasets
from neuromaps import transforms
import nibabel as nb

# This gets all the current maps in neuromaps and puts them into fsaverage164k space

#datasets.fetch_all_atlases()
OUTDIR='./data/preprocessed/fsaverage164k_annots'
#ANNOTS = fetch_annotation(source='all')
ANNOTS=datasets.available_annotations()
for annots in ANNOTS:
    print(annots)
    if annots[2] == 'MNI152':
        annot_dir=datasets.fetch_annotation(source=annots[0],desc=annots[1])
        fsdata = transforms.mni152_to_fsaverage(annot_dir, '164k')
        LH, RH = fsdata
    elif annots[2] == 'fsLR':
        annot_dir=datasets.fetch_annotation(source=annots[0],desc=annots[1])
        fsdata = transforms.fslr_to_fsaverage(annot_dir, '164k',hemi='L')
        LH = fsdata[0]
    elif annots[2] == 'civet':
        annot_dir=datasets.fetch_annotation(source=annots[0],desc=annots[1])
        fsdata = transforms.civet_to_fsaverage(annot_dir, '164k',hemi='L')
        LH = fsdata[0]
    elif annots[2] == 'fsaverage':
        annot_dir=datasets.fetch_annotation(source=annots[0],desc=annots[1])
        fsdata = transforms.fsaverage_to_fsaverage(annot_dir, '164k',hemi='L')
        LH = fsdata[0]
   
       
    nb.save(LH, OUTDIR + '/'+annots[0]+ '_' +annots[1]+'_fsaverage164k_hemi-L.func.gii')