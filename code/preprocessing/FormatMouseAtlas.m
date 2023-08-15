function FormatMouseAtlas()

% Lets format the mouse data into something workable!

load ('./data/preprocessed/AllenGeneDataset_19419.mat','structInfo')

CCFv3_annots = Get_CCFv3_annots;

% First up lets find the CCFv3 ID for the regions in structInfo (i.e.,
% those which correspond to the gene and tract data

[~,idx] = ismember(structInfo.name,CCFv3_annots.name);

% This gives what the ID of the Oh parcellation is in the CCFv3 atlas
mouse_CCFv3_id = CCFv3_annots.ID(idx);

gunzip('./data/ancillary/P56_Annotation.nii.gz','./data/ancillary')

MouseAtlas = double(niftiread('./data/ancillary/P56_Annotation.nii'));

% Lets keep storage down
delete('./data/ancillary/P56_Annotation.nii')
  
% Because the voxel IDs in MouseAtlas are labelled according to the most
% fine grained areas in a region, we need to find all the regions which
% belong to a "meta" region and group them together. Fortunately these
% subregions are (usually) labelled with a prefix that indicates which meta
% region it belongs do.

MouseOhParc = zeros(size(MouseAtlas));

% Number areas according to the Oh ID
for i = 1:length(mouse_CCFv3_id)   
    formatName = structInfo.name{i};
    ids = CCFv3_annots.ID(contains(CCFv3_annots.name,formatName));
    MouseOhParc(ismember(MouseAtlas,ids)) = i;
end

%% Do some manual fixes

% Using info from Table S2 from https://doi.org/10.1016/j.cell.2020.04.007,
% I did some manual corrections

% Subregions of the "Frontal pole" don't share the exact same prefix as the
% meta region

MouseOhParc(ismember(MouseAtlas,7:10))=8;

% The "Posterior parietal association areas" is made up of multiple regions
% which don't share its main name. The 'Anterior area' and 'Rostrolateral
% area' make it up 

MouseOhParc(ismember(MouseAtlas,346:359))=1;

% Fix "Superior colliculus, sensory related" areas
MouseOhParc(ismember(MouseAtlas,809:811))=148;

save('./data/ancillary/MouseOhParc.mat','MouseOhParc')

% The following regions are in the parcellation because we cannot (or more 
% accurately, I could not figure out how to) find a mapping between their 
% structure ID and any of the IDs in the CCFv3 data. But none of these 
% regions are cortex or thalamus so we don't care about them :)

%NotInParc = find(~ismember(1:213,unique(MouseOhParc)));
%structInfo.name(NotInParc)

% {'Subiculum, dorsal part' }
% {'Subiculum, ventral part'}
% {'Central lobule'         }
% {'Culmen'                 }
% {'Ansiform lobule'        }

%% Get the nifti data from the compressed file

gunzip('./data/ancillary/P56_Atlas.nii.gz','./data/ancillary')

MouseBrain = double(niftiread('./data/ancillary/P56_Atlas.nii'));
save('./data/ancillary/MouseBrain.mat','MouseBrain')

% The uncompressed data is deleted to save space
delete('./data/ancillary/P56_Atlas.nii')