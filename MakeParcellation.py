import nibabel as nb
import numpy as np
from fragmenter import Fragment
from fragmenter import adjacency
from fragmenter import RegionExtractor

vertices,faces = nb.freesurfer.io.read_geometry('F:/Documents/MATLAB/GenModelCode/ROIcontinuity/lh.sphere')

E = RegionExtractor.Extractor('F:/Documents/MATLAB/GenModelCode/ROIcontinuity/lh.aparc.annot')
parcels = E.map_regions()

rois=['bankssts','caudalanteriorcingulate','caudalmiddlefrontal','cuneus',
	'entorhinal','fusiform','inferiorparietal','inferiortemporal',
	'isthmuscingulate','lateraloccipital','lateralorbitofrontal',
	'lingual','medialorbitofrontal','middletemporal','parahippocampal',
	'paracentral','parsopercularis','parsorbitalis','parstriangularis',
	'pericalcarine','postcentral','posteriorcingulate','precentral',
	'precuneus','rostralanteriorcingulate','rostralmiddlefrontal',
	'superiorfrontal','superiorparietal','superiortemporal','supramarginal',
	'frontalpole','temporalpole','transversetemporal','insula']

leftcortex = Fragment.Fragment(n_clusters=100)
leftcortex.fit(vertices=vertices, faces=faces,
  parcels=parcels, rois=rois, method = 'k_means') 
  
annot_name = 'F:/Documents/MATLAB/GenModelCode/ROIcontinuity/lh.wholebrain100.annot'
leftcortex.write(annot_name)   

leftcortex = Fragment.Fragment(n_clusters=250)
leftcortex.fit(vertices=vertices, faces=faces,
  parcels=parcels, rois=rois, method = 'k_means') 
  
annot_name = 'F:/Documents/MATLAB/GenModelCode/ROIcontinuity/lh.wholebrain250.annot'
leftcortex.write(annot_name)   

vertices,faces = nb.freesurfer.io.read_geometry('F:/Documents/MATLAB/GenModelCode/ROIcontinuity/rh.sphere')

E = RegionExtractor.Extractor('F:/Documents/MATLAB/GenModelCode/ROIcontinuity/rh.aparc.annot')
parcels = E.map_regions()

rois=['bankssts','caudalanteriorcingulate','caudalmiddlefrontal','cuneus',
	'entorhinal','fusiform','inferiorparietal','inferiortemporal',
	'isthmuscingulate','lateraloccipital','lateralorbitofrontal',
	'lingual','medialorbitofrontal','middletemporal','parahippocampal',
	'paracentral','parsopercularis','parsorbitalis','parstriangularis',
	'pericalcarine','postcentral','posteriorcingulate','precentral',
	'precuneus','rostralanteriorcingulate','rostralmiddlefrontal',
	'superiorfrontal','superiorparietal','superiortemporal','supramarginal',
	'frontalpole','temporalpole','transversetemporal','insula']

rightcortex = Fragment.Fragment(n_clusters=100)
rightcortex.fit(vertices=vertices, faces=faces,
  parcels=parcels, rois=rois, method = 'k_means') 
  
annot_name = 'F:/Documents/MATLAB/GenModelCode/ROIcontinuity/rh.wholebrain100.annot'
rightcortex.write(annot_name)   

rightcortex = Fragment.Fragment(n_clusters=250)
rightcortex.fit(vertices=vertices, faces=faces,
  parcels=parcels, rois=rois, method = 'k_means') 
  
annot_name = 'F:/Documents/MATLAB/GenModelCode/ROIcontinuity/rh.wholebrain250.annot'
rightcortex.write(annot_name)   