import nibabel as nb
import numpy as np
from fragmenter import Fragment
from fragmenter import adjacency
from fragmenter import RegionExtractor

vertices,faces = nb.freesurfer.io.read_geometry('./data/parecllation/lh.sphere')

# The medial wall in lh.aparc.annot doesn't for a single continuous region 
# on the surface. We made a version which fixed this (by just assigning 
# vertices of the medial wall that were not part of the larger group to the
# ROI they were closest to 
E = RegionExtractor.Extractor('./data/parecllation/lh.aparc_connected.annot')
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

leftcortex = Fragment.Fragment(n_clusters=250)
leftcortex.fit(vertices=vertices, faces=faces,
  parcels=parcels, rois=rois, method = 'k_means') 
  
annot_name = './data/parecllation/lh.random500.annot'
leftcortex.write(annot_name)   

vertices,faces = nb.freesurfer.io.read_geometry('./data/parecllation/rh.sphere')

E = RegionExtractor.Extractor('./data/parecllation/rh.aparc_connected.annot')
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

rightcortex = Fragment.Fragment(n_clusters=250)
rightcortex.fit(vertices=vertices, faces=faces,
  parcels=parcels, rois=rois, method = 'k_means') 
  
annot_name = './data/parecllation/rh.random500.annot'
rightcortex.write(annot_name)   