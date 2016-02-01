### IMPORTS ###

#NeuroImage
import nibabel as nib
#Matrix operations
import numpy as np
#Others
import math
import argparse

#returns tissue as int
def tissue2int(tissue_type):
	if tissue_type == 'BG':
		return 0
	elif tissue_type == 'CSF':
		return 1
	elif tissue_type == 'GM':
		return 2
	else:
		return 3
#computes n samples on a unit sphere located at the origin
def random_sampling(n):
	samples = np.zeros((n,3))
	np.random.seed(n)
	zs = np.random.uniform(-1.0,1.0,n)
	np.random.seed(2*n)
	ts = np.random.uniform(0,2*math.pi,n)
	rs = np.sqrt(1-(zs**2))
	xs = rs * np.cos(ts)
	ys = rs * np.sin(ts)
	samples[:,0] = xs
	samples[:,1] = ys
	samples[:,2] = zs

	return samples
	 		
def compute_inclusion_mask(roi,dim=1):
	""" Computes an inclusion mask finding the center of the region and then setting 
	    each voxel in that slice outside the ROI to 0 and all other voxels to 1. 
	    
	    Parameters
 		---------
	    roi_img: numpy 3D-array (mandatory)
	    	Atlas image from which to extract the mask.
		dim: int (Default: 1)
	    	Dimension along which to find the middle slice.
		
		Returns
		-----------
		inclusion_mask: numpy 3D-array
			Inclusion mask.
	"""
	#find the middle slice
	indices = np.where(roi > 0)[dim]
	min_dim, max_dim = np.min(indices), np.max(indices)
	center = (min_dim + max_dim) / 2
	#create inclusion mask
	for x in range(roi.shape[0]):
		for y in range(roi.shape[1]):
			for z in range(roi.shape[2]):
				slice = (((dim==0 and x==center) or (dim==1 and y==center)) or (dim==2 and z==center))
				#if voxel is in middle slice
				if slice: 
				 if roi[x,y,z] > 0: roi[x,y,z] = 1
				#outside the middle slice every voxel is set to 1
				else: roi[x,y,z] = 1
	return roi
	
def normalize_tract(tract):
	""" Normalizes a tract image so that its values sum up to 1.
	    
	    Parameters
 		---------
	    tract: nibabel Nifti1Image (mandatory)
	    	Tract image to normalize.
		
		Returns
		-----------
		norm_tract: nibabel Nifti1Image
			Normalized tract image.
	"""
	data = tract.get_data()
	total_val = np.sum(data)
	if not total_val == 0: data /= total_val
	
	norm_tract = nib.nifti1.Nifti1Image(data,tract.get_affine())
	return norm_tract
