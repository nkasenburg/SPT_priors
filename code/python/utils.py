#!/usr/bin/env python

### IMPORTS ###
#neuro imaging
import nibabel as nib
#Matrix operations
import numpy as np
#Others
import math, argparse, sys, struct

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
	""" Computes an inclusion mask finding the center of the region along the given dimension 
	    and then setting each voxel in that slice outside the ROI to 0 and all other voxels to 1. 
	    
	    Parameters
 		---------
	    roi_img: numpy 3D-array (mandatory)
	    	Atlas image from which to extract the mask.
		dim: int (Default: 1)
	    	Dimension along which to find the middle slice (0=x,1=y,2=z).
		
		Returns
		-----------
		inclusion_mask: numpy 3D-array
			Inclusion mask.
	"""
	#find the middle slice along given dimension
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
	data = tract.get_data() #get tract values
	total_val = np.sum(data) #sum over all values
	if not total_val == 0: data /= total_val #divide by the sum
	#construct new tract image
	norm_tract = nib.nifti1.Nifti1Image(data,tract.get_affine())
	return norm_tract

def read_streamlines(file):
	"""  Transforms the streamlines contained in the given binary file constructed with SPT
	     into an object that can be saved as TrackVis file with nibabel. 
	
		 Parameters
		 ----------
		 file : string
		    Path to the streamline file.
	"""	
	
	streamlines = [] #list of streamlines
	float_size = 4 #size of the floats
	#read the binary file	
	f = open(file)
	data = f.read()
	f.close()
	#initialize loop variables
	counter = 0
	offset = 0			
	while offset < len(data):
		#streamlines are stored in the file as follows (all entries are 32 bit float):  
		#[N, path score, seed point index, <x_1>, <y_1>, <z_1>, ... , <x_N>, <y_N>, <z_N>]
		#read the path length and score
		tmp = struct.unpack('>3f',data[offset:offset+3*float_size]) 
		N = int(tmp[0])
		score = tmp[1]
		offset += 3*float_size #increase offset
		#read in the voxels of the current streamline	
		sl = [] #current streamline
		for i in range(0,N):
			#read in voxel coordinates 
			values = struct.unpack('>3f',data[offset:offset+3*float_size])
			x = values[0]
			y = values[1]
			z = values[2]
			#append current voxel	
			sl.append([x,y,z])
			#increase offset	
			offset += 3*float_size
		#add current streamline to set of streamlines
		sl = np.array(sl)
		streamlines.append((sl, None, None))
		counter += 1
			
		if counter % 100000 == 0:
			sys.stdout.write('%i streamlines read.\n' % counter)
				
	sys.stdout.write('%i streamlines read.\n' % counter)
	
	return streamlines

def main():
	
	parser = argparse.ArgumentParser(
	 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	 description='Converts paths saved in binary format created by SPT into a TrackVis (.trk) file.'
	)
	parser.add_argument('-s','--streamlines', required=True, metavar='<streamline_file>',
					    help='Path to the binary file containing the paths created by SPT.')
	parser.add_argument('-hdr','--header', required=True, metavar='<header_file>',
					    help='Image for header information.')
	parser.add_argument('-o','--output', required=True, metavar='<output>',
					    help='Path to the output file.')
	#parse command line arguments
	args = parser.parse_args()
	#read in data
	hdr_img = nib.load(args.header) #load header image
	#read streamlines
	streamlines = read_streamlines(args.streamlines)
	#create nibabel TrackVis meta information
	hdr_mapping = {'dim':np.array(hdr_img.shape),'voxel_size':np.absolute(np.diagonal(hdr_img.get_affine())[0:3]),'vox_to_ras':hdr_img.get_affine()}
	#save paths as TrackVis streamline file
	nib.trackvis.write(args.output, streamlines, hdr_mapping=hdr_mapping, points_space='voxel')

if __name__ == "__main__":
    main()

