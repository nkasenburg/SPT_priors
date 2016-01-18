#!/usr/bin/env python
#DiPy IO
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.io.pickles import save_pickle, load_pickle
from dipy.data import Sphere
#DiPy reconstruction models
from dipy.reconst.csdeconv import auto_response, ConstrainedSphericalDeconvModel
from dipy.reconst.shm import CsaOdfModel, SphHarmFit
from dipy.reconst.dti import TensorModel, TensorFit
from dipy.reconst.multi_voxel import MultiVoxelFit
from dipy.reconst.peaks import peaks_from_model
#Own imports
from utils import random_sampling
#Others
import sys
from time import time 
import numpy as np
import nibabel as nib
#IO
import argparse

def reconstruction(dwi,bval_file,bvec_file,mask=None,type='dti',b0=0.,order=4):
	""" Uses Dipy to reconstruct an ODF for each voxel.
	    
	    Parameters
 		----------
 	    dwi: numpy array (mandatory)
	    	Holds the diffusion weighted image in a 4D-array (see nibabel).
	    	
	    bval_file: string (mandatory)
	    	Path to the b-value file (FSL format).
	    
	    bvec_file: string (mandatory)
	    	Path to the b-vectors file (FSL format).
	    
	    mask:  numpy array
	    	Holds the mask in a 3D array (see nibabel).
	    	
	    type: string \in {'dti','csd','csa'} (default = 'dti')
	    	The type of the ODF reconstruction.
		
		b0: float (default = 0)
			Threshold to use for defining b0 images.
	    
	    order: int (default = 4)
	    	Order to use for constrained spherical deconvolution (csd) or constant solid angle (csa).
	    	
	    Returns
		-----------
			
		model_fit: Dipy Object (depends on the type)
			Represents the fitted model for each voxel.
	"""	
	
	#b-values and b-vectors
	bvals, bvecs = read_bvals_bvecs(bval_file,bvec_file)
	gtab = gradient_table(bvals, bvecs, b0_threshold=b0)
	
	#reconstruction
	if type == 'dti':
		model = TensorModel(gtab,fit_method='WLS')
	elif type == 'csd':
		response, ratio = auto_response(gtab, dwi, roi_radius=10, fa_thr=0.7)
		model = ConstrainedSphericalDeconvModel(gtab, response, sh_order=order)
	elif type == 'csa':
		model = CsaOdfModel(gtab, order)

	if mask is not None:
		model_fit = model.fit(dwi,mask=mask)
	else:
		model_fit = model.fit(dwi)
	
	return model_fit

def read_reconstruction(model_file,fit_file,mask,type='dti'):
	""" Reads in an ODF for each voxel reconstructed using DiPy.
	    
	    Parameters
 		----------
 		model_file: string (mandatory)
	    	Path to the pickled model used for ODF reconstruction.
	    	
	    fit_file: string (mandatory)
	    	Path to the pickled, fitted ODFs reconstructed by DiPy.
	    
	    mask: numpy array (mandatory)
	    	Logical array defining the brain mask.
	    	
	    type: string \in {'dti','csd','csa'} (default = 'dti')
	    	The type of the ODF reconstruction.
	    	
	    Returns
		-----------
			
		model_fit: Dipy Object (depends on the type)
			Represents the fitted ODF for each voxel.
	"""	
	
	model = load_pickle(model_file)
	if type == 'dti':
		model_params = load_pickle(fit_file)
		model_fit = TensorFit(model,model_params)
	elif type == 'csd':
		fit_array = load_pickle(fit_file)
		model_fit = MultiVoxelFit(model,fit_array,mask)
	elif type == 'csa':
		shm_coeff = load_pickle(fit_file)
		model_fit = SphHarmFit(model,shm_coeff,mask)
	
	return model_fit
	
def main():
	
	#MAIN PARSER
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	 	description='\Script for ODF-reconstruction.'
	)
	subparsers = parser.add_subparsers()
	# ODF RECONSTRUCTION
	parser_odf_reconst = subparsers.add_parser('reconstruct', formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	 									        description='\Given the DWI data, the corresponding b-value and b-vecs file\n\
	              										     and a the type of reconstruction. This scripts\n\
	              										     computes the ODF reconstruction for every voxel.'
	)
	parser_odf_reconst.add_argument('-d', '--dwi', required=True, metavar='<dwi>',
					    help='Path to the diffusion weighted image file (4D in NIFTI format).')
	parser_odf_reconst.add_argument('-bvals', required=True, metavar='<bvals>',
					    help='Path to the b-value file (FSL format).')
	parser_odf_reconst.add_argument('-bvecs', required=True, metavar='<bvecs>',
					    help='Path to the b-vectors file (FSL format).')
	parser_odf_reconst.add_argument('-m','--mask', metavar='<mask_file>',
					    help='Path to the brain mask file.')
	parser_odf_reconst.add_argument('-t', '--type', metavar='<type>', default='dti', choices=['dti','csd','csa'],
					    help='The type of the ODF reconstruction (default dti).')
	parser_odf_reconst.add_argument('-b0', metavar='<b0>', type=float, default=0.,
					    help='Threshold to use for defining b0 images.')
	parser_odf_reconst.add_argument('--order', type=int, metavar='<order>', default=4,
					    help='Order of ODF (not used for DTI).')
	parser_odf_reconst.add_argument('-o','--output', required=True, metavar='<output>',
					    help='Specifies the output file.')
	parser_odf_reconst.set_defaults(command='reconstruct')
	
	args = parser.parse_args()

	if args.command == 'reconstruct':
		#DWI
		img = nib.load(args.dwi)
		data = img.get_data()
	
		#mask
		if args.mask:
			img = nib.load(args.mask)
			mask = img.get_data().astype(bool)
			mask = np.logical_not(mask)
	
			data[mask,:] = 0
	
		sys.stdout.write('Data read.\n')
		
		start = time()
		if args.mask:
			model_fit = reconstruction(data, args.bvals, args.bvecs, mask=np.logical_not(mask), type=args.type, b0=args.b0, order=args.order)
		else:
			model_fit = reconstruction(data, args.bvals, args.bvecs, type=args.type, b0=args.b0, order=args.order)
		elapsed = time() - start
		
		sys.stdout.write('Fitted %s model with order %i in %.2f sec.\n' % (args.type,args.order,elapsed))
		sys.stdout.flush()
	
		if args.type == 'dti':
			save_pickle(args.output+'.fit',model_fit.model_params)
		elif args.type == 'csd':
			save_pickle(args.output+'.fit',model_fit.fit_array)
		elif args.type == 'csa':
			save_pickle(args.output+'.fit',model_fit.shm_coeff)
		
		save_pickle(args.output+'.model',model_fit.model)
		
if __name__ == "__main__":
    main()
