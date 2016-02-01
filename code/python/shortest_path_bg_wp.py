#!/usr/bin/env python

### IMPORTS ###

# IO
import sys, struct
# Neuroimage
import nibabel as nib
# Other
from time import time
import pdb
#numpy
import numpy as np
#math
import math
# Own Imports
from RCSP import Graph, shortest_path_between_rois
from brain_graph_sample import convert_graph
from utils import tissue2int, compute_inclusion_mask
# Parser
import argparse

def find_ROI_nodes(graph,atlas,roi_label,use_wm=[False,False]):
	""" Computes the nodes set for each ROI.
	    
	    Parameters
 		----------
 		graph: RCSP Graph (mandatory)
	    	Graph in which to perform the tractography.
	    atlas: numpy float array (mandatory)
	    	3-D array holding the labels for the ROIs in voxel coordinates.
	    label: [int,int] (mandatory)
			Label of the ROIs.
	    use_wm: [wm,wm] where wm is boolean
	    	Whether to use WM voxels (nodes) of the ROIs. (Default: [False,False])
	    
	    Returns
		-----------
		roi1_nodes: [node1,node2,...]
			List of node IDs for first ROI.
		roi2_nodes: [node1,node2,...]
			List of node IDs for second ROI.
	"""
	
	n = graph.no_of_nodes()
	
	#select target, seed and waypoint voxels
	roi1_nodes = []
	roi2_nodes = []
	for node in range(0,n):
		tissue = graph.node(node).properties['tissue'][0]
		x = int(graph.node(node).properties['position'][0])
		y = int(graph.node(node).properties['position'][1])
		z = int(graph.node(node).properties['position'][2])
		if  tissue == float(tissue2int('GM')):
			if atlas[x,y,z] == roi_label[0]:
				roi1_nodes.append(node)
			if atlas[x,y,z] == roi_label[1]:
				roi2_nodes.append(node)
		elif tissue == float(tissue2int('WM')):
			if atlas[x,y,z] == roi_label[0] and use_wm[0]:
				roi1_nodes.append(node)
			if atlas[x,y,z] == roi_label[1] and use_wm[1]:
				roi2_nodes.append(node)
	
	#make sure to compute shortest paths for smaller region
	if len(roi1_nodes) > len(roi2_nodes):
		temp = roi1_nodes
		roi1_nodes = roi2_nodes
		roi2_nodes = temp
	
	return roi1_nodes, roi2_nodes
	
def SPT(graph,roi_nodes):
	""" Computes the tracts using shortest-path tractography on a graph.
	    
	    Parameters
 		----------
 		graph: RCSP Graph (mandatory)
	    	Graph in which to perform the tractography.
	    roi_nodes: ([roi1_node1,...],[roi2_node2,...],[wp_node1,...])
	    	2-tuple of lists of the node IDs for the two ROIs.
	    	
	    Returns
		-----------
		roi_to_roi: RSCP ROItoROI
			Contains all shortest paths.
	"""	
	
	roi1_nodes = roi_nodes[0]
	roi2_nodes = roi_nodes[1]
	#run tractography
	roi_to_roi = shortest_path_between_rois(graph, roi1_nodes, roi2_nodes)
	
	return roi_to_roi	
	 			
def create_tract_image(graph,roi_to_roi,header):
	""" Computes a tract image from tracts.
	    
	    Parameters
 		---------
 		graph: RCSP Graph (mandatory)
	    	Graph in which tractography was performed.
 		roi_to_roi: RSCP ROItoROI (mandatory)
			Contains all shortest paths.
	    header: nibabel Nifti1Image (mandatory)
	    	Image for header information.
		
		Returns
		-----------
		tract_img: nibabel Nifti1Image
			Image with the tract.
	"""
	
	# This will update the graph so nodes has 'count' and 'weight' properties
	graph.calculate_node_importance(roi_to_roi,'wmprob')
	# We need the max size of the dimensions to create the matrix
	rows, columns, slices = header.get_data().shape
	# This will create an array representing the 3D matrix in row-major order (C order)
	matrix = graph.as_matrix(rows, columns, slices,'weight')
	# Convert to numpy array and reshape
	img_array = np.array(matrix).reshape([rows, columns, slices])
	#convert into nifti image
	tract_img = nib.nifti1.Nifti1Image(img_array,header.get_affine())
	
	return tract_img

def create_streamline_file(roi_to_roi,roi_nodes,graph,out_file):
	""" Computes a tract image from tracts.
	    
	    Parameters
 		---------
 		roi_to_roi: RSCP ROItoROI
			Contains all shortest paths.
		roi_nodes: ([roi1_node1,...],[roi2_node2,...]) (mandatory)
	    	2-tuple of lists of the node IDs for the ROIs.
	    graph: RCSP Graph (mandatory)
	    	Graph in which the tractography was performed.
		out_file: string (mandatory)
			Path to the output file.
	"""
	
	roi1_nodes = roi_nodes[0]
	roi2_nodes = roi_nodes[1]
	
	#output of paths in streamline binary format
	#data is stored like this (all entries are 32 bit float) for each path:  
	#[N, score, seed point index, <x_1>, <y_1>, <z_1>, ... , <x_N>, <y_N>, <z_N>]
	f = open(out_file,'w')
	for n1 in roi1_nodes:
		for n2 in roi2_nodes:
			#get path between the current ROI nodes
			path = roi_to_roi(n1,n2)
			length = len(path)-1 #number of edges on the path
			if length > 0:
				#subtract the log of the square root of the source and target node priors from the path log-score 
				source_weight = math.log(math.sqrt(graph.node(path[0].id).properties['wmprob'][0]))
	 			target_weight = math.log(math.sqrt(graph.node(path[length].id).properties['wmprob'][0]))
				distance = path[length].weight - source_weight - target_weight
				#normalize by path length and transform distance to score 
				avg_length = math.exp(-(distance / float(length)))
				#write path information into file		 
				f.write(struct.pack('>3f',length+1, avg_length, 0.))
				#write coordinates of voxels along the path into file		 
				for node in path:
					ID = int(node.id)
					x = int(graph.node(ID).properties['position'][0])
					y = int(graph.node(ID).properties['position'][1])
					z = int(graph.node(ID).properties['position'][2])
					#update streamline file
					voxel = np.array( [x,y,z], dtype=np.float32)
					f.write(struct.pack('>3f',voxel[0],voxel[1],voxel[2]))		 						
	f.close()

def main():

	parser = argparse.ArgumentParser(
	 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	 description='\Given a brain graph, an atlas files with corresponding label of the ROIs and possible priors or waypoint masks\n\
	              this script computes the tracts and creates an image file\n\
	              and a streamline file in binary format in the given output directory.'
	)
	parser.add_argument('-g','--graph', required=True, metavar='<graph_file>',
					    help='Path to the graph file.')
	parser.add_argument('-a','--atlas', required=True, metavar='<atlas_file>',
					    help='Path to the atlas file.')
	parser.add_argument('-l','--label', required=True, nargs=2, type=int,
					    help='Label of the seed regions of interest in the corresponding atlas.')
	parser.add_argument('-wm','--include-wm', nargs=2, default=[0,0], type=int, choices=[0,1], 
					    help='Whether to use white matter voxels for the corresponding regions.')
	parser.add_argument('-p','--prior', nargs='+', required=False, 
					    help='Additional priors to use (supplied as image files).')
	parser.add_argument('-wpm','--wp-masks', nargs='+', required=False, 
					    help='Additional exclusion masks to use (supplied as image files).')
	parser.add_argument('-aswp','--as_wp', nargs='+', type=int, choices=[0,1], 
					    help='Whether to use the mask as it is or compute the middle slice (needs to be the same number of arguments as masks supplied).')
	parser.add_argument('-dim','--mask-dims', nargs='+', type=int, choices=[0,1,2], 
					    help='Dimensions/axis along which the slices for inclusion/waypoint masks should be chosen (needs to be the same number of arguments as masks supplied).')
	parser.add_argument('-o','--out', required=True, metavar='<outpath>',
					    help='Basename of outputfiles.')
	#read in arguments
	args = parser.parse_args()	
	#whether to use white matter voxels
	use_wm = [False,False]
	use_wm[0] = args.include_wm[0] == 1
	use_wm[1] = args.include_wm[1] == 1	
	#output files
	base_name = args.out
	tract_file = base_name+'.nii.gz'
	streamline_file = base_name+'.Bfloat'
	#atlas for target and seed region
	img = nib.load(args.atlas)
	atlas = img.get_data()
	sys.stdout.write("Atlas read.\n")
	sys.stdout.flush()
	#load graph and add possible priors
	graph = Graph.load_binary(args.graph)
	#create priors by multiplication
	if args.prior:
		prior = np.ones(atlas.shape)
		for file in args.prior: prior *= nib.load(file).get_data()
	else: prior = None
	if args.wp_masks:
		mask = np.ones(atlas.shape)
		index = 0
		for file in args.wp_masks:
			#computes the middle slice of the waypoint based on the given dimension
			if args.as_wp[index]: cur_mask = compute_inclusion_mask(nib.load(file).get_data(),dim=args.mask_dims[index])
			#take mask as it is
			else: cur_mask = nib.load(file).get_data()
			mask *= cur_mask
			index += 1
	else: mask = None
	#add priors into the graph
	if ((prior is not None) or (mask is not None)):  graph = convert_graph(graph,prior=prior,mask=mask) 
	sys.stdout.write("Graph read.\n")
	sys.stdout.flush()
	#collect ROI nodes
	start = time()
	roi1_nodes, roi2_nodes = find_ROI_nodes(graph, atlas, args.label, use_wm)
	sys.stdout.write('ROI voxels identified in %.4f sec\n' % (time()-start) )
	sys.stdout.write('%i roi1 voxels\n' % len(roi1_nodes) )
	sys.stdout.write('%i roi2 voxels\n' % len(roi2_nodes) )
	sys.stdout.flush()
	
	#SPT
	roi_nodes = (roi1_nodes,roi2_nodes)
	start = time()
	roi_to_roi = SPT(graph, roi_nodes)
	sys.stdout.write('Shortest paths computed in %.4f sec\n' % (time()-start) )
	sys.stdout.flush()	
	#create output
	start = time()
	tract_img = create_tract_image(graph, roi_to_roi, img)
	nib.save(tract_img,tract_file)
	create_streamline_file(roi_to_roi, roi_nodes, graph, streamline_file)
	sys.stdout.write('Output created in %.4f sec\n' % (time()-start) )
	sys.stdout.flush()
	
if __name__ == "__main__":
    main()
