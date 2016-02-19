#!/usr/bin/env python

### IMPORTS ###
#neuro imaging
import nibabel as nib
from dipy.data import Sphere
from odf_reconstruction import read_reconstruction
#matrix manipulation
import numpy as np
import scipy.spatial as sp
import math
#logging / status information 
from time import time 
import sys
#IO
import argparse, gzip, cPickle
#own imports
from utils import tissue2int, random_sampling
from RCSP import Property, Node, Edge, Graph

#computes all neighbor direction for a 3x3x3 neighborhood
def neighbour_dir():
	#initialize direction matrix (one row for each neighboring direction)
	dir = np.zeros((26,3))
	#iterate through all directions in the interval [-1,1]s 
	counter = 0	
	for x_dir in range(-1,2):
		for y_dir in range(-1,2):
			for z_dir in range(-1,2):
				if not(((x_dir==0) & (y_dir==0)) & (z_dir==0)): #current location is no neighbor
					dir[counter,:] = [x_dir,y_dir,z_dir]
					counter += 1
	return dir

#returns tissue as string
def tissue2str(tissue_id):
	tissue_id = int(tissue_id)
	if tissue_id == 0:
		return 'BG'
	elif tissue_id == 1:
		return 'CSF'
	elif tissue_id == 2:
		return 'GM'
	else:
		return 'WM'

#computes the cosine distance of the rows of b with vector a
def compute_cosine(a,b):
	return np.apply_along_axis(sp.distance.cosine, 1, b, a)

#returns the edge weights for each direction given an ODF fit for the voxel, samples on a unit sphere 
#and an assignment of which sample directions are closest to which neighbor direction
def compute_edge_weight(model_fit,sphere,closest_dir):	
	#initialize weights
	weights = np.zeros(26)
	#compute fODF values for given sample of unit sphere
	values = model_fit.odf(sphere)
	values[np.where(values < 0)] = 0 #clip negative values
	#number of samples
	N = float(np.count_nonzero(values))
	if N == 0.: return weights #make sure ODF contains values larger than 0
	else:
		#contribution of each sample direction 
		a = math.pi / N
		values = values*a #scale ODF values accordingly
		#for each direction sum up the fODF values for all samples closest to that direction
		for i in range(26): weights[i] = np.sum(values[np.where(closest_dir == i)[0]])
			
	return weights
	
#converts an edge weight (given as the sum of the two corresponding voxels) to a distance
def convert2distance(weight):
	if weight == 0.:
		return weight
	else:
		return -math.log(weight/2.) + 10**-10 #small value added if weight/2 = 1 (largest possible weight) to avoid a distance of zero

def construct_brain_graph(brain_mask,segmentation,model_fit):
	""" onstruct the brain graph based on an fODF and tissue segmentation for every voxel.
	    
	    Parameters
 		----------
 		brain_mask: numpy array (mandatory)
	    	Holds the brain mask in a boolean 3D-array.
	    segmentation: numpy array (mandatory)
	    	Holds segmentation information in a float 3D-array. (0=background, 1=CSF, 2=GM, 3=WM)
	    model_fit: Dipy Object (depends on the type)
			Represents the fitted fODF model for each voxel.
	    	
	    Returns
		-----------
		G: dictionary
			The computed brain graph.
			G['nodes'] - nodes of the graph G['nodes'][n] gives the node attributes of node n
			             G['nodes'][n]['x'] - x-coordinate of the corresponding voxel
			             G['nodes'][n]['y'] - y-coordinate of the corresponding voxel
			             G['nodes'][n]['z'] - z-coordinate of the corresponding voxel
			             G['nodes'][n]['tissue'] - tissue based on segmentation
			             G['nodes'][n]['prior'] - prior of the node
			             G['nodes'][n]['region'] - region label of the node
			G['W'] - edges of the graph G['W'][i][j] gives the edge-weight for the edge between node i and node j 
	"""	
	#image dimensions
	xdim, ydim ,zdim = brain_mask.shape
	#extract voxels which are gray matter (GM) or white matter (WM) and belong to the brain
	no_CSF_voxels = segmentation > 1
	graph_voxels = np.where(np.logical_and(brain_mask,no_CSF_voxels)) 
	#compute mapping between node index and 3D matrix index
	node_hash = (xdim*ydim)*graph_voxels[2]+xdim*graph_voxels[1]+graph_voxels[0] 
	node_dict = {node_hash[x]: x for x in range(0,len(node_hash))}
	#initialize graph
	n = len(graph_voxels[0]) #number of nodes 
	G = {}
	#set node attributes
	nodes = {}
	segment = segmentation[graph_voxels]
	for i in range(0,n): 
		nodes[i] = {'x': graph_voxels[0][i], 'y': graph_voxels[1][i], 'z': graph_voxels[2][i], 'tissue': tissue2str(segment[i]), 'prior': 1, 'region': 0}
	#dictionary of weights for each voxel for each of the discrete neighbor-direction
	n_weights = {}
	#sample directions on unit sphere
	sample = random_sampling(1000)
	sphere = Sphere(xyz=sample)
	#compute discrete directions
	dir = neighbour_dir()
	#compute the closest (smallest angle) direction for each sample 
	closest_dir = np.apply_along_axis(compute_cosine, 1, sample, dir)
	c_dir = np.apply_along_axis(np.argsort, 1, closest_dir)[:,0]
	#pre-compute edge weights for each voxel and each direction
	start = time()
	for i in range(0,n):	
	
		if i % 10000 == 0 and i > 0:
			sys.stdout.write('%i nodes processed in %.4f sec \n' % (i,time()-start))
			sys.stdout.flush()
		#get attributes of current voxel
		#coordinates
		x = nodes[i]['x']
		y = nodes[i]['y']
		z = nodes[i]['z']
		#compute edge weights for all possible neighbors
		weights = compute_edge_weight(model_fit[x,y,z], sphere, c_dir)
		n_weights[i] = weights
	#transform weights into [0,1] interval
	max_weight = 0.
	for weights in n_weights.values():
		cur_max = np.max(weights)
		if cur_max > max_weight: max_weight = cur_max
	for node,weights in n_weights.iteritems():
		n_weights[node] = weights / max_weight

	sys.stdout.write("Weights for nodes constructed in %.2f sec\n" % (time() - start))
	sys.stdout.flush()

	#create edges
	start = time()
	W = {}
	for i in range(0,n):
	
		if i % 10000 == 0 and i > 0:
			sys.stdout.write('%i nodes processed in %.4f sec \n' % (i,time()-start))
			sys.stdout.flush()

		#get attributes of current voxel
		#coordinates
		x = nodes[i]['x']
		y = nodes[i]['y']
		z = nodes[i]['z'] 
    	#tissue
		tissue = nodes[i]['tissue']
		#iterate through directions
		for d in range(0,dir.shape[0]):
			#neighbor coordinate
			x_n = x+dir[d,0]
			y_n = y+dir[d,1]
			z_n = z+dir[d,2]
			#check if neighbor voxel is inside the image
			if ((x_n>=0) & (y_n>=0)) & (z_n>=0):
				#get weight for current direction
				weight = n_weights[i][d]
				#hash value of neighbor voxel	
				hash_value = (xdim*ydim)*z_n+xdim*y_n+x_n		
				if hash_value in node_dict: #check whether voxel is valid (GM or WM)		
					j = node_dict[hash_value] #neighbor node index		
					if not( (tissue == 'GM') and (nodes[j]['tissue'] == 'GM')): #no connections between gray matter voxels  				 				
						if i in W: #adjacency list for i exists 					
							if j in W[i]: #edge between i and j exists (set when j was current node)
								W[i][j]+= weight #weight = w(i->j) + w(j->i)
								W[i][j] = convert2distance(W[i][j]) #set the actual distance ( -log(weight/2) )
								if W[i][j] == 0.: del W[i][j] #delete edge if the distance (weight) is 0
							else:
								W[i][j] = weight #weight = w(i->j)
						else: #initiate adjacency list for node i
							W[i] = {j: weight} #weight = w(i->j)
						if j in W: #adjacency list for j exists
							if i in W[j]: #edge between j and i exists (set when i was current node)
								W[j][i]+= weight #weight = w(i->j) + w(j->i)
								W[j][i] = convert2distance(W[j][i]) #set the actual distance ( -log(weight/2) )
								if W[j][i] == 0.: del W[j][i] #delete edge if the distance (weight) is 0
							else: 					
								W[j][i]=weight #weight = w(j->i)
						else: #initiate adjacency list for node i
							W[j] = {i: weight} #weight = w(j->i)

	sys.stdout.write("Edges constructed in %.2f sec\n" % (time() - start))
	sys.stdout.flush()
	
	#delete unconnected nodes
	nodes_temp = {}
	for node in nodes:
		if node in W:
			if len(W[node]) > 0: nodes_temp[node] = nodes[node]
			else: del W[node]
	
	#set graph
	G['W'] = W
	G['nodes'] = nodes_temp
	
	return G

def convert_graph(graph,prior=None,mask=None,atlas=None):
	""" Convert the brain graph from python to C format, optionally including WM probabilities and ROI labels.\n\
	    Note: edge weights of the given graph are assumed to not contain any prior information!
	    
	    Parameters
 		----------
 		graph: Graph in python or binary format. (mandatory)
	    prior: numpy array
	    	Holds the prior probability information in a float 3D-array.
	    mask: numpy array
	    	Holds binary inclusion/waypoint prior information in a float 3D-array.
	    atlas: numpy array
			Holds a region label for each voxel in a float 3D-array.
	    
	    Returns
		-----------
		new_graph: Graph in binary format.
	"""	
		
	#map from current node IDs to [0, number of nodes) interval
	id_to_idx = {}
    #collect nodes
	nodes = []
	#current IDs of the nodes
	if type(graph) is dict: nodelist = graph['nodes'].keys() #python format
	elif type(graph) is Graph: nodelist = range(0,graph.no_of_nodes()) #binary format
	for i, node_id in enumerate(nodelist):
		id_to_idx[node_id] = i
		#python format
		if type(graph) is dict:
			node = graph['nodes'][node_id]
			#node coordinates
			x = node['x']
			y = node['y']
			z = node['z']
			tissue = tissue2int(node['tissue'])
			cur_roi = node['region']
		#binary format
		elif type(graph) is Graph:
			node = graph.node(node_id)
			#node coordinates
			x = int(node.properties['position'][0])
			y = int(node.properties['position'][1])
			z = int(node.properties['position'][2])
			tissue = node.properties['tissue'][0]
			cur_roi = int(node.properties['region'][0])
			
		#ROI label information
		if atlas is not None: roi = atlas[x,y,z] #overwrite current label with new label
		else: roi = cur_roi #keep current label
		#prior information
		if prior is not None: #new prior
			prob = prior[x,y,z]
			#set probability for GM and WM voxels low but above 0 to avoid unconnected nodes
			if prob == 0. and (tissue == 2 or tissue == 3): prob = 10**-10
		else: prob = 1
		#apply binary mask to prior
		if mask is not None: prob *= mask[x,y,z]
		#set the properties for current node
		properties = {
            'position': Property([ float(x), float(y), float(z) ]),
			'tissue': Property(tissue),
			'region': Property(float(roi)),
			'prior': Property(float(prob))
        	}
		nodes.append(Node(i, properties))
    #collect edges 
	edges = []
	for source_id in nodelist:
		sid = id_to_idx[source_id] #new ID of node
		edge_list = [] #edges of current node
		#neighbor nodes
		if type(graph) is dict: target_list = graph['W'][source_id].keys() #python format
		elif type(graph) is Graph: target_list = [e.node for e in graph.edges(source_id)] #binary format
		for target_id in target_list:
			tid = id_to_idx[target_id] #new ID of neighbor node
			#edge weight
			if type(graph) is dict: ew = graph['W'][source_id][target_id] #python format
			elif type(graph) is Graph: ew = graph.edge(source_id,target_id).weight #binary format
			if ((prior is not None) or (mask is not None)): #edges need to be updated if prior changed or mask was applied
				prob1 = nodes[sid].properties['prior'][0]
				prob2 = nodes[tid].properties['prior'][0]
				if (prob1>0 and prob2>0): #not removed by mask
					#include prior of endnodes in the edge weight
					if prior is not None: ew -= math.log(math.sqrt(prob1) * math.sqrt(prob2))
					edge_list.append(Edge(tid, ew)) #keep old edge weight
			else: edge_list.append(Edge(tid, ew)) #keep old edge weight
		edges.append(edge_list)  #add node's edges to list of edges
	
	#create new graph
	new_graph = Graph(nodes, edges)
	return new_graph
		
def main():
	
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
									description='Given a brain mask, a segmentation file\n\
									             and files containing the fODF (model and fit) for each voxel\n\
									             this scripts computes the voxel based brain graph.')
		
	parser.add_argument('-m','--mask', required=True, metavar='<mask_file>',
					    help='Path to the brain mask file.')
	parser.add_argument('-seg', required=True, metavar='<segmentation>',
					    help='Path to the segmentation file which classifies each voxel into CSF, WM or GM.')
	parser.add_argument('-fit', required=True, metavar='<odf_fit>',
					    help='Pickled file containing the fODF fit for each voxel.')
	parser.add_argument('-model', required=True, metavar='<odf_model>',
					    help='Pickled file containing the fODF model.')
	parser.add_argument('-t', '--type', metavar='<type>', required=False, default='dti', choices=['dti','csd','csa'],
					    help='The type of the fODF reconstruction.')
	parser.add_argument('-o','--output', required=True, metavar='<output>',
					    help='Specifies the output file.')
	
	#parse arguments
	args = parser.parse_args()
	#read images
	img = nib.load(args.mask)
	brain_mask = img.get_data().astype(bool)
	img = nib.load(args.seg)
	segmentation = img.get_data()
	#read fODFs
	model_fit = read_reconstruction(args.model,args.fit,brain_mask,type=args.type)
	sys.stdout.write("Data read\n")
	sys.stdout.flush()
	#construct graph
	G = construct_brain_graph(brain_mask, segmentation, model_fit)
	#output
	fp = gzip.open(args.output + '.gz','wb')
	cPickle.dump(G,fp)
	fp.close()
		
if __name__ == "__main__":
    main()
