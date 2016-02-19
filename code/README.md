# Shortest-path tractography (SPT) with spatial priors - code

**Authors:** Niklas Kasenburg, Silas N. Ørting

## Preceding comments

This folder contains the implementation of the shortest-path tractography framework published by [Kasenburg *et al.* (2016)](http://dx.doi.org/10.1016/j.neuroimage.2016.01.031) and code required to produce the inputs for the SPT method.
The folder [python](https://github.com/nkasenburg/SPT_priors/tree/master/code/python) contains all necessary python scripts to run the code. 
However, since the shortest-path computation is implemented in C++ for runtime improvement, it needs to be compiled to produce the binary files and the pyhton-binding (`RCSP.py`, see also [Compiling the C++ code](https://github.com/nkasenburg/SPT_priors/tree/master/code/README.md#Compiling the C++ code)).
The compiled version contained here should work with any Linux OS, but you need to recompile the code for other operating systems.

## Requirements

To be able to run the code you need to have the following python packages installed (on top of the standard packages):
- nibabel ([http://nipy.org/nibabel/](http://nipy.org/nibabel/))
- Dipy ([http://nipy.org/dipy/](http://nipy.org/dipy/))
- NumPy ([http://www.numpy.org/](http://www.numpy.org/))
- SciPy ([http://www.scipy.org/](http://www.scipy.org/))

## Compiling the C++ code

The C++ code can be compiled using the Makefile in this directory (it will both compile the C++ code in `src` and the corresponding python bindings).
Make sure to first run `make clean` to remove the existing binary files before running `make` to compile the code.
To be able to compile the code you must have a version of the [Boost library](http://www.boost.org/) installed on your machine.

## Running the code

All python files that implement a certain functionality (`odf_reconstruction.py`, `brain_graph.py`, `spt.py` and `utils.py`) contain a main function so that they can be run via the command line via the following to options:
- `python script.py [args]`
- `./script.py [args]` (when using the bash shell; make sure that the first line in `script.py` links to your python installation)

Using the argument `-h` prints a description of the script and possible arguments. The general order in which the scripts needs to be run is the following: 

1. `odf_reconstruction.py` to create the fODF model in each voxel from your data
2. `brain_graph.py` to construct the brain graph
3. `spt.py` to perform shortest-path tractography
4. (optional) `utils.py` to post-process the streamline file

In the following sections I will explain the functionality of each script separately.

### Computing the fODF model

The current implementation to construct the brain graph requires the fODFs to be computed using the [DiPy library](http://nipy.org/dipy/). The script `odf_reconstruction.py` is a wrapper to construct the fODFs using one of the following models implemented in DiPy:

- Diffusion tensor
- Constrained spherical deconvolution (CSD) [1]
- Constant solid angle (CSA) [2]

To construct the fODFs you need the following files which describe the diffusion data:

- A 4D Nifti file containing the diffusion weighted images (e.g. `dwi.nii`)
- A text file in [FSL](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) format describing the b-values for each image (e.g. `bvals.txt`)
- A text file in [FSL](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) format describing the b-vector for each image (e.g. `bvecs.txt`)
- (optional) A Nifti file containing a binary mask to define the brain in the image used to only compute the fODF for voxels within the mask (e.g. `mask.nii`)

Run the following line to construct the fODFs:
```
python odf_reconstruction.py -d dwi.nii -bvals bvals.txt -bvecs bvecs.txt -m mask.nii -t csd -order 8 -o csd_order8
```
Here, `-t` gives the type of the fODF (dti, csd or csa), `-order` gives the order to use for CSD or CSA and `-o` defines the name of the output files. 
Running the script as above will create two files: `csd_order8.model` and `csd_order8.fit` containing the fODF model that can be used for the brain graph construction.

### Brain graph construction

To construct the brain graph you need the following files:

- A Nifti file containing a binary mask to define the brain in the image used to only compute the fODF for voxels within the mask (e.g. `mask.nii`)
- A Nifti file containing the tissue labels (0: background, 1: CSF, 2: gray matter, 3: white matter) in diffusion space (e.g. `tissue_seg.nii`)
- Files describing the fODF model (e.g. `csd_order8.model`) and the corresponding fit for every voxel (e.g. `csd_order8.fit`) 

Run the following line to construct the brain graph:
```
python brain_graph.py -m mask.nii -seg tissue_seg.nii -fit csd_order8.fit -model csd_order8.model -t csd -o brain_graph
```
Here, `-t` gives the type of the fODF (dti, csd or csa), needed to poperly read in the fODFs, and `-o` defines the name of the output file. 
Running the script as above will create the output file `brain_graph.gz` which contains the brain graph in a gzipped, pickled file.
The brain graph can now be used to run the SPT.

### Shortest-path tractography (SPT)

To run SPT you need the following files:

- A file containing the brain graph constructed with `brain_graph.py` (e.g. `brain_graph.gz`)
- A Nifti file containing the labels for the seed and target ROI (`atlas.nii`)

#### Running SPT without any prior

To run SPT without any prior run:
```
python spt.py -g brain_graph.gz -a atlas.nii -l 1 2 -wm 1 1 -o tract
```
Here, `-l` defined the two labels that are contained in `atlas.nii` and should be used to define the target and seed region and `-wm` defines whether the white matter voxels of the corresponding region should be included (1) or not (0).

#### Running SPT with priors and/or waypoint masks

Additionaly you can add any number of prior or waypoint mask as separate Nifti files using the following optional arguments in the command line:

- `-p` followed by a space separated list of files that contain prior information. The prior will be constructed by multiplying the different priors.
- `-wpm` followed by a space separated list of files that contain maypoint mask information. The mask will be multiplied to the existing prior after optional modification of the mask (see below).
- `-aswp` followed by a space separated list of the same length as for `-wpm` of 0s and 1s that describe whether the waypoint mask should be used as it is (0) or whether a waypoint slice should be created (1).
- `-dim` followed by a space separated list of the same length as for `-wpm` of 0s, 1s or 2s describing the dimension (0: x, 1: y, 2: z) along which a waypoint mask should be computed if `-aswp` is set to 1 for the corresponding mask.


#### Output files

Running the SPT script as described above will create two different files:

- `tract.nii.gz` containing the confidence value for each voxel (see [3] for further information)
- `tract.Bfloat`, a binary file containing each computed path (see the documentation in `spt.py` for further information)

To convert the binary file containing the paths into a [TrackVis](http://trackvis.org/) (.trk) file you can use `utils.py`:
```
python utils.py -s tract.Bfloat -hdr mask.nii -o tract.trk
```
Here, `-hdr` needs a file in the diffusion space to set the correct header information for the output file.

## References
[1] J-D. Tournier, F. Calamante and A. Connelly, [Robust determination of the fibre orientation distribution in diffusion MRI: Non-negativity constrained super-resolved spherical deconvolution.](http://dx.doi.org/10.1016/j.neuroimage.2007.02.016) *NeuroImage* 35(4): 1459–1472 (2007)

[2] I. Aganj, C. Lenglet, G. Sapiro, E. Yacoub, K. Ugurbil and N. Harel, [Reconstruction of the Orientation Distribution Function in Single and Multiple Shell Q-Ball Imaging within Constant Solid Angle.](http://dx.doi.org/10.1002/mrm.22365) *Magn Reson Med.* 64(2): 554–566 (2010)

[3] N. Kasenburg, M. Liptrot, N. L. Reislew, S. N. Ørting, M. Nielsen, E. Garde, and A. Feragen, [Training shortest-path tractography: Automatic learning of spatial priors.](http://dx.doi.org/10.1016/j.neuroimage.2016.01.031) *NeuroImage* 130:63-76 (2016).
