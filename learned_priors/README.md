## Learned priors

**Authors:** Niklas Kasenburg,

### Description

This folder contains the learned priors created from shortest-path tractography (SPT) results on data provided by the Human Connectome Project (HCP) [1,2]. Creation and application of the learned priors is described in detail in [Kasenburg *et al.* (2016)](http://dx.doi.org/10.1016/j.neuroimage.2016.01.031). In the following there is a brief overview of the priors included in this folder:
- Cortico spinal tract (CST): The learned prior was created from tractography results that included a subject-specific white matter prior, as described in [3].
- Arcuate fasciculus (AF): The learned prior was created from tractography results that included a subject-specific white matter prior, as described in [3].
- Inferior fronto-occipital fasciculus (IFOF): The learned prior was created from tractography results that included a subject-specific white matter prior and two waypoint priors (occipital lobe and external capsule), as described in [3].
- Fornix: The learned prior was created from tractography results that included a subject-specific white matter prior, a waypoint priors (superior right and left branch of the reference [4]) and white matter tract atlas [5], as described in [3].
  
### HCP subjects

The following 38 HCP subjects were used to create the data: 
- 100307
- 103414
- 105115
- 110411
- 111312
- 113619
- 115320
- 117122
- 118730
- 118932
- 123117
- 124422
- 125525
- 128632
- 129028
- 130013
- 133928
- 135932
- 136833
- 139637
- 149337
- 149539
- 151223
- 151627
- 156637
- 161731
- 192540
- 201111
- 212318
- 214423
- 221319
- 298051
- 397760
- 414229
- 499566
- 654754
- 672756
- 792564

### References

[1] D. C. V. Essen, S. M. Smith, D. M. Barch, T. E. J. Behrens, E. Yacoub and K. Ugurbil for the WU-Minn HCP Consortium, [The WU-Minn Human Connectome Project: An overview.](http://dx.doi.org/10.1016/j.neuroimage.2013.05.041) *NeuroImage* 80: 62–79 (Oct 2013)

[2] M. F. Glasser, S. N. Sotiropoulos, J. A. Wilson, T. S. Coalson, B. Fischl, J. L. Andersson, J. Xu, S. Jbabdi, M. Webster, J. R. Polimeni, D. C. V. Essen and M. Jenkinson for the WU-Minn HCP Consortium, [The minimal preprocessing pipelines for the Human
Connectome Project.](http://dx.doi.org/10.1016/j.neuroimage.2013.04.127) *NeuroImage* 80:105–124 (Oct 2013)

[3] N. Kasenburg, M. Liptrot, N. L. Reislew, S. N. Ørting, M. Nielsen, E. Garde, and A. Feragen, [Training shortest-path tractography: Automatic learning of spatial priors.](http://dx.doi.org/10.1016/j.neuroimage.2016.01.031) *NeuroImage* (2016).

[4] M.T. de Schotten, D. H. Ffytche, A. Bizzi, F. Dell’Acqua, M. Allin, M. Walshe, R. Murray, S. C. Williams, D. G. M. Murphy and M. Catani, [Atlasing location, asymmetry and inter-subject variability of white matter tracts in the human brain with MR diffusiontractography.](http://dx.doi.org/10.1016/j.neuroimage.2010.07.055) *NeuroImage* 54(1):49–59 (Jan 2011).

[5] A. Varentsova, S. Zhang and K. Arfanakis, [Development of a high angular resolution diffusion imaging human brain template.](http://dx.doi.org/10.1016/j.neuroimage.2014.01.009), *NeuroImage* 91:177–186 (May 2014).
