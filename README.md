# NoddiSurfaceMapping
Optimized cerebral cortical surface mapping of neurite properties using diffsion magnetic resonace imaging (dMRI)

### Example NODDI surface maps (NDI and ODI) compared with myelin and cortical thickness
![noddi](https://user-images.githubusercontent.com/16514166/40781643-a102ca94-6517-11e8-81e1-82e9556199de.png)

### Installation and Usage
1. Download NoddiSurfaceMapping.zip and unzip
2. Configure variables in the section of Setup in NoddiSUrfaceMapipng.sh depending on your local settings including AMICO. Note that you need to configure default settings of parallel diffusivity value in AMICO. By default, it is set to 1.7E-3 mm^2/s but you need to change the value to 1.1E-3 mm^2/s for NODDI cortical surface mapping (see Fukutomi et al., Neuroimage 2018).
3. Run NoddiSurfaceMapping.sh - inputs should be \<StudyFolder\> and \<Subject directory(s)\>, which were used to preprocess with HCP pipeline/DiffusionPreprocessing

```
 Calculate NODDI and do surface-mapping for HCP data
 Usage: NoddiSurfaceMapping <StudyFolder> <SubjectID 1> <SubjectID 2> ...
 
 Options:
     -a <num> : species atlas (0. Human [default], 1. Macaque, 2. Marmoset)
     -t <num>,<num>,<num> : b-value upper and lower threshold, and b=0 upper threshold (default: 3100,100,50)
     -M       : REGNAME=MSMAll (default: MSMSulc)
     -s       : do not calculate NODDI but only perform surface mapping
```

### Dependencies
[AMICO][], [HCP pipeline][], [Workbench][], [FSL][]

[AMICO]: https://github.com/daducci/AMICO "AMICO"
[HCP pipeline]: https://github.com/Washington-University/Pipelines "HCP pipeline"
[Workbench]: https://github.com/Washington-University/workbench "Workbench"
[FSL]: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki "FSL"

### References
Fukutomi, H., Glasser, M.F., Zhang, H., Autio, J.A., Coalson, T.S., Okada, T., Togashi, K., Van Essen, D.C., Hayashi, T., 2018. Neurite imaging reveals microstructural variations in human cerebral cortical gray matter. Neuroimage. [DOI][] [BALSA][]

[DOI]: https://doi.org/10.1016/j.neuroimage.2018.02.017
[BALSA]: https://balsa.wustl.edu/study/show/k77v
