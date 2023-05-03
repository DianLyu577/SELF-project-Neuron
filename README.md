# SELF-project-Neuron: Causal Evidence for the Processing of Bodily Self in the Anterior Precuneus
---

The repository contains the group-based ROI of self-hot electrodes and self-hot electrodes in the posteromedial cortex, and some key customized codes for generating the results of paper.


## Description of the data and file structure

The two nifti files are group-based ROIs from the study. With brain stimulation, we found a distinct area in the anterior precuneus that can induce robust subjective changes that are related to self-location displacement (sometimes self-dissociation). These stimulation sites are labelled "hot", while the surrounding posteromedial (PMC, or otherwise addressed as posterior cingulate cortex/precuneus) sites that did not induce such an effect are labelled "cold". The method of extracting the native coordinates of these stimulation sites is presented in the extractSbjNativeCoord.m. By transforming the individual brain to the standardized MNI space, we generated the presented ROIs by summarizing the spatial location of these stimulation sites (using 4mm-radius sphere). The ROIs are plotted in the provided figure ROI_vis.png, where the yellow color indicates hot sites. We hope these ROIs will be helpful for people who are interested in studying self-related processing, PMC heterogeneity and many other. 

Among the *.m files, the fMRINetMapNative.m is the customized class function for doing all the fMRI-related processings, including preprocessing and seed-based FC at the native space as reported in the paper. The usage of this class can be found in the FCconndiff_selfvscold_example.m, which is the code that is used to conduct the seed-based FC of hot/cold stimulation sites and the contrast between them at the individual level. The group-level comparison (including normalizing the first-level stats images) is presented in the Elec_ContrElec_glm2.m. The results of this analysis was reported in the Figure 2 of the paper. 


## Sharing/Access information

Links to other publicly accessible locations of the data:
  * https://github.com/DianLyu577/SELF-project-Neuron/tree/main



## Code/Software
MATLAB >= 2019a
SPM12 
