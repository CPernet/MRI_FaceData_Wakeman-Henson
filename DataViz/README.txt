'Standard' analysis of Wakeman and Henson fMRI face dataset

code/DatViz_fMRI_analysis.m is the code to process the dataset. Note we assume data are stored using the Brain Imaging Data Structure. The analysis differs from the reported one as we included slice timing, changed the number of Gaussians in the segmentation and removed subject 15.

code/DatViz_fMRI_analysis.m is the code used to generate the mapping in subject space (ie being lazy we simply map back the normalized statistical maps to subject space) and to create the average T1w image: subjects_anat_avg.nii
