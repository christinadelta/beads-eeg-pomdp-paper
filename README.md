# beads-eeg-pomdp-paper

## Data and code related to the manuscript: 

### *From Sampling to Stopping: The P300 ERP component and beta power contribute to reward-related decision commitments*

by Christina Dimitriadou and Nicholas Furl

OSF Preregistration DOI: https://doi.org/10.17605/OSF.IO/NFPJM

All avalyses were conducted using Matlab 2021a. For behavioural analyses and modelling we used custom scripts and functions. For EEG preprocessing and analysis we used the SPM toolbox:
https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ Henson, R. N., Abdulrahman, H., Flandin, G., & Litvak, V. (2019). Multimodal integration of M/EEG and f/MRI data in SPM12. Frontiers in neuroscience, 13, 300.

### Requirements to run the analysis codes:
* Matlab_R2016a or later
* SPM 12
* Fieldtrip

### General information about data, and code 
* Raw EEG and behavioural data are provided in *raw_data/behav/* and *raw_data/eeg/* directories.
* For preprocessing of the EEG data we used SPM12 with custom scripts and functions that can be found in *code/EEG/spm_analysis/* directory.
* For Mass Univariate Analyses and for Individual Differences Analyses we used the SPM EEG user interface.
* All EEG figures were generated using the SPM EEG user interface and manually extracted for visualisation. 
* All code related to behavioural analysis, model fitting, parameter recovery and model recovery can be found in *code/behav/* directory
* Analysis code and code to generate Figure 3 can be found in *code/behav/full_model_recovery_beads_2025.m*
* Analysis code and code to generate Figure 4 and Supplementary figure S1 can be found in *code/behav/fit_MDPBeads_fmincon.m*
* Figures 5-8 and Supplementary Figure 7 were generated using the SPM EEG user interface and manually extracted for visualisation.
* Analysis code and code to generate Supplementary Figure S2 and S3 can be found in *code/behav/fit_MDPBeads_paramRec.m*
* Analysis code and code to generate Supplementary Figure S4 can be found in *code/behav/fit_MDPBeadsCost_paramRec.m*
* Analysis code and code to generate Supplementary Figure S5 can be found in *code/behav/fit_MDPBeadsNoise_paramRec.m*
* Analysis code and code to generate Supplementary Figure S6 can be found in *code/behav/sampling_nll_relationship.m*

### General information about the experimental task
* All scripts and functions to run the beads task are available in /information_sampling_eeg_beads_paper/beads_task/
* You should be able to run the task through /information_sampling_eeg_beads_paper/beads_task/ using *startup_script.m*
* The results direcory is created automatically the first time you run the task. All data are store in results, in seperate directories for each participant

### Requirements to run the task:
* The psychophysics toolbox Version 3.0.19 or later (http://psychtoolbox.org)
* If you run the task with EEG via a parallel port, you will need to install a few configuration files. The files and steps of installation can be found here: https://github.com/christinadelta/personal_notes/blob/main/install_io64mex_file_windows.md
