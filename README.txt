This folder contains the analysis pipeline for the manuscript "Ear-EEG Measures of Auditory Attention to Continuous Speech" (https://doi.org/10.3389/fnins.2022.869426).

The corresponding dataset can be found in the OpenNeuro repository https://openneuro.org/datasets/ds004015.

To replicate the analysis use MATLAB R2019b and follow these steps:

1) Download all scripts into one folder.
2) Download the BIDS dataset from https://openneuro.org/datasets/ds004015 into a folder that is on the same level as the folder containing the scripts. 
3) Download EEGLAB v2020.0 and add the respective path within bjh_00_0_main_cEEGrid.m. 
4) Download the EEGLAB plugins "cEEGrid v. 0.9" (https://doi.org/10.5281/zenodo.5946875), "clean_rawdata2.4", and "Fieldtrip-lite20210311" into the plugin folder of EEGLAB v2020.0.
5) Download the mTRF Toolbox v. 2.1 (https://github.com/mickcrosse/mTRF-Toolbox/releases/tag/v2.1) and add the respective path within bjh_00_0_main_cEEGrid.m.
6) Within bjh_00_0_main_cEEGrid.m specify the name of the folder in which you downloaded the data (BIDS_folder_name).
7) Run bjh_00_0_main_cEEGrid.m from within the folder containing all scripts. As a result, a folder called data will be created on the same level as the BIDS and script folder in which intermediate data and figures will be stored (the intermediate data will amount to 17 GB, this process takes several hours).
8) After you completely ran the analysis you can run specific parts again and skip others by changing the config structure within bjh_00_0_main_cEEGrid.m (1 means the analysis step will be run, 0 means it will be skipped).

If there are any remaining questions do not hesitate to contact me at bjoern.holtze[at]uol.de.
