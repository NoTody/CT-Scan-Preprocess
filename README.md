CT-Scan-Preprocess
==============================
Introduction
------------
This repo gives an example of the general pipeline for preprocessing CT scan in NifTi format.

Requirements
------------
Check requirements.txt for all dependencies
To install all dependencies with new environment, run
- python -m venv <path_new_pyvenv>
- source <path_new_pyvenv>/bin/activate
- pip install -r requirements.txt

Files:
------------
- dcm_display_util.py: contains code for displaying dicom
- dcm_metadata_util.py: contains code for analyzing meta info for dicom
- nifti_preprocess_util.py: contains code for preprocessing scans stored with dicom files to nifti format
- usage_example.ipynb: gives an example on how to use the functions to preprocess one scan.
