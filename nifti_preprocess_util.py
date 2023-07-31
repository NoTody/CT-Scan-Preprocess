import skimage.transform as skTrans
import pydicom as dicom
import matplotlib.pylab as plt
import glob

import pandas as pd
import numpy as np
import os
import datetime
import time
import sys
import dask
import pydicom
import copy

import dcm_display_util
import SimpleITK as sitk

# Preprocess one dicom scan folder with nifti
def preprocess(FilePath, OutPath, desired_spacing, desired_size):
    # convert all dicom in one scan to single nifti file
    nifti_img = dicom_to_nifti(FilePath, OutPath)
    
    # Windowing
    window_low, window_high = dcm_display_util.get_window('lung')
    print(f"Window Low: {window_low}, Window High: {window_high}")
    
    # Windowing + Pixel Range Linear Transformation (0-255)
    nifti_img = sitk.IntensityWindowing(nifti_img, window_low, window_high, 0, 255)
    
    # Adjust Voxel Spacing
    nifti_img = spacing_resample(nifti_img, desired_spacing, desired_size)
    
    # Adjust Shape
    new_arr = skTrans.resize(sitk.GetArrayFromImage(nifti_img), desired_size, order=1, preserve_range=True)
    nifti_img = sitk.GetImageFromArray(new_arr)
        
    print(f"Resampled Spacing: {nifti_img.GetSpacing()}, Resampled Size: {nifti_img.GetSize()}")
    
    print(f"Writing NIfTI to file {OutPath} ...")
    sitk.WriteImage(nifti_img, OutPath)
    
    return nifti_img


# read all dicom files into nifti format. Scale/Intercept is converted automatically.
# (https://github.com/InsightSoftwareConsortium/ITK/blob/75babc2b45ba0f8f2f60d6cdb5e8a0ffcc1e77fa/Modules/IO/NIFTI/src/itkNiftiImageIO.cxx#L736)

def dicom_to_nifti(FilePath, OutPath):
    print("Reading Dicom directory:", FilePath)
    reader = sitk.ImageSeriesReader()

    dicom_names = reader.GetGDCMSeriesFileNames(FilePath)
    reader.SetFileNames(dicom_names)

    image = reader.Execute()

    size = image.GetSize()
    
    print("Image size:", size[0], size[1], size[2])
    
    return image


# Create a resampling object with the desired output spacing
def spacing_resample(nifti_img, desired_spacing, desired_size):
    original_spacing = nifti_img.GetSpacing()
    original_size = nifti_img.GetSize()

    new_size = [int(round(original_size[0] * original_spacing[0] / desired_spacing[0])),
                int(round(original_size[1] * original_spacing[1] / desired_spacing[1])),
                int(round(original_size[2] * original_spacing[2] / desired_spacing[2]))]

    # Perform the resampling
    nifti_img = sitk.Resample(nifti_img, new_size, sitk.Transform(), sitk.sitkBSpline, nifti_img.GetOrigin(),
                            desired_spacing, nifti_img.GetDirection(), 0.0, nifti_img.GetPixelID())

    return nifti_img
