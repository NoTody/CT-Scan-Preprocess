# Basic commands for displaying DICOM metadata processing

import pandas as pd
import numpy as np
import os
import datetime
import time
import pydicom
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import matplotlib
import time
import copy
from scipy import ndimage
import re

def orientCoord2viewManual(ImageOrientationPatient,orientationDict):
    if ImageOrientationPatient not in list(orientationDict.keys()):
        view='Other'
    else:
        view=orientationDict[ImageOrientationPatient]
    return view   
    
def orientationAxialManualDictionary():
    orientAxial=[
    "['1.000000' '0.000000' '0.000000' '0.000000' '1.000000' '0.000000']",
    "['1' '0' '0' '0' '1' '0']",
    "['1.00000' '0.00000' '0.00000' '0.00000' '1.00000' '0.00000']",
    "['1.0' '-0.0' '0.0' '-0.0' '1.0' '0.0']",
    "['1.0' '0.0' '0.0' '0.0' '1.0' '0.0']",    
    "['1.0000000' '0.0000000' '0.0000000' '0.0000000' '1.0000000' '0.0000000']",
    "['1' '-0.0' '0.0' '-0.0' '1' '0.0']",
    "['1' '0' '0' '0' '1' '-0']"]
    
    orientProne=[
    "['-1' '0' '0' '0' '-1' '0']",
    "['-1.000000' '0.000000' '0.000000' '0.000000' '-1.000000' '0.000000']",
    "['-1.00000' '0.00000' '0.00000' '0.00000' '-1.00000' '0.00000']",
    "['-1.0000000' '0.0000000' '0.0000000' '0.0000000' '-1.0000000' '0.0000000']"]
    
    orientProneToeHead=[
    "['1' '0' '0' '0' '-1' '0']"]
    
    orientationDict={}
    for x in orientAxial:
        orientationDict[x]='Axial'
    for x in orientProne:
        orientationDict[x]='Prone'
    for x in orientProneToeHead:
        orientationDict[x]='ProneToeHead'

    return orientationDict

def orientationGeneral():
    orientation_dictionary={
    'Coronal':"['1', '0', '0', '0', '0', '-1']",
    'Sagittal':"['0', '1', '0', '0', '0', '-1']",
    'Axial':"['1', '0', '0', '0', '1', '0']"}
    
    return orientation_dictionary

def orientCoord2array(ImageOrientationPatient):
    if pd.isna(ImageOrientationPatient):
        IOP_array=[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
    else:
        split_string=ImageOrientationPatient.split('\'')
        IOP_array_str=[split_string[1],split_string[3],split_string[5],
                   split_string[7],split_string[9],split_string[11]]
        IOP_array=[float(x) for x in IOP_array_str]
    return IOP_array 


def orientVector2view(IOP_array,threshhold=0.05):
    orient_name=['Axial','Prone','ProneToeHead']
    
    orientAxial_vector=[1.0, 0.0,0.0,0.0,1.0,0.0]
    orientProne_vector=[-1.0,0.0,0.0,0.0,-1.0,0.0]
    orientProneToeHead_vector=[1.0,0.0,0.0,0.0,-1.0,0.0]
    
    orient_array=[orientAxial_vector,orientProne_vector,orientProneToeHead_vector]
    distance_array=[np.linalg.norm(np.array(x)-np.array(IOP_array)) for x in orient_array]
    
    #distanceAxial=np.linalg.norm(np.array(orientAxial_vector)-np.array(IOP_array))
    #distanceProne=np.linalg.norm(np.array(orientProne_vector)-np.array(IOP_array))
    # distanceProneToeHead=np.linalg.norm(np.array(orientProneToeHead_vector)-np.array(IOP_array))
    min_dist=np.min(distance_array)
    min_idx=np.argmin(distance_array)
    if min_dist<threshhold:
        view=orient_name[min_idx]
    else:
        view='Other'
    return view

def orientCoord2viewVector(ImageOrientationPatient,threshhold=0.05):
    IOP_array=orientCoord2array(ImageOrientationPatient)
    view=orientVector2view(IOP_array,threshhold=threshhold)
    return view


def get_single_dicom_frame_ordered(FilePath):
    dcm_series  = np.sort(os.listdir(FilePath))
    #fname=FilePath + "/" + dcm_series[0]
    fname=os.path.join(FilePath,dcm_series[0])
    ds=pydicom.dcmread(fname)
    return ds

def get_single_dicom_frame(FilePath):
    fname  = next(os.scandir(FilePath)).name
    ds=pydicom.dcmread(os.path.join(FilePath,fname))
    return ds


def estimate_axial_FOV_header(FilePath):
    slice_location_list=get_dcm_slice_loc(FilePath,'all')
    if len(slice_location_list)>0:
        FOV=max(slice_location_list)-min(slice_location_list)
    else: 
        FOV=None
    return FOV
    
def estimate_axial_FOV_nslice(imaging_df):
    FOV = imaging_df.SliceThickness*(imaging_df.NumberOfSeriesRelatedInstances-1)
    return FOV
    
def get_dcm_slice_loc(FilePath,max_frames='all'):
    '''
    @brief: take filepath
    @input: folder location of where DICOM files are located
    @output: image_list: list of fullfile images locations
            slice_location_list: list of slice locations
            image_list_skip: list of fullfile image locations that don't have slice locations
            skipcount: number of images in folder skipped
    
    To do: program in contingencies for cases where there are no slice locations at all
    '''
    slice_location_list = []

    dcm_series  = np.sort(os.listdir(FilePath))
    num_images = len(dcm_series)
    if max_frames=='all':
        max_frames=num_images

    if max_frames < num_images:
        frame_list_temp=list(np.round(np.linspace(0,max_frames-1,max_frames-1).astype('int')))
        image_list_temp=list(map(dcm_series.__getitem__, frame_list_temp))
    else:
        frame_list_temp=list(np.round(np.linspace(0,num_images-1,num_images-1).astype('int')))
        image_list_temp=dcm_series
    
    for dcm_file in image_list_temp:
        fname=FilePath + "/" + dcm_file
        ds=pydicom.dcmread(fname)
        if hasattr(ds, 'SliceLocation'):
            slice_location_list.append(float(ds.SliceLocation))
    
    if len(slice_location_list)>0:
        slice_location_list=sorted(slice_location_list)

    return slice_location_list

def has_repeated_val(temp_list):
    _,c=np.unique(temp_list,return_counts=True)
    if any(c>1):
        has_repeated=True
    else:
        has_repeated=False
    return has_repeated

def has_repeated_slic_loc(FilePath,max_frames='all'):
    slice_location_list=get_dcm_slice_loc(FilePath,max_frames=max_frames)
    has_repeated=has_repeated_val(slice_location_list)
    return has_repeated

def getSingleAttribute(df,attribute):
    ds=get_single_dicom_frame(df.FilePath)
    if hasattr(ds,attribute):
        return ds[attribute].value
    else:
        return None