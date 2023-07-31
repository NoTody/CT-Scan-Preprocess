# Basic commands for displaying DICOM

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
#from IPython.display import HTML

class dcm_visualization_video:
    def __init__(self,output_location=os.path.abspath(os.path.join('..','Test.avi')),save_file=False,fps=30,maxVideoSize=100):
        self.output_location=output_location
        self.save_file=save_file
        #self.codec=codec
        self.fps=fps
        self.maxVideoSize=maxVideoSize
        
class dcm_video_singlescan(dcm_visualization_video):
    def __init__(self,DCMFolderPath,max_frames=40,DisplayAccessionNumber=True,window='lung',output_location=os.path.abspath(os.path.join('..','Test.avi')),save_file=False,fps=30,maxVideoSize=100):
        super().__init__(output_location=output_location,save_file=save_file,fps=fps)
        self.FilePath=DCMFolderPath
        self.window=window
        self.DisplayAccessionNumber=DisplayAccessionNumber
        self.AccessionNumber=None
        
        matplotlib.rcParams['animation.embed_limit']=self.maxVideoSize
        image_list,slice_location_list,_,_ = return_dcm_list_sorted(DCMFolderPath)
        num_images=len(image_list)
        if max_frames=='all':
            self.max_frames=num_images
        else:
            self.max_frames=max_frames

        self.num_images_folder=num_images
        self.slice_location_list=slice_location_list
        
        if num_images<=self.max_frames:
            self.image_list=image_list
            self.slice_location_list=slice_location_list
        else:
            indices=list(np.round(np.linspace(0,num_images-1,self.max_frames)).astype('int'))
            self.image_list=list(map(image_list.__getitem__, indices))
            self.slice_location_list=list(map(slice_location_list.__getitem__, indices))
            #self.image_list=[image_list[i] for i in indices]
            #self.slice_location_list=[slice_location_list[i] for i in indices]      
        
        if self.window is not None:
            self.window_low,self.window_high=get_window(window)
            ds=pydicom.dcmread(self.image_list[0])
            if hasattr(ds,'RescaleIntercept'):
                self.rescale_intercept=ds.RescaleIntercept
            else:
                self.rescale_intercept=0
        else:
            self.window_low=None
            self.window_high=None
    
        if self.DisplayAccessionNumber:
            ds=pydicom.dcmread(self.image_list[0])
            if hasattr(ds,'AccessionNumber'):
                self.AccessionNumber=ds.AccessionNumber
                
    def generateHTMLvideo(self):
        fig, ax1 = plt.subplots(1,1,figsize=(5,5))
        fig.set_tight_layout(True)
        ds=pydicom.dcmread(self.image_list[0])
        xw,yh=np.shape(ds.pixel_array)
        if self.window is not None:
            im1=ax1.imshow(ds.pixel_array, cmap = 'gray',vmin=self.window_low-self.rescale_intercept,vmax=self.window_high-self.rescale_intercept)
        else:
            im1=ax1.imshow(ds.pixel_array, cmap = 'gray')
        ax1.axes.get_xaxis().set_visible(False)
        ax1.axes.get_yaxis().set_visible(False)
        if self.AccessionNumber is not None:
            ax1.axes.text(xw*0.5,yh*0.9,self.AccessionNumber,fontsize=12.5,fontweight='bold',ha='center',color='forestgreen')

        def init():
            im1.set_data(ds.pixel_array)
            return(im1)

        def animate(i):
            ds=pydicom.dcmread(self.image_list[i])
            im1.set_data(ds.pixel_array)
            return (im1)

        self.anim = manimation.FuncAnimation(fig, animate,init_func=init, frames=len(self.image_list),interval=100)

        if self.save_file:
            writervideo = manimation.ImageMagickWriter(fps=self.fps)
            self.anim.save(self.output_location, writer=writervideo)
            # tested avi with ImageMagickWriter, can also do GIF with pillow
        #HTML(anim.to_js())        


class dcm_video_multiscan(dcm_visualization_video):
    def __init__(self,DCMFolderPathList,max_frames=40,window='lung',order_all_frames=True,outside_frame_list=None,output_location=os.path.abspath(os.path.join('..','Test.avi')),save_file=False,fps=30,maxVideoSize=100):
        super().__init__(output_location=output_location,save_file=save_file,fps=fps)
        self.FilePathList=DCMFolderPathList
        self.window=window
        self.max_frames=max_frames
        self.num_scans=len(DCMFolderPathList)        
        self.order_all_frames=order_all_frames
        #self.presort_images=presort_images
        matplotlib.rcParams['animation.embed_limit']=self.maxVideoSize
        #self.estimateVideoLength()
        self.estimateVideoSize()
        videoFactor=1/np.sqrt(self.estVideoSize/self.maxVideoSize)
        if videoFactor>=1:
            self.upscale_factor=1
        elif 0.3<videoFactor<1:
            self.upscale_factor=videoFactor
        else:
            self.upscale_factor=0.3
        
        
        if self.window is not None:
            self.window_low,self.window_high=get_window(window)
        else:
            self.window_low=None
            self.window_high=None
        
        if outside_frame_list is None:
            self.getFrameList()
        else:
            self.frame_idx_list=outside_frame_list['frame_idx_list']
            self.scan_idx_list=outside_frame_list['scan_idx_list']
            self.frame_loc_list=outside_frame_list['frame_loc_list']
            self.slice_loc_list=outside_frame_list['slice_loc_list']

        
    def loadSingleScan(self,index):
        CT_vis=dcm_video_singlescan(self.FilePathList[index],window=self.window)
        return CT_vis
    
    def estimateVideoLength(self):
        NumSlicesAllVideos=0
        for FilePath in self.FilePathList:
            dcm_series  = np.sort(os.listdir(FilePath))
            num_images = len(dcm_series)
            if self.max_frames=='all':
                NumSlicesAllVideos+=num_images
            else:
                NumSlicesAllVideos+=min(self.max_frames,num_images)
        self.estVideoFrames= NumSlicesAllVideos
    
    def estimateVideoSize(self):
        SizeAllVideos=0
        NumSlicesAllVideos=0
        for FilePath in self.FilePathList:
            dcm_series  = np.sort(os.listdir(FilePath))
            num_images = len(dcm_series)
            fname=FilePath + "/" + dcm_series[0]
            ds=pydicom.dcmread(fname)
            if self.max_frames=='all':
                SizeAllVideos+=num_images*ds.pixel_array.itemsize*ds.pixel_array.size
                NumSlicesAllVideos+=num_images
            else:
                SizeAllVideos+=min(self.max_frames,num_images)*ds.pixel_array.itemsize*ds.pixel_array.size
                NumSlicesAllVideos+=min(self.max_frames,num_images)
        self.estVideoSize= float(SizeAllVideos)/10**6
        self.estVideoFrames= NumSlicesAllVideos
    
    def getFrameList(self):
        scan_idx_list=[]
        frame_idx_list=[]
        frame_loc_list=[]
        slice_loc_list=[]
        
        for idx,FilePath in enumerate(self.FilePathList):
            if self.order_all_frames:
                image_list,slice_location_list_temp,_,_=return_dcm_list_sorted(FilePath)
                #dcm_series  = np.sort(os.listdir(FilePath))
                num_images = len(image_list)
                if (self.max_frames=='all') or (num_images<=self.max_frames):
                    frame_list_temp=list(np.round(np.linspace(0,num_images-1,num_images-1).astype('int')))
                    frame_idx_list.extend(frame_list_temp)
                    scan_idx_list.extend(list((np.ones_like(frame_list_temp)*idx).astype('int')))
                    frame_loc_list.extend(image_list)
                    slice_loc_list.extend(slice_location_list_temp)
                else:
                    frame_list_temp=list(np.round(np.linspace(0,num_images-1,self.max_frames)).astype('int'))
                    frame_idx_list.extend(frame_list_temp)
                    scan_idx_list.extend(list((np.ones_like(frame_list_temp)*idx).astype('int')))
                    frame_loc_list.extend(list(map(image_list.__getitem__, frame_list_temp)))
                    slice_loc_list.extend(list(map(slice_location_list_temp.__getitem__, frame_list_temp)))
            else:
                image_list,slice_location_list_temp=return_dcm_list_subset(FilePath,self.max_frames)
                frame_list_temp=list(np.round(np.linspace(0,len(image_list)-1,len(image_list))).astype('int'))
                frame_idx_list.extend(frame_list_temp)
                scan_idx_list.extend(list((np.ones_like(frame_list_temp)*idx).astype('int')))
                frame_loc_list.extend(image_list)
                slice_loc_list.extend(slice_location_list_temp)

        self.frame_idx_list=frame_idx_list
        self.scan_idx_list=scan_idx_list
        self.frame_loc_list=frame_loc_list
        self.slice_loc_list=slice_loc_list
    
    def exportFrameList(self):
        outside_frame_list={'frame_idx_list':copy.deepcopy(self.frame_idx_list),
                    'frame_loc_list':copy.deepcopy(self.frame_loc_list),
                    'scan_idx_list':copy.deepcopy(self.scan_idx_list),
                    'slice_loc_list':copy.deepcopy(self.slice_loc_list)}
        return outside_frame_list
    
    def generateHTMLvideo(self):        
        # initialize plot
        fig, ax1 = plt.subplots(1,1,figsize=(5,5))
        fig.set_tight_layout(True)
        ds=pydicom.dcmread(self.frame_loc_list[0])
        if self.upscale_factor==1:
            pix_array=ds.pixel_array
        else:
            pix_array=ndimage.zoom(ds.pixel_array,self.upscale_factor)
        xw,yh,rescale_intercept,AccessionNumber=getDisplayMetadata(ds,upscale_factor=self.upscale_factor)
        if self.window is not None:
            im1=ax1.imshow(pix_array, cmap = 'gray',vmin=self.window_low-rescale_intercept,vmax=self.window_high-rescale_intercept)
        else:
            im1=ax1.imshow(pix_array, cmap = 'gray')
        ax1.axes.get_xaxis().set_visible(False)
        ax1.axes.get_yaxis().set_visible(False)
        if AccessionNumber is not None:
            ax1.axes.text(xw*0.5,yh*0.9,AccessionNumber,fontsize=12.5,fontweight='bold',ha='center',color='forestgreen')
        
        ax1.axes.text(xw*0.8,yh*0.1,'Scan: '+str(1),fontsize=12.5,fontweight='bold',ha='center',color='firebrick')
        
        def animate(i):
            ds=pydicom.dcmread(self.frame_loc_list[i])
            if self.upscale_factor==1:
                pix_array=ds.pixel_array
            else:
                pix_array=ndimage.zoom(ds.pixel_array,self.upscale_factor)
            if i>0:
                if self.scan_idx_list[i]!=self.scan_idx_list[i-1]:
                    xw,yh,rescale_intercept,AccessionNumber=getDisplayMetadata(ds,self.upscale_factor)
                    ax1.patches = []
                    ax1.texts = []
                    ax1.axes.text(xw*0.5,yh*0.9,AccessionNumber,fontsize=12.5,fontweight='bold',ha='center',color='forestgreen')
                    ax1.axes.text(xw*0.8,yh*0.1,'Scan: '+str(self.scan_idx_list[i]+1),fontsize=12.5,fontweight='bold',ha='center',color='firebrick')
                    if self.window is not None:
                        im1.set_clim(self.window_low-rescale_intercept,self.window_high-rescale_intercept)
            im1.set_data(pix_array)
            

            return (im1)

        self.anim = manimation.FuncAnimation(fig,animate,frames=len(self.scan_idx_list),interval=100)

        if self.save_file:
            writervideo = manimation.ImageMagickWriter(fps=self.fps)
            self.anim.save(self.output_location, writer=writervideo)
            # tested avi with ImageMagickWriter, can also do GIF with pillow
            

def get_window(window_type):
    '''
    @brief: return window width based on https://radiopaedia.org/articles/windowing-ct?lang=us
    head and neck: brain W:80 L:40, subdural W:130-300 L:50-100, stroke W:8 L:32 or W:40 L:40, temporal bones W:2800 L:600 or W:4000 L:700, soft tissues: W:350–400 L:20–60    
    chest: lungs W:1500 L:-600, mediastinum W:350 L:50
    abdomen: soft tissues W:400 L:50, liver W:150 L:30
    spine: soft tissues W:250 L:50, bone W:1800 L:400
    '''
    
    window_low=None
    window_high=None
    window_dictionary={
        'brain':[80,40],
        'subdural':[200, 75],
        'stroke': [40,40],
        'temporal_bone':[2800,600],
        'head_soft_tissue': [375, 40],
        'lung':[1500,-600],
        'mediastinum':[350,50],
        'abd_soft_tissue':[400,50],
        'liver':[150,30],
        'spine_soft_tissue':[250,50],
        'spine_bone':[1800,400]
    }
    if window_type in list(window_dictionary.keys()):
        width=window_dictionary[window_type][0]
        level=window_dictionary[window_type][1]
        window_low=level-width/2
        window_high=level+width/2
    else:
        print('Warning: Window not found')
    return window_low,window_high

def getDisplayMetadata(ds,upscale_factor=1):
        xw,yh=np.shape(ds.pixel_array)
        if hasattr(ds,'RescaleIntercept'):
            rescale_intercept=ds.RescaleIntercept
        else:
            rescale_intercept=0
        if hasattr(ds,'AccessionNumber'):
            AccessionNumber=ds.AccessionNumber
        else:
            AccessionNumber=None
        xw=xw*upscale_factor
        yh=yh*upscale_factor
        return xw,yh,rescale_intercept,AccessionNumber

def quick_display_dcmseries(FilePath,window=None):
    image_list,slice_location_list,_,_ = return_dcm_list_sorted(FilePath)
    num_images = len(image_list)
        
    if num_images<26:
        num_rows=int(np.ceil(num_images/5))
        num_cols=5
        imlist_ds=image_list
    else:
        num_rows=5
        num_cols=5
        indices=list(np.round(np.linspace(0,num_images-1,25)).astype('int'))
        imlist_ds=list(map(image_list.__getitem__, indices))
    
    
    if window is not None:
        window_low,window_high=get_window(window)
        ds=pydicom.dcmread(imlist_ds[0])
        if hasattr(ds,'RescaleIntercept'):
            rescale_intercept=ds.RescaleIntercept
        else:
            rescale_intercept=0
    scale = 3
    plt.figure(figsize=(num_cols * scale, num_rows * scale) )
    i = 1
        
    for image in imlist_ds: 
        ds=pydicom.dcmread(image)
        plt.subplot(num_rows, num_cols, i)
        if window is not None:
            plt.imshow(ds.pixel_array, cmap = 'gray',vmin=window_low-rescale_intercept,vmax=window_high-rescale_intercept)
            # pydicom has a way to apply specified min / max from header, but in practice it doesn't work well
            # syntax as such: pydicom.pixel_data_handlers.apply_windowing(ds.pixel_array,ds)
            
            #plt.imshow(pydicom.pixel_data_handlers.apply_windowing(ds.pixel_array,ds), cmap = 'gray') 
        else:
            plt.imshow(ds.pixel_array, cmap = 'gray')
        plt.axis('off')
        i+=1

''' deprecated
def quick_display_dcmseries(FilePath):
    dcm_series  = np.sort(os.listdir(FilePath))
    num_images = len(dcm_series)
    image_list = []
    for dcm_file in dcm_series:
        image_list.append( FilePath + "/" + dcm_file )
        
    if num_images<26:
        num_rows=int(np.ceil(num_images/5))
        num_cols=5
        imlist_ds=image_list
    else:
        num_rows=5
        num_cols=5
        indices=list(np.round(np.linspace(0,num_images-1,25)).astype('int'))
        imlist_ds=list(map(image_list.__getitem__, indices))
    scale = 3
    plt.figure(figsize=(num_cols * scale, num_rows * scale) )
    i = 1
    for image in imlist_ds: 
        ds=pydicom.dcmread(image)
        plt.subplot(num_rows, num_cols, i)
        plt.imshow(ds.pixel_array, cmap = 'gray')
        plt.axis('off')
        i+=1
'''
        
def return_dcm_list_sorted(FilePath):
    '''
    @brief: take filepath, checks each dicom slice for slice location, and returns sorted slices
    @input: folder location of where DICOM files are located
    @output: image_list: list of fullfile images locations
            slice_location_list: list of slice locations
            image_list_skip: list of fullfile image locations that don't have slice locations
            skipcount: number of images in folder skipped
    
    To do: program in contingencies for cases where there are no slice locations at all
    '''
    #t0 = time.time()
    dcm_series  = np.sort(os.listdir(FilePath))
    #t1=time.time()
    #print('ListDir:',t1-t0)
    num_images = len(dcm_series)
    image_list = []
    slice_location_list = []
    image_list_skip=[]
    skipcount = 0
    
    for dcm_file in dcm_series:
        if '.dcm' in dcm_file:
            fname=FilePath + "/" + dcm_file
            ds=pydicom.dcmread(fname)
            if hasattr(ds, 'SliceLocation'):
                image_list.append(fname)
                slice_location_list.append(float(ds.SliceLocation))
            else:
                skipcount += 1
                image_list_skip.append(fname)
    #t2=time.time()
    if len(slice_location_list)==0 and len(image_list)==0:
        print('Warning, no slice locations found')
        image_list_outpt=image_list_skip
        slice_location_list=list(range(0,-len(image_list_outpt),-1))
    else:            
        slice_location_list,image_list = zip(*sorted(zip(slice_location_list,image_list),reverse=True))
        slice_location_list=list(slice_location_list)
        image_list_outpt=list(image_list)
    return image_list_outpt,slice_location_list,image_list_skip,skipcount


def return_dcm_list_subset(FilePath,max_frames='all',filename_only=False,fill_nan=False):
    '''
    @brief: take filepath, takes random number of slices
    @input: folder location of where DICOM files are located
    @output: image_list: list of fullfile images locations
            slice_location_list: list of slice locations
            image_list_skip: list of fullfile image locations that don't have slice locations
            skipcount: number of images in folder skipped
    
    To do: program in contingencies for cases where there are no slice locations at all
    '''
    image_list = []
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
    
    if filename_only:
        for dcm_file in image_list_temp:
            fname=FilePath + "/" + dcm_file
            ds=pydicom.dcmread(fname)
            if hasattr(ds, 'SliceLocation'):
                image_list.append(dcm_file)
                slice_location_list.append(float(ds.SliceLocation))
            elif fill_nan:
                image_list.append(dcm_file)
                slice_location_list.append(np.nan)
    else:
        for dcm_file in image_list_temp:
            fname=FilePath + "/" + dcm_file
            ds=pydicom.dcmread(fname)
            if hasattr(ds, 'SliceLocation'):
                image_list.append(fname)
                slice_location_list.append(float(ds.SliceLocation))
            elif fill_nan:
                image_list.append(dcm_file)
                slice_location_list.append(np.nan)

    if len(slice_location_list)>0:
        slice_location_list,image_list = zip(*sorted(zip(slice_location_list,image_list),reverse=True))
        slice_location_list=list(slice_location_list)

    return image_list,slice_location_list