# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 15:55:07 2018

@author: mattjise
Minot modifications by Mohammed R. S. Sunoqrot Oct 2019
"""

import radiomics
import SimpleITK as sitk
import numpy as np
import json
import os
import glob

def getSettings(image_array_in,mask_array_in,nr_bins):
    intensity_range = np.max(image_array_in[mask_array_in == 1])-np.min(image_array_in[mask_array_in == 1])
    settings = {} 
    settings['binWidth'] = intensity_range/64
    settings['correctMask'] = True
    return settings

# Get paths
with open(os.path.join(os.getcwd(),'paths.txt')) as f: 
    flines = f.readlines()
    image_dir = flines[0].strip()
    mask_dir = flines[1].strip()
    results_dir = flines[2].strip()
    # Get the class
    region_class = flines[3].strip()

files = []
os.chdir(mask_dir)
for file in glob.glob("*.mhd"):
    files.append(file)

for patient_nr in files:
        image_name = str(patient_nr[0:18])+'_Normalized.mhd'
        mask_name = str(patient_nr)
        print(image_name, mask_name)
    
        # read in data
        image = sitk.ReadImage(os.path.join(image_dir,image_name))
        mask = sitk.ReadImage(os.path.join(mask_dir,mask_name))
        mask = sitk.Cast( mask, sitk.sitkUInt8)

        mask.SetDirection(image.GetDirection())
        mask.SetOrigin(image.GetOrigin())
        mask_array = sitk.GetArrayFromImage(mask)
        image_array = sitk.GetArrayFromImage(image)

         # whole prostate - shape 3D
        feature_class = 'shape'
        settings = getSettings(image_array,mask_array,64)       
        extractor = radiomics.featureextractor.RadiomicsFeatureExtractor(**settings)
        extractor.disableAllFeatures()
        extractor.enableFeatureClassByName(feature_class)
        featureVector = extractor.execute(image,mask)
        for key in featureVector.keys():
            if type(featureVector[key]) == type(np.array(1)):
                featureVector[key] = float(featureVector[key])
        results_name = str(patient_nr[:-17])+'_'+region_class+'_'+feature_class+'.json'
        with open(os.path.join(results_dir,results_name), 'w') as f:
            f.write(json.dumps(featureVector))
 