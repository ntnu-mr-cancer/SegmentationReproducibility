# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 15:45:35 2018

@author: mohammed R. S. Sunoqrot
"""

import SimpleITK as sitk
import os
import glob
import sys

# Supress output screen
null = open(os.devnull, 'w')
sys.stdout = null
sys.stderr = null

# Get paths
with open(os.path.join(os.getcwd(),'paths.txt')) as f: 
    flines = f.readlines()
    inputdir = flines[0].strip()
    outputdir = flines[1].strip()

reader = sitk.ImageFileReader()
writer = sitk.ImageFileWriter()

files = []
os.chdir(inputdir)
for file in glob.glob("*.gz"):
    files.append(file)

for ii in  files:
    input_name = str(ii) 
    print(input_name)
    
    imagein = os.path.join(inputdir,input_name)
    imageout = os.path.join(outputdir,str(ii[:-7]+'_Segmentation.mhd'))
    
    reader.SetFileName(imagein)
    image = reader.Execute()
    writer.SetFileName(imageout)
    writer.Execute(image)