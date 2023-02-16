#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 15:55:52 2021

@author: eloy
"""
import os
import sys
from glob import glob 
import pandas as pd

from aux import toCompress

info = pd.read_csv('/home/eloy/Escritorio/Oliver/Radboud/newModels/Radboud_phantoms.csv')

BCT_basedir = '/home/eloy/Escritorio/Oliver/Radboud/to_compressed_OLIVER/'
BCT_references = ['BCT_16006_002', 'BCT_16006_004', 'BCT_16006_016',
                  'BCT_16006_021', 'BCT_16006_057', 'BCT_16006_097',
                  'BCT_16006_113', 'BCT_16006_124']

obj_dir = '/home/eloy/Escritorio/Oliver/Radboud/Koen/transfer_1198365_files_56e8c083/'
obj_files = ['BCT_patient01', 'BCT_patient03', 'BCT_patient11',
             'BCT_patient15', 'BCT_patient42', 'BCT_patient61',
             'BCT_patient72', 'BCT_patient80']




for i in range(len(obj_files)):
    temp_info = info[ info['Reference']==BCT_references[i] ]
    thickness = temp_info['Thickness'].values.astype(str)[0]
    
    image_path = BCT_basedir + BCT_references[i] + '/' + BCT_references[i] +'_Segmented.nrrd'
    mesh_path = BCT_basedir + BCT_references[i] + '/' + BCT_references[i] +'_Segmented.vtk'
    output_path = BCT_basedir + BCT_references[i] + '/' + thickness 
    
    toCompress.compressionFunction( BCT_references[i], thickness,
                        image_path, mesh_path, output_path)
    
    text = ['/home/eloy/Escritorio/Oliver/Radboud/Koen/release/ReconstructImage/ReconstructImage',
        output_path + '/'+ BCT_references[i] + '-outputMesh.vtk',
        image_path, 
        output_path + '/'+ BCT_references[i] + '-compressedMesh.vtk',
        output_path + '/'+ BCT_references[i] + '-phantom.nrrd']
    call(text)
