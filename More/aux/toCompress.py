#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 16:42:05 2021

@author: eloy
"""
from subprocess import call
import computeDistances

def compressionFunction(patient, thickness, image_path, mesh_path, output_path):
    ### Create NiftySim Model
    text = ['/home/eloy/Escritorio/Oliver/Radboud/Koen/release/Compression/WriteNiftyModel',
        mesh_path, image_path, output_path + '/'+ patient + '-niftySimModel.xml',
        thickness]
    call(text)
    
    ### Exec Niftysim
    text = ['/home/eloy/Escritorio/Oliver/Libraries/NiftySim-250-release/source/niftysim',
            '-x', 
            output_path + '/'+ patient + '-niftySimModel.xml',
            '-v', '-t', '-sport',
            '-export-mesh',
            output_path+ '/'+ patient + '-outputMesh.vtk']
    call(text)
    
    computeDistances( )