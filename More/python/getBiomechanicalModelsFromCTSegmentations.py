import os
import sys
from subprocess import call
import pandas as pd
import shutil

info = '/home/eloygarcia/Escritorio/Datasets/RadboudThicknes.csv' ## path to the info file. Rows: ['Patient', 'Thickness']
BCT_basedir = '/home/eloygarcia/Escritorio/Datasets/Breast Phantoms Radboud' ## Initial bCT segmentation folder

ItkToCgal_path =  '/home/eloygarcia/Escritorio/Pruebas/Release Compression/Biomechanical Model/New Model Extraction/'

df = pd.read_csv(info)

for index, row in df.iterrows():
    ## Reading information
    patient = row['Patient ']
    thickness = row['Thickness']
    print(patient)

    ## Looking for the corresponding image:
    imagename = os.path.join(BCT_basedir, patient + '_Segmented.nrrd')
    outputdir = os.path.join(BCT_basedir, patient)
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    outputimagename = os.path.join(outputdir, patient + '_Segmented.nrrd')
    outputmeshname = os.path.join(outputdir, patient + '_Segmented.vtk')

    if os.path.exists(imagename) and not os.path.exists(outputimagename):
        shutil.copyfile(imagename, outputimagename )

    if os.path.exists(outputimagename) and not os.path.exists(outputmeshname):
        text = [ItkToCgal_path + 'ItkToCGAL',
                outputimagename,
                outputmeshname,
                '--bc', 'CT',
                '--bc_thickness', '2']
        call(text)
    


