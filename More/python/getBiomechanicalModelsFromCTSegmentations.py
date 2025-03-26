import os
import sys
from subprocess import call
import pandas as pd
import shutil

info = './RadboudThicknes.csv' ## path to the info file. Rows: ['Patient', 'BreastSide', 'Thickness']
BCT_basedir = '/home/eloygarcia/Escritorio/Phantom Validation/phantoms' ## Initial bCT segmentation folder

ItkToCgal_path =  '/home/eloygarcia/Escritorio/Pruebas/Release Compression/Biomechanical Model/New Model Extraction/'

df = pd.read_csv(info)

for index, row in df.iterrows():
    ## Reading information
    patient = row['Patient']
    side = row['BreastSide']
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

    """
    std::cout << " --facet_angle : (default "<< argumentos.facet_angle << ")" << std::endl;
    std::cout << " --facet_size : (default "<< argumentos.facet_size << ")" << std::endl;
    std::cout << " --facet_distance : (default "<< argumentos.facet_distance << ")" << std::endl;
    std::cout << " --cell_radius_edge_ratio : (default "<< argumentos.cell_radius_edge_ratio << ")" << std::endl;
    std::cout << " --cell_size : (default "<< argumentos.cell_size << ")" << std::endl;
    """
    #if side=='L' and os.path.exists(outputimagename):# and not os.path.exists(outputmeshname): ## rewriting biomechanical models
    text = [ItkToCgal_path + 'ItkToCGAL',
            outputimagename,
            outputmeshname,
            '--facet_size', '2',  #default = 2
            '--facet_distance','0.2', #default =1
            '--cell_radius_edge_ratio', '1.5', # default = 2
            '--cell_size', '1.8', #default = 2
            '--bc', 'CT',
            '--bc_thickness', '2',
            '--breast_side', 'R']
    call(text)
    


