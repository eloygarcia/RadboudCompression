import os
import sys
from subprocess import call
import pandas as pd
import shutil


info = '/home/data/Radboud-finalThickness.csv'  ## path to the info file. Rows: ['Patient', 'BreastSide', 'Thickness']
BCT_basedir = '/home/data/'  ## Initial bCT segmentation folder at the docker container

# Aux definitions
ItkToCgal_path = '/home/Model-Extraction/ItkToCGAL' ## exec path at the docker container
##

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
        shutil.copyfile(imagename, outputimagename)

    """
    std::cout << " --facet_angle : (default "<< argumentos.facet_angle << ")" << std::endl;
    std::cout << " --facet_size : (default "<< argumentos.facet_size << ")" << std::endl;
    std::cout << " --facet_distance : (default "<< argumentos.facet_distance << ")" << std::endl;
    std::cout << " --cell_radius_edge_ratio : (default "<< argumentos.cell_radius_edge_ratio << ")" << std::endl;
    std::cout << " --cell_size : (default "<< argumentos.cell_size << ")" << std::endl;
    """
    # if side=='L' and os.path.exists(outputimagename):# and not os.path.exists(outputmeshname): ## rewriting biomechanical models
    if patient==('BCT_16006_056') :#  or patient==('BCT_16006_076'): # or patient==('BCT_16006_084') or patient==('BCT_16006_028'):
        text = [ItkToCgal_path,
                outputimagename,
                outputmeshname,
                '--facet_size', '3',  # default = 2
                '--facet_distance', '2',  # default =1 ## Valor pequeño, más ajustado a la superficie, pero mayor problema para la optimizacion y el resolvedor
                '--cell_radius_edge_ratio', '2',  # default = 2
                '--cell_size', '3',  # default = 2
                '--bc', 'CT',
                '--bc_thickness', '2',
                '--breast_side', 'L']
        call(text)


