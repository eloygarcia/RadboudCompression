import os
import sys
from subprocess import call
import pandas as pd
import shutil


# info = '/home/data/RadboudThicknes.csv'  ## path to the info file. Rows: ['Patient', 'BreastSide', 'Thickness']
BCT_basedir = '/home/data/'  ## Initial bCT segmentation folder at the docker container

# Aux definitions
ItkToCgal_path = '/home/Model-Extraction/ItkToCGAL' ## exec path at the docker container
##

# df = pd.read_csv(info)

#baseDir = '/home/eloygarcia/Escritorio/Phantom Validation/UCM-Results'
baseDir = '/home/data'

##
for root, subdirs, files in os.walk(BCT_basedir):
    patient = os.path.basename(root)
    side = 'R'
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

    # if side=='L' and os.path.exists(outputimagename):# and not os.path.exists(outputmeshname): ## rewriting biomechanical models
    if not os.path.exists(outputmeshname):
        side='R'
        text = [ItkToCgal_path,
                outputimagename,
                outputmeshname,
                '--facet_size', '2.',  # default = 2
                '--facet_distance', '0.5',  # default =1 ## Valor pequeño, más ajustado a la superficie, pero mayor problema para la optimizacion y el resolvedor
                '--cell_radius_edge_ratio', '2',  # default = 2
                '--cell_size', '2.',  # default = 2
                '--bc', 'CT',
                '--bc_thickness', '2',
                '--breast_side', 'R']
        call(text)

