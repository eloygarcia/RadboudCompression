import os
import sys
from subprocess import call
import sys
sys.path.append('./aux/')
import toCompress
import numpy as np
np.bool = np.bool_

import pandas as pd
import shutil

from vtkmodules.vtkIOGeometry import (
    vtkOBJReader
    )
from vtk.util.numpy_support import vtk_to_numpy
from vtk import (
    vtkUnstructuredGridReader
    )
from vtkmodules.vtkFiltersGeometry import (
    #vtkGeometryFilter,
    vtkDataSetSurfaceFilter
    )
from vtkmodules.vtkIOLegacy import (
    vtkPolyDataReader,
    vtkPolyDataWriter
    )

from toCompress import compressionFunction
from compDistances import computeDistances

np.bool = np.bool_

info = '/home/eloygarcia/Escritorio/Datasets/RadboudThicknes.csv' ## path to the info file. Rows: ['Patient', 'Thickness']
BCT_basedir = '/home/eloygarcia/Escritorio/Datasets/Breast Phantoms Radboud' ## Initial bCT segmentation folder
PCA_basedir = '/home/eloygarcia/Escritorio/Pruebas/abreast-generator/Phantoms1'

df = pd.read_csv(info)

for index, row in df.iterrows():
    ## Reading information
    patient = row['Patient ']
    thickness = row['Thickness']
    print(patient)
    print(thickness)
    print(' ')

    patientdir = os.path.join(BCT_basedir, patient)
    meshname = os.path.join(patientdir, patient + '_Segmented.vtk')
    pcafilename = os.path.join(PCA_basedir, 'Thickness-'+str(thickness)+'.vtk')

    if os.path.exists(pcafilename) and os.path.exists(meshname):
        ## reading pca
        #pcareader = vtkPolyDataReader()
        #pcareader.SetFileName(pcafilename)
        #pcareader.Update()

        #print(pcareader.GetOutput().GetPoints().GetBounds())

        ## No necesito leer la malla inicial, solo el nombre para pasarselo al generador de documentos .xml de nifitysim

        ## performing compression
        nifysimmodelname = os.path.join(patientdir, patient +'-niftysim.xml')
        if not os.path.exists(nifysimmodelname):
            compressionFunction(patient, thickness, meshname, nifysimmodelname)


        ## reading vtk compressed mesh
        vtkfilename = os.path.join(patientdir, patient+'-compressedMesh.vtk')
        #reader_mesh = vtkUnstructuredGridReader()
        #reader_mesh.SetFileName(vtkfilename)
        #reader_mesh.Update()

        ## Optimization
        metric = computeDistances(pcafilename, vtkfilename)
        print(metric)

        n_iters = 10
        ## initial params:
        gravity = 150 ## Check
        offset = 20 ## mm
        ## optimization params:
        grav_step = gravity/n_iters
        offset_step = offset/n_iters


        for i in range(n_iters):
            temp_gravity = gravity-(i*grav_step)

            ## compression function
            compressionFunction(patient, thickness, meshname, nifysimmodelname, temp_gravity)

            ## computing distances
            vtkfilename = os.path.join(patientdir, patient +'-'+str(temp_gravity)+'-compressedMesh.vtk')
            temp_metric = computeDistances(pcafilename, vtkfilename)
            print(temp_metric)

            if temp_metric < metric:
                metric = temp_metric
    break

"""
point_array = np.array( pcareader.GetOutput().GetPoints().GetData() )
import matplotlib.pyplot as plt

plt.plot(point_array[0], point_array[1])
plt.show()
"""