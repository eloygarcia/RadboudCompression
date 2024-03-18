import os
from glob import glob as glob
import numpy as np
import pandas as pd
from vtkmodules.vtkIOLegacy import (
    vtkUnstructuredGridReader
)

BCT_basedir = '/home/eloygarcia/Escritorio/Datasets/Breast Phantoms Radboud' ## Initial bCT segmentation folder

listOfModels = glob(os.path.join(BCT_basedir,'*/*_Segmented.vtk'))

names =[]
points = []
cells =[]

for model in listOfModels:
    if not model.endswith('.vtk.vtk'):
        reader = vtkUnstructuredGridReader()
        reader.SetFileName( model )
        reader.Update()

        names.append( os.path.basename(model) )
        points.append( reader.GetOutput().GetNumberOfPoints() )
        cells.append( reader.GetOutput().GetNumberOfCells() )

info = pd.DataFrame( {'Name': names,
                      'Points':points,
                      'Cells':cells}  )
#info.to_csv('./primary-info.csv')
info.to_csv('./newTray-info.csv')