#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 16:44:09 2021

@author: eloy
"""
#!/usr/bin/env python
import numpy as np
import argparse

# noinspection PyUnresolvedReferences
# import vtkmodules.vtkInteractionStyle

# noinspection PyUnresolvedReferences
#import vtkmodules.vtkRenderingOpenGL2
# from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkIOGeometry import (
    # vtkSTLReader,
    vtkOBJReader
    )
from vtkmodules.vtkFiltersCore import vtkImplicitPolyDataDistance
# from vtkmodules.vtkRenderingCore import (
#     vtkActor,
#     vtkDataSetMapper,
#     vtkRenderWindow,
#     vtkRenderWindowInteractor,
#     vtkRenderer
# )
from vtkmodules.vtkCommonCore import (
    #vtkPoints,
    vtkFloatArray
    )
from vtk import (
    vtkUnstructuredGridReader, 
    vtkTransform,
    vtkTransformFilter
    )
from vtkmodules.vtkFiltersGeometry import (
    #vtkGeometryFilter,
    vtkDataSetSurfaceFilter 
    )

def get_program_parameters():
    description = 'Read a .stl file.'
    epilogue = ''''''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', help='42400-IDGH.stl')
    args = parser.parse_args()
    return args.filename

def close_window(iren):
    render_window = iren.GetRenderWindow()
    render_window.Finalize()
    iren.TerminateApp()
    del render_window, iren

def computeDistances(obj_filename, vtk_filename): 
    #colors = vtkNamedColors()
    reader = vtkOBJReader()
    reader.SetFileName(obj_filename)
    reader.Update()
    # points = reader.GetOutput().GetPoints()
    
    reader_mesh = vtkUnstructuredGridReader()
    reader_mesh.SetFileName( vtk_filename)
    reader_mesh.Update()
    
    # surface_filter = vtkGeometryFilter()
    surface_filter = vtkDataSetSurfaceFilter()
    surface_filter.SetInputConnection( reader_mesh.GetOutputPort() )
    surface_filter.Update()
    
    bb1 = np.array(surface_filter.GetOutput().GetBounds())
    
    tr = vtkTransform()
    # tr.Translate(bb1[0]-bb0[0], bb1[2]-bb0[1], bb1[1]-bb0[2])
    tr.RotateX(90)
    ## Checkear la dirección a ver si está bien!!
    
    tp = vtkTransformFilter()
    tp.SetInputConnection(reader.GetOutputPort())
    tp.SetTransform(tr)
    tp.Update()
    
    #points2 = tp.GetOutput()
    bb0 = np.array(tp.GetOutput().GetBounds())
    
    tr2 = vtkTransform()
    tr2.Translate(bb1[0]-bb0[0], 
                  (bb1[3]+bb1[2])/2 - (bb0[3]+bb0[2])/2,
                  bb1[5]-bb0[5])
    # tr.RotateX(-90)
    
    tp2 = vtkTransformFilter()
    tp2.SetInputConnection(tp.GetOutputPort())
    tp2.SetTransform(tr2)
    tp2.Update()
    
    points3 =tp2.GetOutput()
    
    implicitPolyDataDistance = vtkImplicitPolyDataDistance()
    implicitPolyDataDistance.SetInput(surface_filter.GetOutput())
    
    # Add distances to each point
    signedDistances = vtkFloatArray()
    signedDistances.SetNumberOfComponents(1)
    signedDistances.SetName('SignedDistances')
    
    # Evaluate the signed distance function at all of the grid points
    for pointId in range(points3.GetNumberOfPoints()):
        p = points3.GetPoint(pointId)
        signedDistance = implicitPolyDataDistance.EvaluateFunction(p)
        signedDistances.InsertNextValue(signedDistance)
        
    distances = np.array(signedDistances)
    
    dist=np.mean(np.abs(distances))
    return dist