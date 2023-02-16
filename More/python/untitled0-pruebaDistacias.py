#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 11:34:45 2021

@author: eloy
"""


"""
CREO QUE LO TENGO!!!
"""

#!/usr/bin/env python
import numpy as np
import argparse

# noinspection PyUnresolvedReferences
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
#import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkIOGeometry import (
    # vtkSTLReader,
    vtkOBJReader
    )
from vtkmodules.vtkFiltersCore import vtkImplicitPolyDataDistance
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkDataSetMapper,
    vtkPolyDataMapper,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer,
    vtkColorTransferFunction
)
from vtkmodules.vtkCommonCore import (
    #vtkPoints,
    vtkFloatArray,
    vtkIntArray
    )
from vtkmodules.vtkCommonColor import (
    vtkColorSeries,
    vtkNamedColors
)
from vtkmodules.vtkCommonDataModel import (
    # vtkIterativeClosestPointTransform,
    vtkPolyData
)
from vtk import (
    vtkUnstructuredGridReader, 
    vtkTransform,
    vtkTransformFilter,
    vtkVertexGlyphFilter,
    vtkReflectionFilter,
    vtkSimplePointsWriter
    )
from vtkmodules.vtkFiltersGeometry import (
    #vtkGeometryFilter,
    vtkDataSetSurfaceFilter 
    )
from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor

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

colors = vtkNamedColors()

BCT_basedir = '/home/eloy/Escritorio/Oliver/Radboud/to_compressed_OLIVER/'
BCT_references = ['BCT_16006_002', 'BCT_16006_004', 'BCT_16006_016',
                  'BCT_16006_021', 'BCT_16006_057', 'BCT_16006_097',
                  'BCT_16006_113', 'BCT_16006_124']
BCT_thickness = [73, 67, 68, 68, 76, 64, 35, 67]

obj_dir = '/home/eloy/Escritorio/Oliver/Radboud/Koen/transfer_1198365_files_56e8c083/'
obj_files = ['BCT_patient01', 'BCT_patient03', 'BCT_patient11',
             'BCT_patient15', 'BCT_patient42', 'BCT_patient61',
             'BCT_patient72', 'BCT_patient80']

i=0
reflection= False
maximum = False
tipi =''

filename = obj_dir + obj_files[i] + tipi+ '.obj'
mesh_filename = BCT_basedir + BCT_references[i] +'/'+ str( BCT_thickness[i]) + \
    '/MeshPhantom-'+BCT_references[i]+'-'+str( BCT_thickness[i]+1) +'-0.vtk'

reader = vtkOBJReader()
reader.SetFileName(filename)
reader.Update()
points = reader.GetOutput().GetPoints()

# bb0 = points.GetBounds()
# mapper = vtkDataSetMapper()
# mapper.SetInputConnection(reader.GetOutputPort())
# mapper.SetScalarRange(4, 9);
# mapper.Update()

reader_mesh = vtkUnstructuredGridReader()
reader_mesh.SetFileName( mesh_filename)
reader_mesh.Update()

mapper_mesh = vtkDataSetMapper()
mapper_mesh.SetInputConnection(reader_mesh.GetOutputPort())
mapper_mesh.SetScalarRange(4, 9)
mapper_mesh.Update()

# surface_filter = vtkGeometryFilter()
surface_filter = vtkDataSetSurfaceFilter()
surface_filter.SetInputConnection( reader_mesh.GetOutputPort() )
surface_filter.Update()

bb1 = np.array(surface_filter.GetOutput().GetBounds())

# Rotate
tr = vtkTransform()
# tr.Translate(bb1[0]-bb0[0], bb1[2]-bb0[1], bb1[1]-bb0[2])
tr.RotateX(90)
## Checkear la dirección a ver si está bien!!

tp = vtkTransformFilter()
tp.SetInputConnection(reader.GetOutputPort())
tp.SetTransform(tr)
tp.Update()

if reflection:
    # Flip
    reflection = vtkReflectionFilter()
    reflection.SetInputConnection( tp.GetOutputPort())
    reflection.SetPlane(2)
    reflection.CopyInputOff()
    reflection.Update()

"""
Called when mirrorPlane != None
The reflection plane is labeled as follows: From the vtk documentation: 
ReflectionPlane {
USE_X_MIN = 0, USE_Y_MIN = 1, USE_Z_MIN = 2, USE_X_MAX = 3,
USE_Y_MAX = 4, USE_Z_MAX = 5, USE_X = 6, USE_Y = 7,
USE_Z = 8
}
"""

#points2 = tp.GetOutput()
if reflection: bb0 = np.array(reflection.GetOutput().GetBounds())
else: bb0 = np.array(tp.GetOutput().GetBounds())

tr2 = vtkTransform()
if maximum:
    print('maximum')
    tr2.Translate(bb1[0]-bb0[0], 
              bb1[2]-bb0[2],
              bb1[5]-bb0[5])
else:
    print('average')
    tr2.Translate(bb1[0]-bb0[0], 
              (bb1[3]+bb1[2])/2 - (bb0[3]+bb0[2])/2,
              bb1[4]-bb0[4])
# tr.RotateX(-90)

tp2 = vtkTransformFilter()
if reflection: tp2.SetInputConnection(reflection.GetOutputPort())
else: tp2.SetInputConnection(tp.GetOutputPort())
tp2.SetTransform(tr2)
tp2.Update()

points3 =tp2.GetOutput()

# mapper = vtkDataSetMapper()
# mapper.SetInputConnection(reader_mesh.GetOutputPort())
# mapper.SetScalarRange(4, 9)
# mapper.Update()

implicitPolyDataDistance = vtkImplicitPolyDataDistance()
implicitPolyDataDistance.SetInput(surface_filter.GetOutput())

# Add distances to each point
signedDistances = vtkFloatArray()
signedDistances.SetNumberOfComponents(1)
signedDistances.SetName('SignedDistances')

distancesInt = vtkIntArray()
distancesInt.SetNumberOfComponents(1)
distancesInt.SetName('distancesInte')

# Evaluate the signed distance function at all of the grid points
for pointId in range(points3.GetNumberOfPoints()):
    p = points3.GetPoint(pointId)
    signedDistance = implicitPolyDataDistance.EvaluateFunction(p)
    signedDistances.InsertNextValue(signedDistance)
    
    sdist = np.array(signedDistance, dtype='int')
    distancesInt.InsertNextValue(sdist)
    
distances = np.array(signedDistances)

dist=np.mean(np.abs(distances))

# save points
writer = vtkSimplePointsWriter()
writer.SetFileName(obj_dir + obj_files[i] + tipi +'.xyz');
writer.SetInputData(points3)
writer.Write();

# return dist

## Checkear y visualizar!!!
    
# ### check a partir de aquí.
polyData = vtkPolyData()
polyData.SetPoints(points3.GetPoints())
polyData.GetPointData().SetScalars(signedDistances)

vertexGlyphFilter = vtkVertexGlyphFilter()
vertexGlyphFilter.SetInputData(polyData)
vertexGlyphFilter.Update()

signedDistanceMapper = vtkPolyDataMapper()
signedDistanceMapper.SetInputConnection(vertexGlyphFilter.GetOutputPort())
signedDistanceMapper.ScalarVisibilityOn()
   

#####
color_map_idx = 1
# Build a lookup table
color_series = vtkColorSeries()
color_series.SetColorScheme(color_map_idx)
print(f'Using color scheme #: {color_series.GetColorScheme()}, {color_series.GetColorSchemeName()}')

lut = vtkColorTransferFunction()
lut.SetColorSpaceToHSV()

# scalar_range = polyData.GetPointData().GetScalars().GetRange()
scalar_range = distancesInt.GetRange()

# Use a color series to create a transfer function
for i in range(0, color_series.GetNumberOfColors()):
    color = color_series.GetColor(i)
    double_color = list(map(lambda x: x / 255.0, color))
    t = scalar_range[0] + (scalar_range[1] - scalar_range[0]) / (color_series.GetNumberOfColors() - 1) * i
    lut.AddRGBPoint(t, double_color[0], double_color[1], double_color[2])

colors = vtkNamedColors()

signedDistanceMapper.SetScalarModeToUsePointFieldData()
# signedDistanceMapper.SelectColorArray(distancesInt)
signedDistanceMapper.SetScalarModeToUsePointData() 
signedDistanceMapper.SetScalarRange(scalar_range)
signedDistanceMapper.SetLookupTable(lut)

# ###
actor = vtkActor()
# actor.SetMapper(mapper)
actor.SetMapper(signedDistanceMapper)
actor.GetProperty().SetDiffuse(0.8)
actor.GetProperty().SetDiffuseColor(colors.GetColor3d('LightSteelBlue'))
actor.GetProperty().SetSpecular(0.3)
actor.GetProperty().SetSpecularPower(60.0)
actor.GetProperty().SetPointSize(5);

actor_mesh = vtkActor()
# actor.SetMapper(mapper)
actor_mesh.SetMapper(mapper_mesh)
actor_mesh.GetProperty().SetDiffuse(0.8)
actor_mesh.GetProperty().SetDiffuseColor(colors.GetColor3d('LightSteelBlue'))
actor_mesh.GetProperty().SetSpecular(0.3)
actor_mesh.GetProperty().SetSpecularPower(60.0)
actor_mesh.GetProperty().SetOpacity(0.5)
# actor_mesh.GetProperty().SetPointSize(10);

window_width = 800
window_height = 800
# Create a scalar bar
scalar_bar = vtkScalarBarActor()
scalar_bar.SetLookupTable(signedDistanceMapper.GetLookupTable())
# scalar_bar.SetTitle(signedDistances.replace('_', '\n'))
scalar_bar.UnconstrainedFontSizeOn()
scalar_bar.SetNumberOfLabels(5)
scalar_bar.SetMaximumWidthInPixels(window_width // 8)
scalar_bar.SetMaximumHeightInPixels(window_height // 3)

# Create a rendering window and renderer
ren = vtkRenderer()
renWin = vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(window_width, window_height)
#renWin.SetWindowName('ReadSTL')

# Create a renderwindowinteractor
iren = vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# Assign actor to the renderer
ren.AddActor(actor)
ren.AddActor(actor_mesh)
ren.AddActor2D(scalar_bar)
ren.SetBackground(colors.GetColor3d('DarkOliveGreen'))

# Enable user interface interactor
iren.Initialize()
renWin.Render()
iren.Start()

close_window(iren)
del renWin, iren
