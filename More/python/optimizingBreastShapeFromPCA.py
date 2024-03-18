import os
import pandas as pd
import numpy as np
np.bool = np.bool_
os.chdir('/home/eloygarcia/Escritorio/Pruebas/RadboudCompression/More/python')
import sys
sys.path.append('./aux/')

from toCompress import compressionFunction
from compDistances import computeDistances

from scipy import optimize

np.bool = np.bool_

info = '/home/eloygarcia/Escritorio/Datasets/RadboudThicknes.csv' ## path to the info file. Rows: ['Patient', 'Thickness']
BCT_basedir = '/home/eloygarcia/Escritorio/Datasets/Breast Phantoms Radboud' ## Initial bCT segmentation folder
PCA_basedir = '/home/eloygarcia/Escritorio/Pruebas/abreast-generator/Results'

df = pd.read_csv(info)

global meshname
global pcafilname
global nifysimmodelname
global patient
global thickness

class MyTakeStep:
   def __init__(self, stepsize=10, stepratio=[1]): ## incluir un ratio de valores
       self.stepsize = stepsize
       self.stepratio = stepratio
       self.rng = np.random.default_rng()
   def __call__(self, x):
       s = self.stepsize
       for i in range(len(self.stepratio)):
           x[i] += self.rng.uniform(-s*self.stepratio[i], s*self.stepratio[i])
       #x[0] += self.rng.uniform(-4*s, 4*s) ## al incluir el ratio, puedo hacer que estas lineas se resuelvan con un solo bucle
       #x[1] += self.rng.uniform(-s, s)
       return x

# def optimizationFunction( thickness, gravity, offset):
def optimizationFunction( variable):
    gravity = variable[0]
    offset = variable[1]

    ## compression function
    idx = compressionFunction(patient, thickness, meshname, nifysimmodelname, gravity, offset)
    if idx==0:
        return ValueError('No files!')


    ## computing distances
    vtkfilename = os.path.join(patientdir,
                               patient + '-' + str(gravity) + '-' + str(offset) + '-compressedMesh.vtk')

    metric = computeDistances(pcafilename, vtkfilename)

    print(metric)
    return(metric)

for index, row in df.iterrows():
    ## Reading information
    patient = row['Patient ']
    thickness = row['Thickness']
    print(patient)
    print(thickness)
    print(' ')

    patientdir = os.path.join(BCT_basedir, patient)
    meshname = os.path.join(patientdir, patient + '_Segmented.vtk')
    pcafilename = os.path.join(PCA_basedir, patient+'.vtk')

    if os.path.exists(pcafilename) and os.path.exists(meshname):
        ## reading pca
        #pcareader = vtkPolyDataReader()
        #pcareader.SetFileName(pcafilename)
        #pcareader.Update()

        #print(pcareader.GetOutput().GetPoints().GetBounds())

        ## No necesito leer la malla inicial, solo el nombre para pasarselo al generador de documentos .xml de nifitysim

        ## performing compression
        nifysimmodelname = os.path.join(patientdir, patient +'-niftysim.xml')
        #if not os.path.exists(nifysimmodelname):
        #    compressionFunction(patient, thickness, meshname, nifysimmodelname)


        ## reading vtk compressed mesh
        vtkfilename = os.path.join(patientdir, patient+'-compressedMesh.vtk')
        #reader_mesh = vtkUnstructuredGridReader()
        #reader_mesh.SetFileName(vtkfilename)
        #reader_mesh.Update()

        ## Optimization
        # initial params:
        gravity = 100 ## Check
        offset = 20 ## mm
        x0 = [gravity, offset]
        ## bounds
        bounds =((0,(2*gravity)), (0,(2*offset)))
        ratio = x0/np.min(x0)

        mytakestep = MyTakeStep(stepsize=np.min(x0), stepratio=ratio)

        res = optimize.basinhopping(optimizationFunction, x0=x0 , T=100,
                                    niter=10, disp=True, take_step=mytakestep)#, bounds=bounds)
        print(f"Minima found at {res[0][0]} with function value {res[1]}")
        # plot_fitness(f, np.array(lst))
        # plt.vlines(x=res[0][0], lw=3, ymin=300, ymax=-400, color='k', linestyle='--')


        """
        #for i in range(n_iters+1):
        while stop==False:
            for j in range(2):
                if j==0:
                    temp_gravity = gravity - grav_step
                else:
                    temp_gravity = gravity + grav_step

                temp_offset = offset #-(i*offset_step)

                ## compression function
                compressionFunction(patient, thickness, meshname, nifysimmodelname, temp_gravity, temp_offset)

                ## computing distances
                vtkfilename = os.path.join(patientdir, patient +'-'+str(temp_gravity)+'-'+str(temp_offset)+'-compressedMesh.vtk')
                temp_metric[j] = computeDistances(pcafilename, vtkfilename)
                print(temp_metric[j])

                if temp_metric[j] < metric:
                    metric = temp_metric[j]
                    gravity = temp_gravity
                if
        """
    break

"""
point_array = np.array( pcareader.GetOutput().GetPoints().GetData() )
import matplotlib.pyplot as plt

plt.plot(point_array[0], point_array[1])
plt.show()
"""