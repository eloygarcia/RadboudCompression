import os
from subprocess import call

#baseDir = '/home/eloygarcia/Escritorio/Phantom Validation/UCM-Results'
baseDir = '/home/data'

## aux definitions
reconstructionPath = '/home/ReconstructImage/release/ReconstructImage'
##

for root, subdirs, files in os.walk(baseDir):
    patient = os.path.basename(root)

    imagename = os.path.join(root, patient +'_Segmented.nrrd')
    initialmesh = os.path.join(root, patient +'_Segmented.vtk')
    # compressedmesh = os.path.join(root, patient +'-compressedMesh.vtk')
    compressedmesh = initialmesh

    outputfile =  os.path.join(root, 'Phantom-'+patient+'.nrrd')
    print(compressedmesh)
    print(os.path.exists(compressedmesh))
    if os.path.exists(imagename) and os.path.exists(initialmesh) and os.path.exists(compressedmesh) and not os.path.exists(outputfile):
        text = [reconstructionPath,
                initialmesh, compressedmesh, imagename, outputfile]
        print('Beggining Image Reconstruction')

        value = call(text)
        print(value)
        #if not value == 0:
        #    return ValueError