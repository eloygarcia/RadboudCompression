## Hi Eloy!
# Fist of all:
# - check swap memory
# - Release swap memory
# - run from terminal
#
#
# TO DO:
# - Write a main file to directly use the terminal, including input and output filenames

import SimpleITK as sitk
# import pylab
import imageio.v3 as iio
import numpy as np

####
inputfilename = '/home/eloygarcia/Escritorio/UCM/models/hetero/hetero_50um.tif'
outputfilename = '/home/eloygarcia/Escritorio/UCM/models/hetero/hetero_50um.nrrd'

#### SCATTERED
# inputfilename = '/home/eloygarcia/Escritorio/UCM/models/scattered/scattered_50um.tif'
# outputfilename = '/home/eloygarcia/Escritorio/UCM/models/scattered/scattered_50um.nrrd'

sitk.WriteImage( sitk.GetImageFromArray(iio.imread(inputfilename)), outputfilename, useCompression=True)
