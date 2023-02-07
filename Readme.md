# Breast Phantom Generation 

This document summarizes the whole pipeline, from CT images to compressed breast phantoms, by means of a FEM simulation.

The algorithm is divided into 3 main parts:

-  Biomechanical model construction
- FEM simulation - i.e. breast compression
- Image (i.e. phantom) - reconstruction.

### Dependencies

This software was written using c++ and, mainly, the [itk](https://itk.org/) and [vtk](https://vtk.org/) libraries.

The biomechanical model extraction is performed using the library [iso2mesh](https://github.com/fangq/iso2mesh) and, therefore, we need to use either Matlab$^{TM}$ or Octave$^{TM}$. We are trying a code based on the [cgal library](https://www.cgal.org/) to directly transform itk images into 3D meshes.

### Current situation

The current code was not directly written to create bCT phantoms but as a part of a multimodal registration pipeline. Therefore, it may fail in several parts.

Currently, I'm Checking and sharing the code.