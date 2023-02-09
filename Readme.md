# Breast Phantom Generation 

This document summarizes the whole pipeline, from CT images to compressed breast phantoms, by means of a FEM simulation.

The algorithm is divided into 3 main parts:

-  Biomechanical model construction
- FEM simulation - i.e. breast compression
- Image (i.e. phantom) - reconstruction.

### Dependencies

This software was written using c++ and, mainly, the [itk](https://itk.org/) and [vtk](https://vtk.org/) libraries. 
> :warning: Currently, the software uses vtk v.7 (i.e. datafile version v.2). 
> However, the latest vtk version v.9 uses another file datafile version (v.5.1), and it may affect to the performance. 
> **I'm checking this issue**


The biomechanical model extraction is performed using the library [iso2mesh](https://github.com/fangq/iso2mesh) and, therefore, we need to use either Matlab$^{TM}$ or Octave$^{TM}$.
We are trying a code based on the [cgal library](https://www.cgal.org/) to directly transform itk images into 3D meshes.

The FE compression is performed using [NiftySim v.2.5.1](https://sourceforge.net/projects/niftysim/).
> :warning: This is the most important but the most outdated part.


### Current situation

The current code was not directly written to create bCT phantoms but as a part of a multimodal registration pipeline. Therefore, it may fail in several parts.

Currently, I'm Checking and sharing the code.

## Pipeline

- Biomechanical Model/Rewrite Information --> Convert `breast.tiff` images into `breast.nrrd`images

```bash
./RewriteInformation breast_segmentation.tiff output_path.nrrd voxelSizeX voxelSizeY voxelSizeZ

```
  - output:
    - `breast_image.nrrd` original breast segmentation image with spacial information
    - `resampled_image.nrrd` breast mask with $0.273\times0.273\times0.273~mm^3$

- Biomechanical Model/Matlab Code --> Biomechanical model extraction using `iso2mesh`.

> :warning: **Check the folder!**

  - output:
    - `mesh.vtk` Contains the FE breast model (with boundary Conditions and element labels???) 

- FE Simulation/Write Nifty Model
