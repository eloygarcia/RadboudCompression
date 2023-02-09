# Breast Phantom Generation 

This document summarizes the whole pipeline, from CT images to compressed breast phantoms, by means of a FEM simulation.

The algorithm is divided into 3 main parts:

- Biomechanical model construction
- FEM simulation - i.e. breast compression
- Image (i.e. phantom) reconstruction.


### Dependencies

This software was written using c++ and, mainly, the [itk](https://itk.org/) and [vtk](https://vtk.org/) libraries. 
> :warning: Currently, the software uses vtk v.7 (i.e. datafile version v.2). 
> However, the latest vtk version v.9 uses another file datafile version (v.5.1), and it may affect to the performance. 
> **I'm checking this issue**


The biomechanical model extraction is performed using the library [iso2mesh](https://github.com/fangq/iso2mesh) and, therefore, we need to use either Matlab$^{TM}$ or Octave$^{TM}$.
We are trying a code based on the [cgal library](https://www.cgal.org/) to directly transform itk images into 3D meshes.

The FE compression is performed using [NiftySim v.2.5.1](https://sourceforge.net/projects/niftysim/).
> :warning: This is the most important but the most outdated part.

NiftySim and the reconstruction algorithm uses CUDA, to speed up the process.


### Current situation

The current code was not directly written to create bCT phantoms but as a part of a multimodal registration pipeline. Therefore, it may fail in several parts.

Currently, I'm Checking and sharing the code.



## Pipeline

### Biomechanical Model/Rewrite Information 

Convert `breast.tiff` images into `breast.nrrd`images

```bash
./RewriteInformation breast_segmentation.tiff output_path.nrrd voxelSizeX voxelSizeY voxelSizeZ

```
  - output:
    - `breast_image.nrrd` original breast segmentation image with spacial information
    - `resampled_image.nrrd` breast mask with $0.273\times0.273\times0.273~mm^3$


### Biomechanical Model/Matlab Code 

Biomechanical model extraction using `iso2mesh`.

> :warning: **Check the folder!**

  - output:
    - `mesh.vtk` Contains the FE breast model (with boundary Conditions and element labels???) 


### FE Simulation/Write Nifty Model

Function `WriteNiftyModel` writes a `.xml` file to perform the FE simulation using *NiftySim*.
Some arguments are defined as default. 
**I should improve this**

```bash
./WriteNiftyModel <biomechanical-mesh.vtk> <breast-image.nrrd> <path-to-nifitysim-model.xml> <breast-thickness>
```

- inputs:
  - `biomechanical-mesh.vtk` or `mesh.vtk`obtained from the previous step
  - `breast-image.nrrd` obtained from the first step
  - `path-to-niftysim-model.xml` 
  - `breast-thickness` 

- output:
  - `niftysim-model.xml`

> Support +20 mm
> 
> Paddle BB-20-thickness (check!)


NiftySim usage:

```bash
./niftysim -x <path-to-niftysim-model.xml> -v -t -sport -export-mesh <path-to-output-mesh.vtk>
```
- arguments:
  - `-v` verbose
  - `-t` print time
  - `-sport` GPU (**MANDATORY**)
  - `-export-mesh`

- output:
  - `output-mesh.vtk`

The output `.vtk` mesh contains the node displacement as cell data. 
**The newest vtk version may affect this step and the following**.


### Phantom Reconstruction/Reconstruct Image

Function `ReconstructImage` perform the final ...

```bash
./ReconstructImage <path-to-niftysim-output-mesh.vtk> <path-to-breast-image.nrrd> <output-compressed-breast-mesh.vtk> <output-phantom-image.nrrd> 
```

- output
  - `compressed-breast-mesh.vtk`
  - `phantom-image.nrrd`