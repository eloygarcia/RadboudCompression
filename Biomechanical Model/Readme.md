# Breast biomecanical model extraction

First step is to build a biomechanical model of the breast, from CT segmented images.


## Rewrite information

```bash
./RewriteInformation breast.tiff breast.nrrd voxelSizeX voxelSizeY voxelSizeZ
```
 
In our case, we start using `.tiff` images, without physical (i.e. voxel size) information. 
To obtain a realistic breast shape, we need to know the actual physical shape.
So, we use the itk file format (i.e `.nrrd`) in order to preserve this information.
This function was created to use isotropic and anisotropic voxel size, preserving the information.
However, to reduce the complexity, we use isotropic voxels to extract the biomechanical model.
Therefore, we obtain both, the original breast segmentation and a resampled version.

- Input 
  - `breast.tiff` path to `.tiff` image.
  - `breast.nrrd` path to the output image.
  - `voxelSizeX voxelSizeY voxelSizeZ` original physical shape of the `.tiff` image.

- Output
  - `breast.nrrd` original segmentated image, in `.nrrd` format
  - `breast-resampled.nrrd` resampled to isotropic voxels ($0.273\times0.273\times0.273~mm^3$) due to historical issues.



### Matlab code

> > :warning: **REQUIRED [iso2mesh](https://github.com/fangq/iso2mesh)!!**
  
```
addpath( $iso2mesh )

im = uint8( logical( nrrdread( breast_mask.nrrd ));
im = permute( im, [2,1,3]);

spacing = 0.273; 

opt = varargin2struc('radbound', 10, 'distbound', 10, 'maxnode', 20); % to Koen
% opt = varargin2struc('radbound', 15, 'distbound', 15, 'maxnode', 30); % to Marco (why?)

[nodes, elements, faces] = cgalv2m(im, opt, 10);
[nodes, elemnents] = clearmesh( nodes, elements(:,1:4) );
nodes = nodes * spacing;

boundaryConditions = nodes(:,1) < (min(nodes(:,1))+2);

writeVTK( `path_to_mesh.vtk`, elements(:,1:4), nodes(:,1:3), boundaryConditions); 
```

The Matlab code is used to obtain the biomechanical models from the `.nrrd` images. 
If the original image has isotropic voxels, we could use the original `.tiff` images.

All features were optimized to obtain a suitable number of elements using isotropic voxels.

- Inputs & options
  - `breast.nrrd` path to the `.nrrd` image (I think)
  - `spacing` due to historical issues, it is defined to $0.273\times0.273\times0.273~mm^3$
  - `opt` arguments to define the mesh. I've found two versions (to Koen and to Marco) but I don't know which is the best
  - `cgalv2m` the `iso2mesh`function used to create the mask 
  - `clearmesh` we need to avoid isolated nodes on the mesh
  - `boundaryConditions` depends on the breast orientation
  - `writeVTK` is a function used to write the `.vtk` unstructured mesh with its information.

- Output
  - `phantom.vtk` the phantom mesh, which will be compressed using FE

The exposed pipeline need to be repeated for each image. 



### New Model Extraction

```bash
./ItkToCGAL inputfilename outputfilename arguments
```

This functions should be used to extract the biomechanical model without using `iso2mesh`  and `Matlab`.
However, we need the [`CGAL` Library](https://www.cgal.org/)

- Input 
  - `ìnputfilename`  path to the `.nrrd` or `.mhd` image (both correspond to itk image formats)
  - `outputfilename` path to the mesh file. I usually get a `.vtk` file, but (I do not why) this function is ready  to write `.inp` files (Abacus$^{TM}$ and ANSYS$^{TM]$ format)
  - `arguments` several arguments to modify the number of elements and quality of the final mesh. I need to check what each argument does.
    - `--value_ouside` (default 0)
    - `--facet_angle` (default 30)
    - `--facet_size` (default 5)
    - `--facet_distance` (default 1)
    - `--cell_radius_edge_ratio` (default 1)
    - `--cell_size` (default 5)

> :warning: **First try:** --facet_size 1 --facet_distance 0.5 --cell_size 0.75 --cell_radius_edge_ratio 2
> 
> no. points 39190
> 
> no. elements 221386
>
> VTK WRITER DOESN'T WORK!!!
> 
 

