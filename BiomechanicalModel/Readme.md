# Breast biomecanical model extraction

First step is to build a biomechanical model of the breast, from CT segmented images.


## Rewrite information

```bash
./RewriteInformation breast.tiff breast.nrrd voxelSizeX voxelSizeY voxelSizeZ
```

> :warning: **This fuction is NOT testes, because we do not have the original `.tiff` images
> **: Be very careful here!
> 
In our case, we start using `.tiff` images, without physical (i.e. voxel size) information. 
To obtain a realistic breast shape, we need to know the actual physical shape.
So, we use the itk file format (i.e `.nrrd`) in order to preserve this information.
This function transform was created to transform the images.

- Input 
  - `breast.tiff` path to `.tiff` image.
  - `breast.nrrd` path to the output image.
  - `voxelSizeX voxelSizeY voxelSizeZ` original physical shape of the `.tiff` image.

- Output
  - `breast.nrrd` original segmentated image, in `.nrrd` format
  - `breast-resampled.nrrd` breast mask (binary, background vs. breast, I think), resampled to isotropic voxels ($0.273\times0.273\times0.273~mm^3$) due to historical issues.

### Matlab code

```
I need to check
```

The Matlab code is used to obtain the biomechanical models from the `.nrrd` images. 
If the original image has isotropic voxels, we could use the original `.tiff` images.

All features were optimized to obtain a suitable number of elements using isotropic voxels.

- Input
  - `breast.nrrd` path to the `.nrrd` image (I think)

- Output
  - `phantom.vtk` the phantom mesh, which will be compressed using FE


### New Model Extraction

```bash
./ItkToCGAL inputfilename outputfilename arguments
```

> :warning: **This function was found during the current work and it was NOT testes** 
> Be very careful here!

This functions should be used to extract the biomechanical model without using `iso2mesh`  and `Matlab`.
However, we need the [`CGAL` Library](https://www.cgal.org/)

- Input 
  - `ìnputfilename`  path to the `.nrrd` or `.mhd` image (both correspond to itk image formats)
  - `outputfilename` path to the mesh file. I usually get a `.vtk` file, but (I do not why) this function is ready  to write `.inp` files (Abacus$^{TM}$ and ANSYS$^{TM]$ format)
  - `arguments` several arguments to modify the number of elements and quality of the final mesh. I need to check what each argument does.
    - `--value_ouside`
    - `--facet_angle`
    - `--facte_size`
    - `--facet_distance`
    - `--cell_radius_edge_ratio`
    - `--cell_size`