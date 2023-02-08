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
  - `breast-resampled.nrrd` breast mask (binary, background vs. breast, I think), resampled to $0.273~mm^{3}$ due to historical issues.

### Matlab code

```
I need to check
```

The Matlab code is used to obtain the biomechanical models from the `.nrrd` images. 
If the original image has isotropic voxels, 