# Breast biomecanical model extraction

First step is to build a biomechanical model of the breast, from CT segmented images.


## Rewrite information

```bash
./RewriteInformation breast.tiff breast.nrrd voxelSizeX voxelSizeY voxelSizeZ
```

In our case, we start using `.tiff` images, without physical (i.e. voxel size) information. 
To obtain a realistic breast shape, we need to know the actual physical shape.
So, we use the itk file format (i.e `.nrrd`) in order to preserve this information.
This function transform was created to transform the images.

- Input 
  - `breast.tiff` path to `.tiff` image.
  - `breast.nrrd` path to the output image.
  - `voxelSizeX voxelSizeY voxelSizeZ` original physical shape of the `.tiff` image.

- Output
  - `breast.nrrd` original segmentation image, in `.nrrd` format
  - `breast-resampled.nrrd` breast mask, resampled to 