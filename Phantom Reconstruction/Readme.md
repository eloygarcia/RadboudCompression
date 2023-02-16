# Image reconstruction

After performing the FE compression, we have to obtain the final compressed images

> CUDA REQUIRED!
> 
> Compilation issues! Change to:
> 
> NIFTYSIM_CUDA_ARCH = compute_80, NIFTYSIM_CUDA_CODE sm_80 (GPU RTX3090)
> 
> CUDA_SDK_ROOT_DIR : NVIDIA_CUDA_11.3_Samples (**Eloy: Ordenador grande. En el laptop es CUDA 12**)
> 
> CUDA_SDK_COMMON_INCLUDE_DIR : NVIDIA_CUDA_11.3_Samples/Common
> 
> CUDA_SEPARABLE_COMPILATION: FALSE
> 
> USE_GPU_GP_CONTACT : TRUE
> 
> USE_NAMESPACE_STD

> Check gcc --version to compile Niftsim (**Eloy: Sobremesa 7.5.0)
> 
> Use Docker?

## Extract Compressed Mesh

Function `extractCompressedMesh` is used to obtain the compressed mesh after the niftysim simulation

```bash
./getCompressedMesh <path-to-niftysim-output-mesh.vtk> <output-compressed-breast-mesh.vtk> 
```
- input:
  - `niftysim-output-mesh.vtk` mesh after compression
  - `<output-compressed-breast-mesh.vtk>` path to save the compressed breast mesh

## Reconstruct Image

Function `ReconstructImage` is used to obtain the compressed image.

> :warning: This function was tested using vtk filedata version 2.0 but using vtk v9.0

```bash
./ReconstructImage <iso2mesh-breast-mesh.vtk> <output-compressed-breast-mesh.vtk> <path-to-breast-image.nrrd>  <output-phantom-image.nrrd> 
```

- inputs: 
  - `<iso2mesh-breast-mesh.vtk>` path to the original mesh
  - `<output-compressed-breast-mesh.vtk>` path to save the compressed breast mesh
  - `breast-image.nrrd` original breast image
  - `<output-phantom-image.nrrd>`path to save the phantom

- output:
  - `compressed-breast-mesh.vtk`
  - `phantom-image.nrrd`