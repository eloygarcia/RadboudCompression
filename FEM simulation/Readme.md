# Finite Element Simulation

> :warning: This is the most important but the most outdated part.

We perform the breast compression using Finite Elements algorithm. 
The current version use [NiftySim v.2.5.1](https://sourceforge.net/projects/niftysim/).

We need to write the FE model and, later, perform the simulation.


### Some parameters to take into account

`auxiliarDefinitions.h` contains the material parameters (glandular and adipose tissue only). 
** I'm not sure if this is the latest version**

`bulk modulus`, `shear modulus`, and `Lame perameters` can be computed using the file `MechanicalProperties.h`.

> **Ref:** Wellman, P., Tactile Imaging, PhD thesis, Cambridge, MA: Harvard University’s Division of Engineering and Applied Sciences (1999)
- Young modulus, glandular tissue: $15.1~KPa$
- Young modulus, adipose tissue: $4.46~KPa$
- Poisson ratio: $0.499$ in the two cases

Element labels (i.e glandular or fatty) are obtained using just the barycenter of the element. 
Their physical position, within the original segmented image, define the label.

> **Check!** 	
> 
>disp.clear() ; disp.push_back(0); disp.push_back(+20); disp.push_back(0); 
> en NiftySimEjecutable.cpp (line 161). El soporte sube 2 cm mientras que el superior hace todo el trabajo debido al efecto del técnico.




## Write Nifty Model

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


## NiftySim


> To compile niftysim:
> 
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
