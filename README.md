# hebec_mag_trap
**[David K. Shin](https://github.com/spicydonkey) [Bryce M. Henson](https://github.com/brycehenson), [Kieran F. Thomas](https://github.com/KF-Thomas)**  
Magnetostatic field calculator for trapping neutral atoms
- Multiple magnetic feild sources
  - arbitrary path numeric solver
  - loop (analytic)
  - line (analytic)
  - helix (numeric)
- Trap Charaterization
  - gradient decent to find the trap center.
  - finds trap frequency
  - Finds the trap depth via jacobian magnitude decent to find the critial points of the poential.
- Plots
  - 2d slice
  - 3d potential
- Traps implmeneted
  - Biquic https://doi.org/10.1016/j.optcom.2006.09.031
 
   ![mag trap potential](/plots/biquic.png "Potential")
   
  ![mag trap potential](/plots/potential_xz.png "Potential")
  
 ![The sum of the absolute magnitudes of the jacobians](/plots/jacobian_xz.png "Jacobinan Landscape")

## Install
``` 
git clone --recurse-submodules -j8 https://github.com/spicydonkey/hebec_mag_trap.git
```
then to update 
```
git submodule update --remote --merge
```

## To Do
- find the principle axes of the trap (dont just assume cartesian)
- implement other trap geometeries
- coherent output struct
- magnetic field source code
  - flat spiral/ pancake coil 
- add in multiple potential sources
  - gravity
  - optical

## Contributions
This project would not have been possible without the many open source tools that it is based on.
* ***John D'Errico*** [Adaptive Robust Numerical Differentiation](https://au.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation)
* ***Ander Biguri*** [Perceptually uniform colormaps](https://au.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)
* ***Georg Stillfried** [mArrow3.m - easy-to-use 3D arrow](https://au.mathworks.com/matlabcentral/fileexchange/25372-marrow3-m-easy-to-use-3d-arrow)

