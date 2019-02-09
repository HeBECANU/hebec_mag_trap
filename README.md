# hebec_mag_trap
Magnetostatic field calculator for trapping neutral atoms
- Does gradient decent to find the trap center.
- Finds the trap depth via jacobian magnitude decent to find the critial points of the poential.

This code is to calculate the magnetic feild for a magnetic trap
 https://doi.org/10.1016/j.optcom.2006.09.031
 
  ![mag trap potential](/plots/potential_xz.png "Potential")
  
 ![The sum of the absolute magnitudes of the jacobians](/plots/jacobian_xz.png "Jacobinan Landscape")
 
 ## Bryce Done
 -
 
## To Do
- 3d minimizations
- coherent output struct
- more field sources
  - Current line
  - current spiral
- Add all open source contributions to the ack.

- add in multiple potential sources
  - gravity
  - optical
- full spiral integeration
- coil geometery ploting


## Contributions
This project would not have been possible without the many open source tools that it is based on.
* ***John D'Errico*** [Adaptive Robust Numerical Differentiation](https://au.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation)
* ***Ander Biguri*** [Perceptually uniform colormaps](https://au.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)
* ***Georg Stillfried** [mArrow3.m - easy-to-use 3D arrow](https://au.mathworks.com/matlabcentral/fileexchange/25372-marrow3-m-easy-to-use-3d-arrow)

