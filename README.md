# hebec_mag_trap
Magnetostatic field calculator for trapping neutral atoms
- Does gradient decent to find the trap center.
- Finds the trap depth via jacobian magnitude decent to find the critial points of the poential.

This code is to calculate the magnetic feild for a magnetic trap
 https://doi.org/10.1016/j.optcom.2006.09.031
 
  ![mag trap potential](/plots/potential_xz.png "Potential")
  
 ![The sum of the absolute magnitudes of the jacobians](/plots/jacobian_xz.png "Jacobinan Landscape")
 
 
##To Do
- Add all open source contributions to the ack.
- 3d minimizations
- add in multiple potential sources
  - gravity
  - optical
- full spiral integeration
- coil geometery ploting
