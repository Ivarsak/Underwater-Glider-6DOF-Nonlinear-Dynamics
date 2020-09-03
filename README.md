# MATLAB simulation:  Spiral Motion of an Underwater Glider

This is a matlab script that simulates the dynamics of an underwater glider in three-dimensional space. The equations of motion are numerically solved using the ODE45 (Runge Kutta) matlab solver.
The main goal of the script is to simulate a steady spiral motion.
 The dynamic model and the hydrodynamic coefficients is based on the Seawing glider presented in [Zhang et al.](https://www.researchgate.net/publication/256817942_Spiraling_motion_of_underwater_gliders_Modeling_analysis_and_experimental_results) 

### Characteristics of steady spiral motion

1. Pitch and Roll angles are constant, while the yaw angle changes at a constant rate. Given these conditions the glider trajectory becomes a circular helix shown in the figure below. 
2. The change rate of translational and angular velocitiy converges towards 0.

-Simulations are carried out with fixed internal actuators and net bouyancy.  

<img src="https://user-images.githubusercontent.com/59923925/78356409-cc90f500-75af-11ea-978f-1beea733b0b1.png" width="500" height="500">

