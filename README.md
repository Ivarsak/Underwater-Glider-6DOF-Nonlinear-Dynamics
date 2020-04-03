# MATLAB simulation:  Spiral Motion of an Underwater Glider

This is a matlab script that simulates the dynamics of an underwater glider in three-dimensional space. The equations of motion are numerically solved using the ODE45 (Runge Kutta) matlab solver.
The main goal of the script is to simulate a steady spiral motion.
 The dynamic model and the hydrodynamic coefficients is based on the Seawing glider presented in [Zhang et al.](https://www.researchgate.net/publication/256817942_Spiraling_motion_of_underwater_gliders_Modeling_analysis_and_experimental_results) 
The paper also presents an analytical approach to the spiraling equilibria soultions trough a recrusive estimation. 

## Running the script
Run the Glider_main.m file with the Glider_function.m file open. 

## Plotting
The script plots the following

* xyz position in earth frame coordinates
* Pitch, Roll and Yaw angles (Theta, Phi, Psi)
* Velocity

### Characteristics of steady spiral motion

1. Pitch and Roll angles are constant, while the yaw angle changes at a constant rate. Meaning that the moving mass block is fixed for a desired (equilibria point) pitch and roll angle.
2. Given that there is no control inputs, the initial conditions should remain the same for pitch and roll over time 
3. The change rate of translational and angular velocitiy V_dot and Omega_dot respectively, should converge towards 0.
4. The ballast mass is fixed for the equilibria mass, which implies that the ballast rate is zero -> mb_dot = 0

![Spiralmotion](https://user-images.githubusercontent.com/59923925/78355221-a5d1bf00-75ad-11ea-8c20-72d2629f16de.png)
<img src="https://user-images.githubusercontent.com/59923925/78355221-a5d1bf00-75ad-11ea-8c20-72d2629f16de.png" width="100" height="100">

Note: This only validates for either a downward or upward spiral motion. It does not account for the transition between the two states, which require the ballast mass and the linear moving mass to change.
