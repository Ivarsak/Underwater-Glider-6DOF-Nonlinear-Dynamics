# MATLAB simulation:  Spiral Motion of an Underwater Glider

This is a matlab script that simulates the dynamics of an underwater glider in three-dimensional space. The equations of motion are numerically solved using the ODE45 (Runge Kutta) matlab solver.
The main goal of the script is to simulate a steady spiral motion.
 The dynamic model and the hydrodynamic coefficients is based on the Seawing glider presented in [Zhang et al.](https://www.researchgate.net/publication/256817942_Spiraling_motion_of_underwater_gliders_Modeling_analysis_and_experimental_results) 
The paper also presents an analytical approach to the spiraling equilibria soultions trough an recrusive estimation. 

## Running the script
Run the Glider_main.m file with the Glider_function.m file open. 

## Plotting

* xyz position in earth frame cooardinates
* Pitch, Roll and Yaw angles (Theta, Phi, Psi)

### Characteristics of steady spiral motion

1. Pitch and Roll angles are constant.
2. Given that there is no control inputs, the initial conditions should remain the same for pitch and roll 
3. The change rate of the velocitiy V_dot should converge to 0 (at equilibria)
4. The ballast mass is held constant, ballast rate is zero -> mb_dot = 0

