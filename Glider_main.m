
clear;
clc;

%See the readme file for info regarding the script

%Inertial frame

N = [1 0 0]';                     %Unit vector along x
E = [0 1 0]';                     %Unit vector along y
D = [0 0 1]';                     %Unit vector along z

I3 = eye(3);                      %Identity matrix (3x3)

%% Mass terms

mrb = 54.28;                     %Hull mass [Kg]
mp = 11;                         %Longitiudal and lateral moving mass [kg]
m = 65.28;                       %Displacement of fluid mass [Kg]
mb1 = 0.5;                       %Net bouyancy [Kg]
mv = mrb + mp + mb1;             %Total vehicle mass [Kg]
m0 = mv - m;                     %Buoyancy [kg]
g = 9.81;                        %Gravity accleration

%Added mass terms

mf1 = 1.48;                      %Added mass term 1
mf2 = 49.58;                     %Added mass term 2
mf3 = 65.92;                     %Added mass term 3

M_A = diag([mf1 mf2 mf3]);       %Added mass matrix

r_rb1 = -0.0814;                 %Position of static mass block e1 direction
r_rb2 = 0;                       %Position of static mass block e2 direction
r_rb3 = 0.0032;                  %Position of static mass block e3 direction

rrb = [r_rb1 r_rb2 r_rb3]';      %Position of statick block [m]
Rr = 0.014;                      %Movable mass offset position [m]

%% Inerta terms

I_x = 0.53;                      %Added inertia term 1 
I_y = 7.88;                      %Added inertia term 2
I_z = 10.18;                     %Added inertia term 3

I_A = diag([I_x I_y I_z]);       %Added inertia matrix

I_rb = diag([0.60 15.27 15.32]); %Inertia of static block/hull
I_rm = diag([0.02 10.16 0.17]);  %Inertia of movable block mass

%% Coupling term

C_1 = 2.57;                      %Added coupling term 1
C_2 = 3.61;                      %Added coupling term 2


C_A = [0 0 0;                    %Added coupled matrix
       0 0 C_2;
       0 C_1 0];
   
   
%% Hydrodynamic coefficients

KD0 = 7.19;                     %Coefficient of Drag force [D]
KD = 386.29;                    %Coefficient of Drag force [D]
K_beta = -115.65;               %Coefficient of Sideforce  [SF]
KL0 = -0.36;                    %Coefficient of Lift force [L]
K_alpha = 440.99;               %Coefficient of Lift force [L]
KMR = -58.27;                   %Coefficient of moment [MDL_1]
Kp = -19.83;                    %Coefficient of moment [MDL_1]
KM0 = 0.28;                     %Coefficient of moment [MDL_2]
KM = -65.84;                    %Coefficient of moment [MDL_2]
Kq = -205.64;                   %Coefficient of moment [MDL_2]
KMY = 34.10;                    %Coefficient of moment [MDL_3]
Kr = -389.30;                   %Coefficient of moment [MDL_3]


%% Equilibrium values for steady state spiral motion 

%These values are based on the simulation given in Zhang et al.

Gamma_d = deg2rad(45);              % Servo Angle [Rad]
Beta_d = deg2rad(-1.283);           % Sideslip anlge [Rad]
Alpha_d = deg2rad(1.267);           % Angle of attack [Rad]
Theta_d = deg2rad(-35.641);         % Pitch angle [Rad]
Phi_d = deg2rad(-13.703);           % Roll angle [Rad]   
rp1_d = 0.4216;                     % Moving mass block position [m]
mb_d = 0.3;                         % Ballast mass [Kg]
V_d = 0.490;                        % Velocity [m/s]
Omega_3_d = 0.0039;                 % Turn rate [rad/s]



v1_d = V_d*cos(Alpha_d)*cos(Beta_d);    %Initial velocity in x direction
v2_d = V_d*sin(Beta_d);                 %Initial velocity in y direction
v3_d = V_d*sin(Alpha_d)*cos(Beta_d);    %Initial velocity in z direction

       
%Save workspace for the glider function
save('Glider_variables.mat');

%% Runge Kutta


        %Initial conditions, this corresponds to the glidig Equlibria for a
        %steady spiral motion.

        y0 = [
               [Phi_d Theta_d 0]'     % Roll, Pitch, Yaw 
               [0 0 0]'               %XYZ position in earth frame [u v w]
               [0 0 Omega_3_d]'       %Angular velocities [p q r]
               [v1_d v2_d v3_d]'             %Velocity [v1 v2 v3]
               [0 Gamma_d 0]'         %Position of moving mass block
               [0 0 0]'               %Position of the ballast mass
               [rp1_d 0 0]'           %Position of rolling mass block 
                mb_d ];               %Ballast mass
           
   
    
    tspan = [0 100];                  %Interval for the ODE solver
    
    
    
    Radius = (V_d*cos(Theta_d-Alpha_d)/Omega_3_d);
    
    fprintf('Integration interval: Tspan = %f \n', tspan(2));
    fprintf('\n');
    fprintf('\n');
    fprintf('Values at equilibria (Initial conditions)\n');
    fprintf('\n');
    fprintf('Velocity = %2.3f m/s \n', V_d);
    fprintf('Theta = %2.1f degrees, Phi = %2.1f degrees, Psi = %2.1f degrees \n', rad2deg(Theta_d), rad2deg(Phi_d), 0);
    fprintf('Alpha = %2.1f degrees, Beta = %2.1f degrees \n', rad2deg(Alpha_d), rad2deg(Beta_d));
    fprintf('Omega3 = %f rad/s \n', Omega_3_d);
    fprintf('rp_x = %2.3f m \n', rp1_d);
    fprintf('Gamma = %2.1f degrees \n', rad2deg(Gamma_d));
    fprintf('mb = %2.1f kg \n', mb_d);
    fprintf('Turning radius = %2.1f m \n', Radius);
    
    %Runge Kutta method
    [E,I] = ode45(@Glider,tspan,y0);
     
              


%% Plotting

subplot(5,1,1)                     %3D plot with XYZ cooardinates in earth frame
plot3(I(:,4), I(:,5), I(:,6))
grid on
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')


subplot(5,1,2)                      %Plotting Pitch angle in degrees
plot(E, I(:,2)*180/pi)
ylabel('Theta [deg]')
grid on


subplot(5,1,3)                       %Plotting Roll angle in degrees
plot(E, I(:,1)*180/pi)
ylabel('Phi [deg]')
grid on


subplot(5,1,4)
plot(E, I(:,3)*180/pi)
ylabel('Psi [deg]')                   %Plotting Yaw angle in degrees
grid on


subplot(5,1,5)
plot(E, I(:,10))
ylabel('m/s')                          %Plotting angular velocity
grid on

hold on
plot(E, I(:,11))
plot(E, I(:,12))

legend('v1 [m/s]','v2 [m/s]','v3 [m/s]');
hold off




