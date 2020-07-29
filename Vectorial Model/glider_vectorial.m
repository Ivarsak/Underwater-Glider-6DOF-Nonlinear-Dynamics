clear;
clc;


e1 = [1 0 0]';                     %Body frame unit vector along x
e2 = [0 1 0]';                     %Body frame unit vector along y
e3 = [0 0 1]';                     %Body frame unit vector along z

I3 = eye(3);                      %Identity matrix (3x3)

%% Masses

m_s = 54.28;                     %Hull mass [Kg]
m_p = 11;                        %Longitiudal and lateral moving mass [kg]
m = 65.28;                       %Displacement of fluid mass [Kg]
m_b = 0.5;                       %Net buoyancy mass [kg]
m_t = m_s + m_p;                 %Total mass
g = 9.81;                        %Gravity accleration

I_s = diag([0.60 15.27 15.32]);  %Inertia of static block/hull
I_rm = diag([0.02 10.16 0.17]);  %Inertia of movable block mass
Rr = 0.014;                      %Movable mass offset [m]

%% Added mass Matrix M_A


M_A = [diag([1.48 49.58 10.18]) zeros(3);
       zeros(3) diag([0.53 7.88 10.18]) ];

   
%% Hydrodynamic coefficients for forces and moments

KD0 = 7.19;                     %Coefficient of Drag force [X]
KD = 386.29;                    %Coefficient of Drag force [X]
K_beta = -115.65;               %Coefficient of Sideforce  [Y]
KL0 = -0.36;                    %Coefficient of Lift force [Z]
K_alpha = 440.99;               %Coefficient of Lift force [Z]
KMR = -58.27;                   %Coefficient of moment [K]
Kp = -19.83;                    %Coefficient of moment [K]
KM0 = 0.28;                     %Coefficient of moment [M]
KM = -65.84;                    %Coefficient of moment [M]
Kq = -205.64;                   %Coefficient of moment [M]
KMY = 34.10;                    %Coefficient of moment [N]
Kr = -389.30;                   %Coefficient of moment [N]

%% RK ODE solver


%Initalizing ballast mass and mass block positions

Gamma_d = deg2rad(40);              %Servo Angle [Rad]
rp1_d = 0.4216;                     %Moving mass block position [m]
mb_d = 0.3;                         %Ballast mass [Kg]

%Vectors to store ODE data in
t_data = [];
r_data = [];

      
         y0 = [
               [0 0 0]'       %Roll, Pitch, Yaw 
               [0 0 0]'       %XYZ position in earth frame
               [0 0 0]'       %Angular velocities
               [0.01 0 0]'    %Velocity 
               [0 Gamma_d 0]' %Position of rolling mass block
               [rp1_d 0 0]'   %Position of moving mass block 
               mb_d];            %Ballast mass
               
 
 
    tspan = [0 100];                  %Interval for the ODE solver
    
       
    fprintf('Integration interval: Tspan = %2.1f \n', tspan(2));
    fprintf('\n');
    fprintf('rp_x = %2.3f m \n', rp1_d);
    fprintf('Gamma = %2.1f degrees \n', rad2deg(Gamma_d));
    fprintf('mb = %2.1f kg \n', mb_d);
   

  save('Glider_variables.mat');

    %Runge Kutta solver
    [E,I] = ode45(@Glider_function_vectorial,tspan,y0);
  

    t_data = [t_data;E];  %ODE time - 1 variable 
    r_data = [r_data;I];  %ODE data - 19 variables  

%% Plotting
 %Plotting surge velocity u, heave velocity v and sway velocity w [m/s]
 
figure(1)
subplot(5,1,1)                     
plot(t_data, r_data(:,10))
ylabel('\it [m/s]','FontSize',9)                         
grid on
title('Glider velocities \bf \it [u v w]') 

hold on
plot(t_data, r_data(:,11))
plot(t_data, r_data(:,12))

legend('\it u [m/s]','\it v [m/s]','w [m/s]','FontWeight','bold');
hold off
%-------------------------------------------------------------------%
%Plotting roll angle [deg]
subplot(5,1,2)                      
plot(t_data, r_data(:,1)*180/pi)
ylabel(' [\circ]','FontSize',10)
grid on
title('Roll angle \bf \Phi','FontSize',10) 
%-------------------------------------------------------------------%
%Plotting pitch angle [deg]
subplot(5,1,3)                       
plot(t_data, r_data(:,2)*180/pi)
ylabel(' [\circ]','FontSize',10)
grid on
title('Pitch angle \bf \Theta','FontSize',10) 
%-------------------------------------------------------------------%
%Plotting Yaw angle [deg]
subplot(5,1,4)
plot(t_data, r_data(:,3)*180/pi)
ylabel('[\circ]','FontSize',10)                   
grid on
title('Heading \bf \Psi','FontSize',10) 

%-------------------------------------------------------------------%
%Plotting angular velocites - Omega1, Omega2, Omega3 [rad/s]
subplot(5,1,5)
plot(t_data, r_data(:,7))
ylabel('\it [rad/s]','FontSize',10)                       
grid on
title('Angular velocities \Omega_{1}, \Omega_{2}, \Omega_{3}','FontSize',10) 

hold on
plot(t_data, r_data(:,8))
plot(t_data, r_data(:,9))

legend('\Omega_{1} [rad/s]','\Omega_{2} [rad/s]','\Omega_{3} [rad/s]');
hold off
%-------------------------------------------------------------------%

figure(2)
plot3(r_data(:,4),r_data(:,5), r_data(:,6))
ylabel('Northings [m]','FontSize',10)
xlabel('Eastings [m]','FontSize',10)
zlabel('Down [m]','FontSize',10)

grid on

