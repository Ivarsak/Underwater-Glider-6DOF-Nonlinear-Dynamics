


%Parameteres definition
ms = 54.28;                     %Static block mass   [Kg]
mr = 11;                        %Movable block mass  [Kg]
Rr = 0.014;                     %Movable block offset [m]
g = 9.81;                       %Gravity

rs1 = -0.0814;                  %Position of static block x
rs2 = 15.27;                    %Position of static block y
rs3 = 15.32;                    %Position of static block z

MA1 = 1.48;                     %Added mass term 1 
MA2 = 49.58;                    %Added mass term 2
MA3 = 65.92;                    %Added mass term 3

IA1 = 0.53;                     %Added inertia term 1
IA2 = 7.88;                     %Added inertia term 2
IA3 = 10.18;                    %Added inertia term 3

Isx = 0.60;                     %Inertia term static block 1
Isy = 15.27;                    %Inertia term static block 2
Isz = 15.32;                    %Inertia term static block 3

Irx = 0.02;                     %Inertia term movable block 1
Iry = 10.16;                    %Inertia term movable block 2
Irz = 0.17;                     %Inertia term movable block 3

I3 = eye(3);
M_t = ms + mr;

mt1 = M_t;
mt2 = M_t;
mt3 = M_t
%Input to nonlinear solver

V = 0.490                       %Total velcotiy [m/s]
Alpha = deg2rad(1.267)          %Angle of Attack [rad]
Beta = deg2rad(-1.283);         %Sideslip angle [rad];


%Hydrodynamic coefficients

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



D = (KD0 + (KD*Alpha^2))*V^2;                %Drag
SF = K_beta*Beta*V^2;                        %Sideforce
L = (KL0 + K_alpha*Alpha)*V^2;               %Lift
   
MDL1 = (KMR*Beta)*V^2;                      %Roll Moment
MDL2 = (KM0 + KM*Alpha)*V^2;                %Pitch Moment
%MDL_3 = (KMY*Beta + Kr*x(1))*V^2;          %Yaw Moment


%omega3  =  x(1)
%Theta   =  x(2)
%Phi     =  x(3)
%mb      =  x(4)
%rp      =  x(5)
%Gamma   =  x(6)


%Steady state function

F = @(x) [(V*x(1)*(mt2*sin(Beta)*cos(x(3))*cos(x(2)) - mt3*sin(Alpha)*cos(Beta)*sin(x(3))*cos(x(2))) + ...
    mr*x(1)^2*(x(5)*cos(x(2)).^2 + Rr*((sin(2*x(2)))/2)*cos(x(3) + x(6))) - x(4)*g*sin(x(2)) + ...
    ms*x(1)^2*(cos(x(2)).^2*rs1 + (sin(2*x(2))/2)*cos(x(3))*rs3) - ...
    D*cos(Alpha)*cos(Beta) - SF*cos(Alpha)*sin(Beta) + L*sin(Alpha));
    
    (V*x(1)*(-mt3*sin(Alpha)*cos(Beta)*sin(x(2)) - mt1*cos(Alpha)*cos(Beta)*cos(x(3))*cos(x(2)))...
    + x(4)*g*sin(x(3))*cos(x(2)) - D*sin(Beta) + SF*cos(Beta)...
    + ms*x(1)^2*(((sin(x(3))*sin(2*x(2))*rs1)/2) - ((sin(2*x(3))*cos(x(2)).^2*rs3)/2))...
    + mr*x(1)^2*(((x(5)*sin(2*x(2))*sin(x(3)))/2) - ((Rr*sin(2*x(3))*cos(x(2)).^2*cos(x(6)))/2))...
    + mr*x(1)^2*(-Rr*sin(x(6))*(cos(x(3)).^2*cos(x(2)).*2 + sin(x(2)).^2)));
    
    ((Iry - Irz - mr*Rr^2)*x(1)^2*((sin(2*x(6))*sin(x(3))*sin(2*x(2)))/4) + (MA3 - MA1)*V^2*cos(Beta).^2*((sin(2*Alpha))/2)...
    -mr*Rr*x(5)*x(1)^2*((sin(2*x(3))*cos(x(2)).^2*sin(x(6)))/2) - mr*V*x(1)*x(5)*(cos(Alpha)*cos(Beta)*sin(x(3))*cos(x(2)) + ...
    sin(x(2))*sin(Beta)) + (Isx-Isz+Irx+IA1-IA3-Iry*sin(x(6)).^2+(mr*Rr^2-Irz)*cos(x(6)).^2-mr*x(5)^2)*x(1)^2*cos(x(3))*((sin(2*x(2)))/2)...
    -mr*g*x(5)*cos(x(3))*cos(x(2))+ MDL1*sin(Beta) + MDL2*cos(Beta) + ms*x(1)^2*((rs3^2-rs1^2)*cos(x(3))*((sin(2*x(2)))/2)...
    + rs3*rs1*(cos(x(3)).^2*cos(x(2)).^2 - sin(x(2)).^2)) + ms*x(1)*V*sin(Beta)*(cos(x(3))*cos(x(2))*rs3-sin(x(2))*rs1)...
    -ms*g*(rs3*sin(x(2))+rs1*cos(x(3))*cos(x(2))) - ms*V*x(1)*sin(x(3))*cos(x(2))*cos(Beta)*(rs3*sin(Alpha) - rs1*cos(Alpha))...
    + mr*V*Rr*x(1)*cos(x(6))*cos(x(2))*(sin(Beta)*cos(x(3))-sin(x(3))*sin(Alpha)*cos(Beta))...
    + cos(x(6))*(-mr*g*Rr*sin(x(2)) + mr*Rr*x(5)*x(1)^2*(cos(x(3)).^2*cos(x(2)).^2-sin(x(2)).^2)))]

    
 
%Inital guess    
x0 = [0.003 , deg2rad(-30), deg2rad(-13), 0.3, 0.4, deg2rad(45)];

%Levenberg-Marquardt algorithm
[x, fval] = fsolve(F, x0)


omega3 = x(1)            %Angular velocity around z
Theta  = rad2deg(x(2))   %Pitch angle
Phi    = rad2deg(x(3))   %Roll angle
mb     = x(4)            %Buoyancy
rp     = x(5)            %Movable block position
gamma  = rad2deg(x(6))   %Movable block rotating position


