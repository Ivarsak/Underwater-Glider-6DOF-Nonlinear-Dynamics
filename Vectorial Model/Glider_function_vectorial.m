function T=Glider_function_vectorial(t,y)

load('Glider_variables');

R =         y(1:3);                 %R(1) = Phi (roll), R(2) = Theta (pitch), R(3) = Psi (yaw)
b =         y(4:6);                 %b(1) = x, b(2) = y, b(3) = z
omega =     y(7:9);                 %omega(1) = p, omega(2) = q, omega(3) = r
V =         y(10:12);               %V(1) = V1, V(2) = V2, V(3) = V3 in body frame
r_gamma =   y(13:15);               %rp(1) = rp1, rp(2) = rp2, rp(3) = rp3 position of moving mass block
rp =        y(16:18);               %r_gamma(1) = 0, r_gamma(2) = Gamma, r_gamma(3) = 0 position of rolling mass block
mb =        y(19);                  %mb -> ballast mass

u = V(1);     %Sway velocity 
v = V(2);     %Heave velocity
w = V(3);     %Surge velocity 

p = omega(1); %Angular velocity around roll angle, Phi
q = omega(2); %Angular velocity around pitch angle, Theta
r = omega(3); %Angular velocity around yaw angle, Psi

u_gamma = [0 0 0]';
u_mp    = [0 0 0]';
u_mb    = 0;

tau = [0 0 0 u_gamma(2) u_mp(1) 0]';  %Control inputs

m_t = mb + m_p + m_s;
 
%% Kinematics

% Rotation matrix to map the Euler Angle attitudes
n_1 = [1 sin(R(1))*tan(R(2)) cos(R(1))*tan(R(2));
       0    cos(R(1))            -sin(R(1));
       0 sin(R(1))*sec(R(2))    cos(R(1))*sec(R(2))];
   

%Rotation matrix to map the translational velocities in the body frame
n_2 = [cos(R(3))*cos(R(2))    (-sin(R(3))*cos(R(1)) + cos(R(3))*sin(R(2))*sin(R(1))) (sin(R(3))*sin(R(1)) + cos(R(3))*cos(R(1))*sin(R(2)));
     sin(R(3))*cos(R(2))    (cos(R(3))*cos(R(1)) + sin(R(1))*sin(R(2))*sin(R(3)))  (-cos(R(3))*sin(R(1)) + sin(R(2))*sin(R(3))*cos(R(1)));
     -sin(R(2))            cos(R(2))*sin(R(1))                                 cos(R(2))*cos(R(1))];


%Inertia matrix for the rolling mass
R_x = [1 0 0;
       0 cos(r_gamma(2)) -sin(r_gamma(2));
       0 sin(r_gamma(2)) cos(r_gamma(2))];
   
%% Computing Hydrodynamic Damping Matrix D(v)
  
   v_r = sqrt(u^2 + v^2 + w^2);
   
   Alpha = atan(w/u);                     %Angle of attack
   Beta  = asin(v/v_r);                   %Sideslip angle
   
   X = (KD0 + (KD*(Alpha^2)))*v_r^2;       %Drag force
   Y = (K_beta*Beta)*v_r^2;                %Sideforce
   Z = (KL0 + K_alpha*Alpha)*v_r^2;        %Lift force
   
   K = (KMR*Beta + Kp*p)*v_r^2;            %Roll Moment
   M = (KM0 + KM*Alpha + Kq*q)*v_r^2;      %Pitch Moment
   N = (KMY*Beta + Kr*r)*v_r^2;            %Yaw Moment
   
   F_hy = [X -Y Z]';                        %Total Hydrodynamic force in flow frame
   M_hy = [-K -M -N]';                      %Total hydrodynamic moment in flow frame

   
 %Rotation matrix from the flow frame to body frame
    R_BF = [cos(Alpha)*cos(Beta) -cos(Alpha)*sin(Beta) -sin(Alpha);
            sin(Beta) cos(Beta) 0;
            sin(Alpha)*cos(Beta) -sin(Alpha)*sin(Beta) cos(Alpha)];
     
  % Damping Matrix
   D_v = [R_BF*F_hy; R_BF*M_hy];
   

 %% Computing Restoring Forces
 
F_G = [0 0 m_t*g]';      %Gravitational force in NED {n}
F_B = [0 0 -m*g]';       %Buoyant force in NED {n}

%Computing the CG and CB vectors in body-frame {b}
r_mp = rp(1)*e1 + Rr*(cos(r_gamma(2) + pi/2)*e2 + sin(r_gamma(2) + pi/2)*e3); %Position of moving mass/rolling mass block
r_mb = [0 0 0]';                                                              %Position of internal ballast tank (fixed)
r_s = [-0.0814 0 0.0032]';                                                    %Position of static mass (fixed)

r_CB = [0 0 0]';                             %CB vector. Coincides with body-fixed origin O_{b} = [0,0,0]^T
r_CG = (r_mp*m_p + r_mb*m_b + r_s*m_s)/m_t   %CG vector. Which is offset from CB = [0,0,0] in the body frame {b}

rx_CG = r_CG(1);                             %CG in x-axis
ry_CG = r_CG(2);                             %CG in y-axis
rz_CG = r_CG(3);                             %CG in z-axis

%Computing the restoring terms g(n) in the body-frame {b}
 G_n = -[n_2'*(F_G + F_B);
        cross(r_CG, n_2'*F_G) + cross(r_CB, n_2'*F_B)];
           
%% Computing Rigid-Body Dynamics
 %-------------------------------------------------------------------------------------------------------------------
 %------------------------------------------------------------------------------------------------------------------- 
 %Computing the Total System Inertia Matrix - M
 
 I_mp = (R_x'*I_rm*R_x);                    %Inertia of the moving/rolling mass block m_p about CB/CO
 I_rb = I_s + I_mp;                         %Total inertia about CB
 
 M_rb = [m_t*I3 -m_t*S(rx_CG,ry_CG,rz_CG);  %Rigid-Body system Inertia Matrix
         m_t*S(rx_CG,ry_CG,rz_CG) I_rb];
     
 M = M_rb + M_A;                            %Total System Inertia Matrix  
 %-------------------------------------------------------------------------------------------------------------------
 %-------------------------------------------------------------------------------------------------------------------
 % Computing the Coriolis and Centripetal Matrix C(v)
 
 IO = I_rb*omega;      
 
 C_v = [zeros(3) (-m_t*S(u,v,w) - m_t*S(p,q,r)*S(rx_CG,ry_CG,rz_CG));
        (-m_t*S(u,v,w) + m_t*S(rx_CG,ry_CG,rz_CG)*S(p,q,r)) -S(IO(1), IO(2), IO(3))];
    

%% Computing the Equations of Motion

 v_dot = inv(M)*(-D_v - C_v*[u v w p q r]' - G_n + tau); %6DOF dynamics
 
R_dot =     n_1*omega;                         %Kinematic equation  
b_dot =     n_2*V;                             %Kinematic equation 
V_dot =     [v_dot(1) v_dot(2) v_dot(3)]'      %Dynamic equation 
omega_dot = [v_dot(4) v_dot(5) v_dot(6)]'      %Dyamic equation 

%Control inputs

rp1_dot =   u_mp;         %change rate = [0 0 0]'
Gamma_dot = u_gamma;      %change rate = [0 0 0]'
mb_dot =    u_mb;         %change rate = 0


%Return vector for the ODE45 solver
T =[R_dot' b_dot' omega_dot' V_dot' Gamma_dot' rp1_dot' mb_dot']'; 

    
end

