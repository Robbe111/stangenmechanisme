%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Voorbeeldanalyse van een vierstangenmechanisme.
%
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data initialization (all data is converted to SI units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% program data
fig_kin_4bar = 1;        % draw figures of kinematic analysis if 1
fig_dyn_4bar = 1;        % draw figures of dynamic analysis if 1

% kinematic parameters (link lengths)
L1 = 1; 
L2 = 1;
L3 = 0.3;
l3 = 0.15;
L4 = 0.35;
l4 = 0.15;
L5 = 2;
l5 = 1.5;
l7 = 1.2;
L6 = 1.2;
L8 = L5 - l5;
L7 = 1.5;

phi1 = 0;



%dynamic parameters, defined in a local frame on each of the bars.


m2 = L2*1.76;
m3 = 3;
m4 = 3;
m5 = L5*1.76;
m6 = L6*1.76;
m7 = L7*1.76;
m8 = L8*1.76;

a3 = L3*sin(deg2rad(10));
h3 = L3*cos(deg2rad(10));
b3 = l3;

a4 = -L4*sin(deg2rad(10));
h4 = L4*cos(deg2rad(10));
b4 = b3;

J2 = m2*L2^2/12;
J3 = m3*(h3*b3^3+h3*a3*b3^2+h3*a3^2*b3+b3*h3^3)/12;
J5 = m5*L5^2/12;
J4 = m4*(h4*b4^3+h4*a4*b4^2+h4*a4^2*b4+b4*h4^3)/12;
J6 = m6*L6^2/12;
J7 = m7*L7^2/12;
J8 = m8*L8^2/12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial angles
phi2_init = 0;
phi4_init = deg2rad(337);
phi5_init = deg2rad(25);
phi6_init = deg2rad(140);
phi7_init = deg2rad(140);
phi8_init = deg2rad(25);



t_begin = 0;                   % start time of simulation
t_end = pi;                    % end time of simulation
Ts = 0.01;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = 0.5;
A = 1;
phi3= 2.5 + omega * t;
dphi3= omega * ones(315,1) ;
ddphi3=zeros(315,1);

% calculation of the kinematics (see kin_4bar.m)
[phi2,phi4,phi5,phi6,phi7,phi8,dphi2,dphi4,dphi5,dphi6,dphi7,dphi8,ddphi2,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8] = ... 
                                    kinematics_4bar(Ts,L1,L2,L3,L4,l3,l4,L5,l5,L6,L7,l7,L8,...
                                            phi2_init,phi4_init,phi5_init,phi6_init,phi7_init,phi8_init,phi3,dphi3,ddphi3,t,fig_kin_4bar,phi1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculation of the dynamics (see dyn_4bar.m)
[F_P_x,F_Q_x,F_R_x,F_S_x,F_T_x,F_U_x,F_V_x,F_W_x,F_X_x,F_P_y,F_Q_y,F_R_y,F_S_y,F_T_y,F_U_y,F_V_y,F_W_y,F_X_y,M_P] = dynamics_4bar(phi2,phi3,phi4,phi5,phi6,phi7,phi8,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8 ...
                                                                                                    ,L1,L2,L3,L4,l3,l4,L5,l5,L6,L7,l7,L8, ...
                                                                                                        m2,m3,m4,m5,m6,m7,m8,J2,J3,J4,J5,J6,J7,J8,t,fig_dyn_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% load fourbar_movie Movie
% movie(Movie)

