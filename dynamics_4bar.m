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


function [F_P_x,F_Q_x,F_R_x,F_S_x,F_T_x,F_U_x,F_V_x,F_W_x,F_X_x,F_P_y,F_Q_y,F_R_y,F_S_y,F_T_y,F_U_y,F_V_y,F_W_y,F_X_y,M_P] = ...
dynamics_4bar(phi2,phi4,phi5,phi6,phi7,phi8,dphi2,dphi4,dphi5,dphi6,dphi7,dphi8,ddphi2,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8 ...
                                                                                                    ,L1,L2,L3,L4,l3,l4,L5,l5,L6,L7,l7,L8, ...
                                                                                                        m2,m3,m4,m5,m6,m7,m8,J2,J3,J4,J5,J6,J7,J8,t,fig_dyn_4bar)


% a lot of definitions to make the matrix A and B a bit clear.
% skip the definitions for now (move down to "force analysis")
% and check them when you need them.


% cogi_P_x, cogn_P_y = vector from the centre of gravity of bar i to point P
cog3_x = (L3*cos(phi3)+l3*cos(theta3))/3;
cog3_y = (L3*sin(phi3)+l3*sin(theta3))/3;
cog4_x = (L4*cos(phi4)+l4*cos(theta4))/3;
cog4_y = (L4*sin(phi4)+l4*sin(theta4))/3;


cog3_P_x = -cog3_x;
cog3_P_y = -cog3_y;
cog3_Q_x = L3*cos(phi3)-cog3_x;
cog3_Q_y = L3*sin(phi3)-cog3_y;
cog3_X_x = l3*cos(theta3)-cog3_x;
cog3_X_y = l3*sin(theta3) - cog3_y;
cog5_Q_x = -(L5/2)*cos(theta5);
cog5_Q_y = -(L5/2)*sin(theta5);
cog5_R_x = (l5-L5/2)*cos(theta5);
cog5_R_y = (l5-L5/2)*sin(theta5);
cog5_S_x = -cog5_Q_x;
cog5_S_y = -cog5_Q_y;
cog6_R_x = (L6/2)*cos(phi6);
cog6_R_y = (L6/2)*sin(phi6);
cog6_U_x = -(L6/2)*cos(phi6);
cog6_U_y = -(L6/2)*sin(phi6;
cog7_S_x = L7/2*cos(phi7);
cog7_S_y = L7/2*sin(phi7);
cog7_T_x = -L7/2*cos(phi7);
cog7_T_y = -L7/2*sin(phi7);
cog8_U_x = -L8/2*cos(phi8);
cog8_U_y = -L8/2*sin(phi8);
cog8_T_x = L8/2*cos(phi8);
cog8_T_y = L8/2*sin(phi8);
cog4_U_x = L4*cos(phi4)-cog4_x;
cog4_U_y = L4*sin(phi4)-cog4_y; 
cog4_V_x = l4*cos(theta4) - cog4_x; 
cog4_V_y = l4*sin(theta4) - cog4_y;
cog4_W_x = -cog4_x; 
cog4_W_y = cog4_y; 
cog2_X_x = -(L2/2)*cos(phi2); 
cog2_X_y = -(L2/2)*sin(phi2); 
cog2_W_x = (L2/2)*cos(phi2); 
cog2_W_y = (L2/2)*sin(phi2); 

% 3D omega (dphi) and alpha (ddphi) vectors)
omega2 = [zeros(size(phi3)) zeros(size(phi3)) dphi2];
omega3 = [zeros(size(phi3)) zeros(size(phi3)) dphi3];
omega4 = [zeros(size(phi3)) zeros(size(phi3)) dphi4];
omega5 = [zeros(size(phi3)) zeros(size(phi3)) dphi5];
omega6 = [zeros(size(phi3)) zeros(size(phi3)) dphi6];
omega7 = [zeros(size(phi3)) zeros(size(phi3)) dphi7];
omega8 = [zeros(size(phi3)) zeros(size(phi3)) dphi8];
alpha2 = [zeros(size(phi3)) zeros(size(phi3)) ddphi2];
alpha3 = [zeros(size(phi3)) zeros(size(phi3)) ddphi3];
alpha4 = [zeros(size(phi3)) zeros(size(phi3)) ddphi4];
alpha5 = [zeros(size(phi3)) zeros(size(phi3)) ddphi5];
alpha6 = [zeros(size(phi3)) zeros(size(phi3)) ddphi6];
alpha7 = [zeros(size(phi3)) zeros(size(phi3)) ddphi7];
alpha8 = [zeros(size(phi3)) zeros(size(phi3)) ddphi8];

% 3D model vectors
P_cog3_vec = [-cog3_P_x    -cog3_P_y    zeros(size(phi3))];
Q_cog3_vec = [-cog3_Q_x    -cog3_Q_y    zeros(size(phi3))];
X_cog3_vec = [-cog3_X_x    -cog3_X_y    zeros(size(phi3))];
X_cog2_vec = [-cog2_X_x    -cog2_X_y    zeros(size(phi3))];
W_cog2_vec = [-cog2_W_x    -cog2_W_y    zeros(size(phi3))];
U_cog4_vec = [-cog4_U_x    -cog4_U_y    zeros(size(phi3))];
V_cog4_vec = [-cog4_V_x    -cog4_V_y    zeros(size(phi3))];
W_cog4_vec = [-cog4_W_x    -cog4_W_y    zeros(size(phi3))];
Q_cog5_vec = [-cog5_Q_x    -cog5_Q_y    zeros(size(phi3))];
R_cog5_vec = [-cog5_R_x    -cog5_R_y    zeros(size(phi3))];
S_cog5_vec = [-cog5_S_x    -cog5_S_y    zeros(size(phi3))];
R_cog6_vec = [-cog6_R_x    -cog6_R_y    zeros(size(phi3))];
U_cog6_vec = [-cog6_U_x    -cog6_U_y    zeros(size(phi3))];
S_cog7_vec = [-cog7_S_x    -cog7_S_y    zeros(size(phi3))];
T_cog7_vec = [-cog7_T_x    -cog7_T_y    zeros(size(phi3))];
T_cog8_vec = [-cog8_T_x    -cog8_T_y    zeros(size(phi3))];
U_cog8_vec = [-cog8_U_x    -cog8_U_y    zeros(size(phi3))];

PQ_vec = [L3*cos(phi3) L3*sin(phi3) zeros(size(phi3))];
PX_vec = [l3*cos(theta5) l3*sin(theta5) zeros(size(phi3))];
VU_vec = [L4*cos(phi4) L4*sin(phi4) zeros(size(phi4))];
UT_vec = [L3*cos(phi3) L3*sin(phi3)   zeros(size(phi3))];

% acceleration vectors
acc_3 =       cross(omega3,cross(omega3,P_cog3_vec))+cross(alpha3,P_cog3_vec);
acc_Q =       cross(omega3,cross(omega3,PQ_vec    ))+cross(alpha3,PQ_vec    );
acc_5 = acc_Q+cross(omega5,cross(omega5,Q_cog5_vec))+cross(alpha5,Q_cog5_vec);
acc_X =       cross(omega3,cross(omega3,PX_vec    ))+cross(alpha3,PX_vec    );
acc_2 = acc_X+cross(omega3,cross(omega3,PX_vec    ))+cross(alpha3,PX_vec    );
acc_4 =       cross(omega4,cross(omega4,V_cog4_vec))+cross(alpha4,V_cog4_vec);
acc_U =       cross(omega4,cross(omega4,VU_vec    ))+cross(alpha4,VU_vec    );
acc_6 = acc_U+cross(omega6,cross(omega6,U_cog6_vec))+cross(alpha6,U_cog6_vec);
acc_8 = acc_U+cross(omega8,cross(omega8,U_cog8_vec))+cross(alpha8,U_cog8_vec);
acc_T = acc_U+cross(omega8,cross(omega8,UT_vec    ))+cross(alpha8,UT_vec    );
acc_7 = acc_T+cross(omega7,cross(omega7,T_cog7_vec))+cross(alpha7,T_cog7_vec);


% **********************
% *** force analysis ***
% **********************

% allocate matrices for force (F) and moment (M)
F_P_x = zeros(size(phi3));
F_P_y = zeros(size(phi3));
F_Q_x = zeros(size(phi3));
F_Q_y = zeros(size(phi3));
F_R_x = zeros(size(phi3));
F_R_y = zeros(size(phi3));
F_S_x = zeros(size(phi3));
F_S_y = zeros(size(phi3));
F_T_x = zeros(size(phi3));
F_T_x = zeros(size(phi3));
F_U_x = zeros(size(phi3));
F_U_x = zeros(size(phi3));
F_V_x = zeros(size(phi3));
F_V_x = zeros(size(phi3));
F_W_x = zeros(size(phi3));
F_W_x = zeros(size(phi3));
F_X_x = zeros(size(phi3));
F_X_y = zeros(size(phi3));

M_P = zeros(size(phi3));

% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
  A = [ 1           0            1            0            0            0            0           0           0;
        0           1            0            1            0            0            0           0           0;
        0           0           -1            0           -1            0            0           0           0;
        0           0            0           -1            0           -1            0           0           0;
        0           0            0            0            1            0            1           0           0;
        0           0            0            0            0            1            0           1           0;
       -cog2_P_y(k) cog2_P_x(k) -cog2_Q_y(k)  cog2_Q_x(k)  0            0            0           0           1;
        0           0            cog3_Q_y(k) -cog3_Q_x(k)  cog3_R_y(k) -cog3_R_x(k)  0           0           0;
        0           0            0            0           -cog4_R_y(k)  cog4_R_x(k) -cog4_S_y(k) cog4_S_x(k) 0];
    
  B = [ m2*acc_2x(k);
        m2*acc_2y(k);
        m3*acc_3x(k);
        m3*acc_3y(k);
        m4*acc_4x(k);
        m4*acc_4y(k);
        J2*ddphi2(k);
        J3*ddphi3(k);
        J4*ddphi4(k)];
    
    x = A\B;
    
    % save results
    F_P_x(k) = x(1);
    F_P_y(k) = x(2);
    F_Q_x(k) = x(3);
    F_Q_y(k) = x(4);
    F_R_x(k) = x(5);
    F_R_y(k) = x(6);
    F_S_x(k) = x(7);
    F_S_y(k) = x(8);
    M_P(k)   = x(9);
end



% **********************
% *** plot figures ***
% **********************

if fig_dyn_4bar
    
    figure
    subplot(221)
    plot(F_P_x,F_P_y),grid
    xlabel('F_P_x [N]')
    ylabel('F_P_y [N]')
    axis tight
    subplot(222)
    plot(F_Q_x,F_Q_y),grid
    xlabel('F_Q_x [N]')
    ylabel('F_Q_y [N]')
    axis tight
    subplot(223)
    plot(F_R_x,F_R_y),grid
    xlabel('F_R_x [N]')
    ylabel('F_R_y [N]')
    axis tight
    subplot(224)
    plot(F_S_x,F_S_y),grid
    xlabel('F_S_x [N]')
    ylabel('F_S_y [N]')
    axis tight
    
    figure
    plot(t,M_P)
    ylabel('M_P [N-m]')
    xlabel('t [s]')
    
end


