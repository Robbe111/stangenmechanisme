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


function [F_P_x, F_P_y, F_X_x, F_X_y, F_Q_x ,F_Q_y, F_W_x, F_W_y, F_V_x, F_V_y ,F_U4_x, F_U4_y, F_U6_x, F_U6_y, F_U8_x, F_U8_y, F_T_x ,F_T_y ,F_R_x, F_R_y, F_S_x, F_S_y, M_P, ...
                vel_3x, vel_3y, vel_5x, vel_5y, vel_2x, vel_2y, vel_4x, vel_4y, vel_6x, vel_6y, vel_8x, vel_8y, vel_7x, vel_7y, ...
                acc_3x, acc_3y, acc_5x, acc_5y, acc_2x, acc_2y, acc_4x, acc_4y, acc_6x, acc_6y, acc_8x, acc_8y, acc_7x, acc_7y, ...
                omega2, omega3, omega4, omega5, omega6, omega7, omega8, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8] = ...
                dynamics_4bar(phi2,phi3,phi4,phi5,phi6,phi7,phi8,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8 ...
                                                                                                    ,L1,L2,L3,L4,l3,l4,L5,l5,L6,L7,l7,L8, ...
                                                                                                        m2,m3,m4,m5,m6,m7,m8,J2,J3,J4,J5,J6,J7,J8,t,fig_dyn_4bar)

theta3 = phi3 + deg2rad(80);
theta4 = phi4 - deg2rad(100);
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
cog5_Q_x = -(L5/2)*cos(phi5);
cog5_Q_y = -(L5/2)*sin(phi5);
cog5_R_x = (l5-L5/2)*cos(phi5);
cog5_R_y = (l5-L5/2)*sin(phi5);
cog5_S_x = -cog5_Q_x;
cog5_S_y = -cog5_Q_y;
cog6_R_x = (L6/2)*cos(phi6);
cog6_R_y = (L6/2)*sin(phi6);
cog6_U_x = -(L6/2)*cos(phi6);
cog6_U_y = -(L6/2)*sin(phi6);
cog7_S_x = l7/2*cos(phi7);
cog7_S_y = l7/2*sin(phi7);
cog7_T_x = -l7/2*cos(phi7);
cog7_T_y = -l7/2*sin(phi7);
cog8_U_x = -L8/2*cos(phi8);
cog8_U_y = -L8/2*sin(phi8);
cog8_T_x = L8/2*cos(phi8);
cog8_T_y = L8/2*sin(phi8);
cog4_U_x = L4*cos(phi4)-cog4_x;
cog4_U_y = L4*sin(phi4)-cog4_y; 
cog4_V_x = -cog4_x;  
cog4_V_y = -cog4_y; 
cog4_W_x = l4*cos(theta4) - cog4_x;
cog4_W_y = l4*sin(theta4) - cog4_y;
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
PX_vec = [l3*cos(theta3) l3*sin(theta3) zeros(size(phi3))];
VU_vec = [L4*cos(phi4) L4*sin(phi4) zeros(size(phi4))];
UT_vec = [L8*cos(phi8) L8*sin(phi8)   zeros(size(phi3))];


% velocity vectors
vel_3 = cross(omega3,P_cog3_vec);
vel_Q = cross(omega3,PQ_vec);
vel_5 = vel_Q + cross(omega5,Q_cog5_vec);
vel_X = cross(omega3,PX_vec);
vel_2 = vel_X + cross(omega2,X_cog2_vec);
vel_4 = cross(omega4,V_cog4_vec);
vel_U = cross(omega4,VU_vec);
vel_6 = vel_U + cross(omega6,U_cog6_vec);
vel_8 = vel_U + cross(omega8,U_cog8_vec);
vel_T = vel_U + cross(omega8,UT_vec);
vel_7 = vel_T + cross(omega7,T_cog7_vec);
vel_3x = vel_3(:,1);
vel_3y = vel_3(:,2);
vel_5x = vel_5(:,1);
vel_5y = vel_5(:,2);
vel_2x = vel_2(:,1);
vel_2y = vel_2(:,2);
vel_4x = vel_4(:,1);
vel_4y = vel_4(:,2);
vel_6x = vel_6(:,1);
vel_6y = vel_6(:,2);
vel_8x = vel_8(:,1);
vel_8y = vel_8(:,2);
vel_7x = vel_7(:,1);
vel_7y = vel_7(:,2);



% acceleration vectors
acc_3 =       cross(omega3,cross(omega3,P_cog3_vec))+cross(alpha3,P_cog3_vec);
acc_Q =       cross(omega3,cross(omega3,PQ_vec    ))+cross(alpha3,PQ_vec    );
acc_5 = acc_Q+cross(omega5,cross(omega5,Q_cog5_vec))+cross(alpha5,Q_cog5_vec);
acc_X =       cross(omega3,cross(omega3,PX_vec    ))+cross(alpha3,PX_vec    );
acc_2 = acc_X+cross(omega2,cross(omega2,X_cog2_vec))+cross(alpha2,X_cog2_vec);
acc_4 =       cross(omega4,cross(omega4,V_cog4_vec))+cross(alpha4,V_cog4_vec);
acc_U =       cross(omega4,cross(omega4,VU_vec    ))+cross(alpha4,VU_vec    );
acc_6 = acc_U+cross(omega6,cross(omega6,U_cog6_vec))+cross(alpha6,U_cog6_vec);
acc_8 = acc_U+cross(omega8,cross(omega8,U_cog8_vec))+cross(alpha8,U_cog8_vec);
acc_T = acc_U+cross(omega8,cross(omega8,UT_vec    ))+cross(alpha8,UT_vec    );
acc_7 = acc_T+cross(omega7,cross(omega7,T_cog7_vec))+cross(alpha7,T_cog7_vec);
acc_3x = acc_3(:,1);
acc_3y = acc_3(:,2);
acc_5x = acc_5(:,1);
acc_5y = acc_5(:,2);
acc_2x = acc_2(:,1);
acc_2y = acc_2(:,2); 
acc_4x = acc_4(:,1);
acc_4y = acc_4(:,2);
acc_6x = acc_6(:,1);
acc_6y = acc_6(:,2);
acc_8x = acc_8(:,1);
acc_8y = acc_8(:,2);
acc_7x = acc_7(:,1);
acc_7y = acc_7(:,2);


% **********************
% *** force analysis ***
% **********************
% x = [F13x F13y F32x F32y F35x F35y F24x F24y F14x F14y F46x F46y F68x
         % F68y F84x F84y F87x F87y F56x F56y F57x F57y M]
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
F_T_y = zeros(size(phi3));
F_U4_x = zeros(size(phi3));
F_U4_y = zeros(size(phi3));
F_U6_x = zeros(size(phi3));
F_U6_y = zeros(size(phi3));
F_U8_x = zeros(size(phi3));
F_U8_y = zeros(size(phi3));
F_V_x = zeros(size(phi3));
F_V_y = zeros(size(phi3));
F_W_x = zeros(size(phi3));
F_W_y = zeros(size(phi3));
F_X_x = zeros(size(phi3));
F_X_y = zeros(size(phi3));
M_P = zeros(size(phi3));

% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps¨
% x = [F13x F13y F32x F32y F35x F35y F24x F24y F14x F14y F46x F46y F68x
         % F68y F84x F84y F87x F87y F56x F56y F57x F57y M]
% 
for k=1:t_size
        % F_P_x(k)   F_P_y(k)    F_X_x(k)    F_X_y(k)    F_Q_x(k)     F_Q_y(k)    F_W_x(k)    F_W_y(k)     F_V_x(k)   F_V_y(k)   F_U4_x(k)    F_U4_y(k)    F_U6_x(k)  F_U6_y(k)   F_U8_x(k)   F_U8_y(k)   F_T_x(k)    F_T_y(k)    F_R_x(k)    F_R_y(k)    F_S_x(k)    F_S_y(k) 
    
    A = [ 0          0           -1            0            0            0            1           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0   0;%stang 2 
        0           0            0           -1            0            0            0           1           0           0           0           0           0           0           0           0           0           0           0           0           0           0   0;
        0           0            cog2_X_y(k)  -cog2_X_x(k) 0            0           -cog2_W_y(k) cog2_W_x(k) 0           0           0           0           0           0           0           0           0           0           0           0           0           0   0;
        
        1           0            1            0            1            0            0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0   0;%lichaam 3
        0           1            0            1            0            1            0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0   0;
       -cog3_P_y(k) cog3_P_x(k) -cog3_X_y(k)  cog3_X_x(k) -cog3_Q_y(k)  cog3_Q_x(k)  0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0   1;
        
        0           0            0            0            0            0           -1           0           1            0            1           0           0           0           0           0           0           0           0           0           0           0   0;%lichaam 4
        0           0            0            0            0            0            0          -1           0            1            0           1           0           0           0           0           0           0           0           0           0           0   0;
        0           0            0            0            0            0           cog4_W_y(k) -cog4_W_x(k) -cog4_V_y(k) cog4_V_x(k) -cog4_U_y(k) cog4_U_x(k) 0           0           0           0           0           0           0           0           0           0   0; 
        
        0           0            0            0            0            0            0           0           0           0           0           0           0           0           1           0           -1            0           0           0           0         0   0;%lichaam 8
        0           0            0            0            0            0            0           0           0           0           0           0           0           0           0           1            0           -1           0           0           0         0   0;
        0           0            0            0            0            0            0           0           0           0           0           0           0           0          -cog8_U_y(k) cog8_U_x(k)  cog8_T_y(k) -cog8_T_x(k) 0           0           0         0   0;
        
        0           0            0            0           -1            0            0           0           0           0           0           0           0           0           0           0           0           0           1           0            1           0           0;%lichaam 5
        0           0            0            0            0           -1            0           0           0           0           0           0           0           0           0           0           0           0           0           1            0           1           0;
        0           0            0            0            cog5_Q_y(k) -cog5_Q_x(k)  0           0           0           0           0           0           0           0           0           0           0           0          -cog5_R_y(k) cog5_R_x(k) -cog5_S_y(k) cog5_S_x(k) 0;
       
        0           0            0            0            0            0            0           0           0           0           0           0           1           0           0           0           0           0          -1            0           0           0   0;%lichaam 6
        0           0            0            0            0            0            0           0           0           0           0           0           0           1           0           0           0           0           0           -1           0           0   0;
        0           0            0            0            0            0            0           0           0           0           0           0          -cog6_U_y(k) cog6_U_x(k) 0           0           0           0           cog6_R_y(k) -cog6_R_x(k) 0           0   0;
       
        0           0            0            0            0            0            0           0           0           0           0           0           0           0           0           0           1           0           0           0          -1            0           0;%lichaam 7
        0           0            0            0            0            0            0           0           0           0           0           0           0           0           0           0           0           1           0           0           0           -1           0;
        0           0            0            0            0            0            0           0           0           0           0           0           0           0           0           0          -cog7_T_y(k) cog7_T_x(k) 0           0           cog7_S_y(k) -cog7_S_x(k) 0;   
       
        0           0            0            0            0            0            0           0           0           0           1           0           1           0           1           0           0           0           0           0           0           0              0;                        
        0           0            0            0            0            0            0           0           0           0           0           1           0           1           0           1           0           0           0           0           0           0              0 ];
  B = [ m2*acc_2x(k);
        m2*acc_2y(k);
        J2*ddphi2(k);
        m3*acc_3x(k);
        m3*acc_3y(k);
        J3*ddphi3(k);
        m4*acc_4x(k);
        m4*acc_4y(k);
        J4*ddphi4(k);
        m8*acc_8x(k);
        m8*acc_8y(k);
        J8*ddphi8(k);
        m5*acc_5x(k);
        m5*acc_5y(k);
        J5*ddphi5(k);
        m6*acc_6x(k);
        m6*acc_6y(k);
        J6*ddphi6(k);
        m7*acc_7x(k);
        m7*acc_7y(k);
        J7*ddphi7(k);
        0;
        0];
    
 
    
    x = A\B;
    % x = [F13x F13y F32x F32y F35x F35y F24x F24y F14x F14y F46x F46y F68x
         % F68y F84x F84y F87x F87y F56x F56y F57x F57y M]
% 
    % save results
    F_P_x(k) = x(1);
    F_P_y(k) = x(2);
    F_X_x(k) = x(3);
    F_X_y(k) = x(4);
    F_Q_x(k) = x(5);
    F_Q_y(k) = x(6);
    F_W_x(k) = x(7);
    F_W_y(k) = x(8);
    F_V_x(k) = x(9);
    F_V_y(k) = x(10);
    F_U4_x(k) = x(11);
    F_U4_y(k) = x(12);
    F_U6_x(k) = x(13);
    F_U6_y(k) = x(14);
    F_U8_x(k) = x(15);
    F_U8_y(k) = x(16);
    F_T_x(k) = x(17);
    F_T_y(k) = x(18);
    F_R_x(k) = x(19);
    F_R_y(k) = x(20);
    F_S_x(k) = x(21);
    F_S_y(k) = x(22);
    M_P(k) = x(23);


 
    
    
    
    
    
    
    
%     F_P_x(k) = x(1);
%     F_P_y(k) = x(2);
%     F_Q_x(k) = x(3);
%     F_Q_y(k) = x(4);
%     F_R_x(k) = x(5);
%     F_R_y(k) = x(6);
%     F_S_x(k) = x(7);
%     F_S_y(k) = x(8);
%     M_P(k)   = x(9);
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
    plot(F_X_x,F_X_y),grid
    xlabel('F_X_x [N]')
    ylabel('F_X_y [N]')
    axis tight
    subplot(223)
    plot(F_Q_x,F_Q_y),grid
    xlabel('F_Q_x [N]')
    ylabel('F_Q_y [N]')
    axis tight
    
    figure
    subplot(221)
    plot(F_R_x,F_R_y),grid
    xlabel('F_R_x [N]')
    ylabel('F_R_y [N]')
    axis tight
    subplot(222)
    plot(F_S_x,F_S_y),grid
    xlabel('F_S_x [N]')
    ylabel('F_S_y [N]')
    axis tight
    subplot(223)
    plot(F_T_x,F_T_y),grid
    xlabel('F_T_x [N]')
    ylabel('F_T_y [N]')
    axis tight
    
   
    figure
    subplot(311)
    plot(t,F_P_x,'r')
    xlabel('t[s]')
    ylabel('F_P_x[N]')
    
    subplot(312)
    plot(t,F_X_x,'r')
    xlabel('t[s]')
    ylabel('F_X_x[N]')
    
    subplot(313)
    plot(t,F_Q_x,'r')
    xlabel('t[s]')
    ylabel('F_Q_x[N]')
    
    figure 
    subplot(311)
    plot(t,F_P_y,'r')
    xlabel('t[s]')
    ylabel('F_P_y[N]')
    
    subplot(312)
    plot(t,F_X_y,'r')
    xlabel('t[s]')
    ylabel('F_X_y[N]')
    
    subplot(313)
    plot(t,F_Q_y,'r')
    xlabel('t[s]')
    ylabel('F_Q_y[N]')
    
    
    
    
    
    figure
    subplot(311)
    plot(t,F_R_x,'r')
    xlabel('t[s]')
    ylabel('F_R_x[N]')
    
    subplot(312)
    plot(t,F_S_x,'r')
    xlabel('t[s]')
    ylabel('F_S_x[N]')
    
    subplot(313)
    plot(t,F_T_x)
    xlabel('t[s]')
    ylabel('F_T_x[N]')    
    
    figure 
    subplot(311)
    plot(t,F_R_y,'r')
    xlabel('t[s]')
    ylabel('F_R_y[N]')
    
    subplot(312)
    plot(t,F_S_y,'r')
    xlabel('t[s]')
    ylabel('F_S_y[N]')
    
    subplot(313)
    plot(t,F_T_y)
    xlabel('t[s]')
    ylabel('F_T_y[N]')  
    
    figure
    subplot(311)
    plot(t,F_U4_x,'r')
    xlabel('t[s]')
    ylabel('F_U4_x[N]')
    
    subplot(312)
    plot(t,F_U6_x,'r')
    xlabel('t[s]')
    ylabel('F_U6_x[N]')
    
    subplot(313)
    plot(t,F_U8_x)
    xlabel('t[s]')
    ylabel('F_U8_x[N]')    
    
    figure 
    subplot(311)
    plot(t,F_U4_y,'r')
    xlabel('t[s]')
    ylabel('F_U4_y[N]')
    
    subplot(312)
    plot(t,F_U6_y,'r')
    xlabel('t[s]')
    ylabel('F_U6_y[N]')
    
    subplot(313)
    plot(t,F_U8_y)
    xlabel('t[s]')
    ylabel('F_U8_y[N]')
    
    figure
    subplot(221)
    plot(t,F_V_x)
    xlabel('t[s]')
    ylabel('F_V_x[N]')
    
    subplot(222)
    plot(t,F_V_y)
    xlabel('t[s]')
    ylabel('F_V_y[N]')
    
    
    subplot(223)
    plot(t,F_W_x)
    xlabel('t[s]')
    ylabel('F_W_x[N]')
    
    
    subplot(224)
    plot(t,F_W_y)
    xlabel('t[s]')
    ylabel('F_W_y[N]')
    
    
    
    
    
    
    figure
    plot(t,M_P)
    ylabel('M_P [N-m]')
    xlabel('t [s]')
    
end


