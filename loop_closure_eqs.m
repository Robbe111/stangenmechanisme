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


function F=loop_closure_eqs(phi_init,phi3,L1,L2,L3,L4,l3,l4,L5,l5,L6,L7,l7,L8,phi1)

% first argument: the initial values of the unknown angles phi3 and phi4
% argument phi2: input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
% arguments a1 ... phi1: constants


% copy initial values of unknown angles phi3 and phi4
phi2=phi_init(1);
phi4=phi_init(2);
phi5=phi_init(3);
phi6=phi_init(4);
phi7=phi_init(5);
phi8=phi_init(6);

theta3 = phi3 + deg2rad(80);
theta4 = phi4 - deg2rad(100);
% loop closure equations:
F(1)= L3*cos(phi3) + l5*cos(phi5) - L6*cos(phi6) - L4*cos(phi4) - L1;
F(2)=L3*sin(phi3) + L5*sin(phi5) - L6*sin(phi6) - L4*sin(phi4);
F(3)=l3*cos(theta3) + L2 * cos(phi2) - l4*cos(theta4) - L1;
F(4)=l3 *sin(theta3) + L2* sin(phi2) - l4*sin(theta4);
F(5)=L3*cos(phi3) + L5*cos(phi5) - l7*cos(phi7) - L8*cos(phi8) - L4*cos(phi4) - L1;
F(6)=L3*sin(phi3) + L5*sin(phi5) - l7*sin(phi7) - L8*sin(phi8) - L4*sin(phi4);
end 

