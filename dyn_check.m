%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Analyse 8 stangenmechanisme
% Wynn Geenen <wynn.geenen@student.kuleuven.be>
% Jan Fissette <janfissette@student.kuleuven.be>
%
% Naar een voorbeeld van:
%
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dyn_check(vel_2x,vel_2y,vel_3x,vel_3y,vel_4x,vel_4y,vel_5x,vel_5y,vel_6x,vel_6y,vel_7x,vel_7y,vel_8x,vel_8y,...
                   acc_2x,acc_2y,acc_3x,acc_3y,acc_4x,acc_4y,acc_5x,acc_5y,acc_6x,acc_6y,acc_7x,acc_7y,acc_8x,acc_8y,...
                   M_P,omega2,omega3,omega4,omega5,omega6,omega7,omega8,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,m2,m3,m4,m5,m6,m7,m8,...
                   J2,J3,J4,J5,J6,J7,J8,t,fig_dyn_check,g)

dE2_pot = m2*g*vel_2y;
dE3_pot = m3*g*vel_3y;
dE4_pot = m4*g*vel_4y;
dE5_pot = m5*g*vel_5y;
dE6_pot = m6*g*vel_6y;
dE7_pot = m7*g*vel_7y;
dE8_pot = m8*g*vel_8y;

% variation in kinetic energy of each bar
dE_kin2 = m2*(vel_2x.*acc_2x+vel_2y.*acc_2y) + J2*omega2(:,3).*alpha2(:,3);
dE_kin3 = m3*(vel_3x.*acc_3x+vel_3y.*acc_3y) + J3*omega3(:,3).*alpha3(:,3);
dE_kin4 = m4*(vel_4x.*acc_4x+vel_4y.*acc_4y) + J4*omega4(:,3).*alpha4(:,3);
dE_kin5 = m5*(vel_5x.*acc_5x+vel_5y.*acc_5y) + J5*omega5(:,3).*alpha5(:,3);
dE_kin6 = m6*(vel_6x.*acc_6x+vel_6y.*acc_6y) + J6*omega6(:,3).*alpha6(:,3);
dE_kin7 = m7*(vel_7x.*acc_7x+vel_7y.*acc_7y) + J7*omega7(:,3).*alpha7(:,3);
dE_kin8 = m8*(vel_8x.*acc_8x+vel_8y.*acc_8y) + J8*omega8(:,3).*alpha8(:,3);


% total kinetic energy
dEtot = dE_kin2+dE_kin3+dE_kin4+dE_kin5+dE_kin6+dE_kin7+dE_kin8+dE2_pot+dE3_pot+dE4_pot+dE5_pot+dE6_pot+dE7_pot+dE8_pot;


Mcontr = dEtot./omega3(:,3);

% error between driving force of kin_8bar and dyn_check
Err = M_P - Mcontr;


%plot if fig_dyn_check = 1
if fig_dyn_check
    
  figure
  subplot(311)
  plot(t,M_P),grid
  xlabel('t [s]')        
  ylabel('M_P [N]')
  subplot(312)
  plot(t,Mcontr),grid
  xlabel('t [s]')
  ylabel('Mcontr [N]')
  subplot(313)
  plot(t,Err),grid
  xlabel('t [s]')
  ylabel('Error [N]')
  
end
