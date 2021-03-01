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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phi2,phi4,phi5,phi6,phi7,phi8,dphi2,dphi4,dphi5,dphi6,dphi7,dphi8,ddphi2,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8] = ... 
                                    kinematics_4bar(Ts,L1,L2,L3,L4,l3,l4,L5,l5,L6,L7,l7,L8,...
                                            phi2_init,phi4_init,phi5_init,phi6_init,phi7_init,phi8_init,phi3,dphi3,ddphi3,t,fig_kin_4bar,phi1);


% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
% angles
theta3 = zeros(size(t));
theta4 = zeros(size(t)); 
phi2 = zeros(size(t));
phi4 = zeros(size(t));
phi5 = zeros(size(t));
phi6 = zeros(size(t));
phi7 = zeros(size(t));
phi8 = zeros(size(t));
%velocities
dphi2 = zeros(size(t));
dphi4 = zeros(size(t));
dphi5 = zeros(size(t));
dphi6 = zeros(size(t));
dphi7 = zeros(size(t));
dphi8 = zeros(size(t));
%accelerations
ddphi2 = zeros(size(t));
ddphi4 = zeros(size(t));
ddphi5 = zeros(size(t));
ddphi6 = zeros(size(t));
ddphi7 = zeros(size(t));
ddphi8 = zeros(size(t));


% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off');

% *** loop over positions ***
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
    % *** position analysis ***
    
    % fsolve solves the non-linear set of equations
    % loop closure equations: see loop_closure_eqs.m
    % argument loop_closure_eqs: file containing closure equations
    % argument [..]': initial values of unknown angles phi3 and phi4
    % argument optim options: parameters for fsolve
    % argument phi2(k): input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
    % argument a1 ... phi1: constants
    % return value x: solution for the unknown angles phi3 and phi4
    % return exitflag: indicates convergence of algorithm
    [x, fval, exitflag]=fsolve('loop_closure_eqs',[phi2_init phi4_init phi5_init phi6_init phi7_init phi8_init]',optim_options,phi3(k),L1,L2,L3,L4,l3,l4,L5,l5,L6,L7,l7,L8,phi1);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    
    % save results of fsolve
    phi2(k)=x(1);
    phi4(k)=x(2);
    phi5(k)=x(3);
    phi6(k)=x(4);
    phi7(k)=x(5);
    phi8(k)=x(6);
    theta3(k) = phi3(k) + deg2rad(80);
    theta4(k) = phi4(k) - deg2rad(100);
    % *** velocity analysis ***
    
    A = [0,  -L5*sin(phi5(k)), 0, l7*sin(phi7(k)), L8*sin(phi8(k)), L4*sin(phi4(k));
         0, L5*cos(phi5(k)), 0, -l7*cos(phi7(k)), -L8*cos(phi8(k)), -L4*cos(phi4(k));
         0, -l5*sin(phi5(k)), L6*sin(phi6(k)), 0,       0,        L4*sin(phi4(k));
         0, l5*cos(phi5(k)), -L6*cos(phi6(k)), 0,       0,        -L4*cos(phi4(k));
         -L2*sin(phi2(k)),0,                0,  0,       0,        l4*sin(phi4(k)); 
         L2*cos(phi2(k)),0,                0,  0,       0,        -l4*cos(phi4(k)); 
         ];
    B = [ L3*sin(phi3(k))*dphi3(k);
          -L3*cos(phi3(k))*dphi3(k);
          L3*sin(phi3(k))*dphi3(k);
          -L3*cos(phi3(k))*dphi3(k);
          l3*sin(theta3(k))*dphi3(k);
          -l3*cos(theta3(k))*dphi3(k)
        ];
     
    x = A\B;
    
    % save results
    dphi2(k)=x(1);
    dphi5(k)=x(2);
    dphi6(k)=x(3);
    dphi7(k)=x(4);
    dphi8(k)=x(5);
    dphi4(k)=x(6);
   
    
%     % *** acceleration analysis ***
%     
%     A = [-r3*sin(phi3(k)),  r4*sin(phi4(k));
%          r3*cos(phi3(k)), -r4*cos(phi4(k))];
%     B = [r2*cos(phi2(k))*dphi2(k)^2+r2*sin(phi2(k))*ddphi2(k)+r3*cos(phi3(k))*dphi3(k)^2-r4*cos(phi4(k))*dphi4(k)^2;
%          r2*sin(phi2(k))*dphi2(k)^2-r2*cos(phi2(k))*ddphi2(k)+r3*sin(phi3(k))*dphi3(k)^2-r4*sin(phi4(k))*dphi4(k)^2];
%     
%     x = A\B;
%     % save results
%     ddphi3(k) = x(1);
%     ddphi4(k) = x(2);
    
    
    % *** calculate initial values for next iteration step ***
    phi2_init = phi2(k)+Ts*dphi2(k);
    phi4_init = phi4(k)+Ts*dphi4(k);
    phi5_init = phi5(k)+Ts*dphi5(k);
    phi6_init = phi6(k)+Ts*dphi6(k);
    phi7_init = phi7(k)+Ts*dphi7(k);
    phi8_init = phi8(k)+Ts*dphi8(k);

    
end % loop over positions



% *** create movie ***
% 
% % point P = fixed
% P = 0;
% % point S = fixed
% S = r1*exp(j*phi1);
% % define which positions we want as frames in our movie
% frames = 40;    % number of frames in movie
% delta = floor(t_size/frames); % time between frames
% index_vec = [1:delta:t_size]';
% 
% % Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% % This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% % axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% % plots.
% x_left = -1.5*r2;
% y_bottom = -1.5*max(r2,r4);
% x_right = r1+1.5*r4;
% y_top = 1.5*max(r2,r4);
% 
% figure(10)
% hold on
% plot([x_left, x_right], [y_bottom, y_top]);
% axis equal;
% movie_axes = axis;   %save current axes into movie_axes
% 
% % draw and save movie frame
% for m=1:length(index_vec)
%     index = index_vec(m);
%     Q = P + r2 * exp(j*phi2(index));
%     R1 = Q + r3 * exp(j*phi3(index));
%     R2 = S + r4 * exp(j*phi4(index));
%     
%     loop1 = [P Q R1 R2 S];
%     
%     figure(10)
%     clf
%     hold on
%     plot(real(loop1),imag(loop1),'-o')
%     
%     axis(movie_axes);     % set axes as in movie_axes
%     Movie(m) = getframe;  % save frame to a variable Film
% end
% 
% % save movie
% save fourbar_movie Movie
% close(10)
% 
% 
% % *** plot figures ***
% 
% if fig_kin_4bar
%     
%     %plot assembly at a certain timestep 
%     index = 1; %select 1st timestep
%     P = 0;
%     S = r1*exp(j*phi1);
%     Q = P + r2 * exp(j*phi2(index));
%     R = Q + r3 * exp(j*phi3(index));
%     
%     figure
%     assembly=[P, Q, R, S];
%     plot(real(assembly),imag(assembly),'ro-')
%     xlabel('[m]')
%     ylabel('[m]')
%     title('assembly')
%     axis equal
%     
%     figure
%     subplot(311)
%     plot(t,phi2)
%     ylabel('\phi_2 [rad]')
%     subplot(312)
%     plot(t,phi3)
%     ylabel('\phi_3 [rad]')
%     subplot(313)
%     plot(t,phi4)
%     ylabel('\phi_4 [rad]')
%     xlabel('t [s]')
%     
%     figure
%     subplot(311)
%     plot(t,dphi2)
%     ylabel('d\phi_2 [rad/s]')
%     subplot(312)
%     plot(t,dphi3)
%     ylabel('d\phi_3 [rad/s]')
%     subplot(313)
%     plot(t,dphi4)
%     ylabel('d\phi_4 [rad/s]')
%     xlabel('t [s]')
%     
%     figure
%     subplot(311)
%     plot(t,ddphi2)
%     ylabel('dd\phi_2 [rad/s^2]')
%     subplot(312)
%     plot(t,ddphi3)
%     ylabel('dd\phi_3 [rad/s^2]')
%     subplot(313)
%     plot(t,ddphi4)
%     ylabel('dd\phi_4 [rad/s^2]')
%     xlabel('t [s]')
% end



