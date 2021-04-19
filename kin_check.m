% Emma Wessel& Robbe Vermeiren: kinematic check
function [ derphi2,derphi3,derphi4,derphi5,derphi6,derphi7,derphi8,...
           dderphi2,dderphi3,dderphi4,dderphi5,dderphi6,dderphi7,dderphi8 ] = ...
kin_check( phi2,phi3,phi4,phi5,phi6,phi7,phi8, ...
           dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8, ...
           ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8, ...
           Ts, fig_kin_check, t)
       



% f'(x) = (f(x+1)-f(x-1))/(2*Ts)

derphi2 = (phi2(3:size(phi2)) - phi2(1:size(phi2)-2)) / (2*Ts);
derphi3  = (phi3(3:size(phi3)) - phi3(1:size(phi3)-2)) / (2*Ts);
derphi4  = (phi4(3:size(phi3)) - phi4(1:size(phi3)-2)) / (2*Ts);
derphi5  = (phi5(3:size(phi3)) - phi5(1:size(phi3)-2)) / (2*Ts);
derphi6  = (phi6(3:size(phi3)) - phi6(1:size(phi3)-2)) / (2*Ts);
derphi7  = (phi7(3:size(phi3)) - phi7(1:size(phi3)-2)) / (2*Ts);
derphi8  = (phi8(3:size(phi3)) - phi8(1:size(phi3)-2)) / (2*Ts);


% f"(x) = (f(x+1)+f(x-1)-2*f(x))/(Ts^2)

dderphi2  = (phi2(3:size(phi2)) + phi2(1:size(phi2)-2) - 2*phi2(2:size(phi2)-1)) / (Ts^2);
dderphi3  = (phi3(3:size(phi3)) + phi3(1:size(phi3)-2) - 2*phi3(2:size(phi3)-1)) / (Ts^2);
dderphi4  = (phi4(3:size(phi3)) + phi4(1:size(phi3)-2) - 2*phi4(2:size(phi3)-1)) / (Ts^2);
dderphi5  = (phi5(3:size(phi3)) + phi5(1:size(phi3)-2) - 2*phi5(2:size(phi3)-1)) / (Ts^2);
dderphi6  = (phi6(3:size(phi3)) + phi6(1:size(phi3)-2) - 2*phi6(2:size(phi3)-1)) / (Ts^2);
dderphi7  = (phi7(3:size(phi3)) + phi7(1:size(phi3)-2) - 2*phi7(2:size(phi3)-1)) / (Ts^2);
dderphi8  = (phi8(3:size(phi3)) + phi8(1:size(phi3)-2) - 2*phi8(2:size(phi3)-1)) / (Ts^2);

% errors

errordphi2  = derphi2 - dphi2(2:size(phi3)-1,:);
errordphi3  = derphi3 - dphi3(2:size(phi3)-1,:);
errordphi4  = derphi4 - dphi4(2:size(phi3)-1,:);
errordphi5  = derphi5 - dphi5(2:size(phi3)-1,:);
errordphi6  = derphi6 - dphi6(2:size(phi3)-1,:);
errordphi7  = derphi7 - dphi7(2:size(phi3)-1,:);
errordphi8  = derphi8 - dphi8(2:size(phi3)-1,:);

errorddphi2  = dderphi2 - ddphi2(2:size(phi3)-1,:);
errorddphi3  = dderphi3 - ddphi3(2:size(phi3)-1,:);
errorddphi4  = dderphi4 - ddphi4(2:size(phi3)-1,:);
errorddphi5  = dderphi5 - ddphi5(2:size(phi3)-1,:);
errorddphi6  = dderphi6 - ddphi6(2:size(phi3)-1,:);
errorddphi7  = dderphi7 - ddphi7(2:size(phi3)-1,:);
errorddphi8  = dderphi8 - ddphi8(2:size(phi3)-1,:);


% plot if fig_kin_check = 1

if fig_kin_check
figure
subplot(311)
    plot(t(1:size(errordphi2),:),derphi2)
    ylabel('d\phi_2 approximation [rad/s]')
    xlabel('t [s]')
subplot(312)
    plot(t,dphi2)
    ylabel('d\phi_2 exact [rad/s]')
    xlabel('t [s]')
subplot(313)
    plot(t(1:size(errordphi2),:),errordphi2)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad/s]')
    xlabel('t [s]')

    
figure

subplot(331)
    plot(t(1:size(errordphi3),:),derphi3)
    ylabel('d\phi_3 approximation [rad/s]')
    xlabel('t [s]')
subplot(332)
    plot(t,dphi3)
    ylabel('d\phi_3 exact [rad/s]')
    xlabel('t [s]')
subplot(333)
    plot(t(1:size(errordphi3),:),errordphi3)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad/s]')
    xlabel('t [s]')
subplot(334)
    plot(t(1:size(errordphi4),:),derphi4)
    ylabel('d\phi_4 approximation [rad/s]')
    xlabel('t [s]')
subplot(335)
    plot(t,dphi4)
    ylabel('d\phi_4 exact [rad/s]')
    xlabel('t [s]')
subplot(336)
    plot(t(1:size(errordphi4),:),errordphi4)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad/s]')
    xlabel('t [s]')
subplot(337)
    plot(t(1:size(errordphi5),:),derphi5)
    ylabel('d\phi_5 approximation [rad/s]')
    xlabel('t [s]')
subplot(338)
    plot(t,dphi5)
    ylabel('d\phi_5 exact [rad/s]')
    xlabel('t [s]')
subplot(339)
    plot(t(1:size(errordphi5),:),errordphi5)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad/s]')
    xlabel('t [s]')
    
figure

subplot(331)
    plot(t(1:size(errordphi6),:),derphi6)
    ylabel('d\phi_6 approximation [rad/s]')
    xlabel('t [s]')
subplot(332)
    plot(t,dphi6)
    ylabel('d\phi_6 exact [rad/s]')
    xlabel('t [s]')
subplot(333)
    plot(t(1:size(errordphi6),:),errordphi6)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad/s]')
    xlabel('t [s]')
subplot(334)
    plot(t(1:size(errordphi7),:),derphi7)
    ylabel('d\phi_7 approximation [rad/s]')
    xlabel('t [s]')
subplot(335)
    plot(t,dphi7)
    ylabel('d\phi_7 exact [rad/s]')
    xlabel('t [s]')
subplot(336)
    plot(t(1:size(errordphi7),:),errordphi7)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad/s]')
    xlabel('t [s]')
subplot(337)
    plot(t(1:size(errordphi8),:),derphi8)
    ylabel('d\phi_8 approximation [rad/s]')
    xlabel('t [s]')
subplot(338)
    plot(t,dphi8)
    ylabel('d\phi_8 exact [rad/s]')
    xlabel('t [s]')
subplot(339)
    plot(t(1:size(errordphi8),:),errordphi8)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad/s]')
    xlabel('t [s]')
    
figure

subplot(331)
    plot(t(1:size(errorddphi3),:),dderphi3)
    ylabel('dd\phi_3 approximation [rad/s^2]')
    xlabel('t [s]')
subplot(332)
    plot(t,ddphi3)
    ylabel('dd\phi_3 exact [rad/s^2]')
    xlabel('t [s]')
subplot(333)
    plot(t(1:size(errorddphi3),:),errorddphi3)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad/s^2]')
    xlabel('t [s]')
subplot(334)
    plot(t(1:size(errorddphi4),:),dderphi4)
    ylabel('dd\phi_4 approximation [rad/s^2]')
    xlabel('t [s]')
subplot(335)
    plot(t,ddphi4)
    ylabel('dd\phi_4 exact [rad/s^2]')
    xlabel('t [s]')
subplot(336)
    plot(t(1:size(errorddphi4),:),errorddphi4)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad/s^2]')
    xlabel('t [s]')
subplot(337)
    plot(t(1:size(errorddphi5),:),dderphi5)
    ylabel('dd\phi_5 approximation [rad/s^2]')
    xlabel('t [s]')
subplot(338)
    plot(t,ddphi5)
    ylabel('dd\phi_5 exact [rad/s^2]')
    xlabel('t [s]')
subplot(339)
    plot(t(1:size(errorddphi5),:),errorddphi5)
    axis([0 10 -0.3 0.3])
    ylabel('error [rad/s^2]')
    xlabel('t [s]')
    
figure

subplot(331)
    plot(t(1:size(errorddphi6),:),dderphi6)
    ylabel('dd\phi_6 approximation [rad/s^2]')
    xlabel('t [s]')
subplot(332)
    plot(t,ddphi6)
    ylabel('dd\phi_6 exact [rad/s^2]')
    xlabel('t [s]')
subplot(333)
    plot(t(1:size(errorddphi6),:),errorddphi6)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad/s^2]')
    xlabel('t [s]')
subplot(334)
    plot(t(1:size(errorddphi7),:),dderphi7)
    ylabel('dd\phi_7 approximation [rad/s^2]')
    xlabel('t [s]')
subplot(335)
    plot(t,ddphi7)
    ylabel('dd\phi_7 exact [rad/s^2]')
    xlabel('t [s]')
subplot(336)
    plot(t(1:size(errorddphi7),:),errorddphi7)
    axis([0 10 -0.3 0.3])
    ylabel('error [rad/s^2]')
    xlabel('t [s]')
subplot(337)
    plot(t(1:size(errorddphi8),:),dderphi8)
    ylabel('dd\phi_8 approximation [rad/s^2]')
    xlabel('t [s]')
subplot(338)
    plot(t,ddphi8)
    ylabel('dd\phi_8 exact [rad/s^2]')
    xlabel('t [s]')
subplot(339)
    plot(t(1:size(errorddphi8),:),errorddphi8)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad/s^2]')
    xlabel('t [s]')
  

end

