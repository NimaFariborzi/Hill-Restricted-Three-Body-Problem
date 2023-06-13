close all
clearvars
clc

%Initial conditions(all planar)

%X01=[(1/3)^(1/3)-0.03; 1.e-11;0;0]; %IC around L2

xi1 = .28350;
J1 = 4.49999;
eta_dot1 = sqrt(3*xi1^2 + 2/xi1 - J1);

xi2 = .56291;
J2 = 4.300000;
eta_dot2 = sqrt(3*xi2^2 + 2/xi2 - J2);

X01=[xi1;0;0;eta_dot1];
X02=[xi2;0;0;eta_dot2];

options = odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',@xcross_event); 
tspan = [0,1500]; 

[T1,X1,Tevent1,Xevent1] = ode45(@(t,x) eom_hR3bp_2d(x),tspan,X01,options); 
[T2,X2,Tevent2,Xevent2] = ode45(@(t,x) eom_hR3bp_2d(x),tspan,X02,options); 

x1=Xevent1(:,1);
vx1=Xevent1(:,3);
x2=Xevent2(:,1);
vx2=Xevent2(:,3);

figure(1)
plot(x1,vx1,'k.',x1,-vx1,'k.','MarkerSize',15)
title('Poincare Surface of Section(y=0,ydot<0)')
xlabel('x')
ylabel('xdot')
grid on;


figure(2)
plot(x2,vx2,'r.',x2,-vx2,'r.','MarkerSize',15)
title('Poincare Surface of Section(y=0,ydot<0)')
xlabel('x')
ylabel('xdot')
grid on;

%% functions

function sys_of_1st_ordereqns=eom_hR3bp_3d(X)
x=X(1);
y=X(2);
z=X(3);
dxdt=X(4); %vx
dydt=X(5); %vy
dzdt=X(6); %vz
r=(x^2+y^2+z^2)^0.5;
dvxdt=2*dydt-(x/(r^3))+3*x; %ax
dvydt=-2*dxdt-(y/(r^3)); %ay
dvzdt=-(z/(r^3))-z; %az
sys_of_1st_ordereqns=[dxdt;dydt;dzdt;dvxdt;dvydt;dvzdt];
end

function sys_of_1st_ordereqns=eom_hR3bp_2d(X)
x=X(1);
y=X(2);
dxdt=X(3); %vx
dydt=X(4); %vy
r=(x^2+y^2)^0.5;
dvxdt=2*dydt-(x/(r^3))+3*x; %ax
dvydt=-2*dxdt-(y/(r^3)); %ay
sys_of_1st_ordereqns=[dxdt;dydt;dvxdt;dvydt];
end

function [value, isterminal, direction]=xcross_event(~,X) %event when crosses poincare' section
value=X(2);
isterminal=0;
direction=-1;
end