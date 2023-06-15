close all
clearvars
clc

pertur=0.0001;

%planar periodic initial condition 
x01 = .28350; 
J1 = 4.49999;
ydot1 = sqrt(3*x01^2 + 2/x01 - J1);
X01=[x01;0;0;ydot1];

x01=0.2835;
xdot1=0;
y01=0;
ydot1=1.6721;
J1=0.5*(xdot1^2+ydot1^2)-(1/(sqrt(x01^2+y01^2)))-0.5*(3*x01^2);

% %perturbed IC
% x01p = .28350+pertur; 
% J1p = 4.49999;
% ydot1p = sqrt(3*x01p^2 + 2/x01p - J1p);
% X01p=[x01p;0;0;ydot1p];


options = odeset('AbsTol',1e-12,'RelTol',1e-12); 
tspan = [0,10]; 

[T1,X1] = ode45(@(t,x) eom_hR3bp_2d(x),tspan,X01,options);
[T1p,X1p] = ode45(@(t,x) eom_hR3bp_2d(x),tspan,X01p,options);

X1pinter=interp1(T1p,X1p,T1,'spline');
deltaX=zeros(size(T1,1),1);
r=zeros(size(T1,1),1);
LE=zeros(size(T1,1),1);
for i=1:size(T1,1)
    deltaX(i)= sqrt((X1(i,1)-X1pinter(i,1))^2+(X1(i,2)-X1pinter(i,2))^2);
    %deltaX(i)= sqrt((X1(i,1)-X1p(i,1))^2+(X1(i,2)-X1p(i,2))^2);
    r(i)=log(deltaX(i)/deltaX(1));
    LE(i)= 1/size(T1,1) *sum(r,"all");
    i
end




figure(3)
plot(T1,LE)
%legend()
xlabel('time');
ylabel('LCE');
title('LCE vs time');





%% functions
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

