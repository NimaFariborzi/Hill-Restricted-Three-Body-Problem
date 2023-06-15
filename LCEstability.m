close all
clearvars
clc

pertur=0.01;

%planar periodic IC
x01 = .28350; 
J1 = 4.49999;
ydot1 = sqrt(3*x01^2 + 2/x01 - J1);
X01=[x01;0;0;ydot1];

%perturbed IC
X01p=[x01+pertur;0;0;ydot1];

options = odeset('AbsTol',1e-14,'RelTol',1e-14); 
tspan = [0,200]; 

[T1,X1] = ode113(@(t,x) eom_hR3bp_2d(x),tspan,X01,options);
[T1p,X1p] = ode113(@(t,x) eom_hR3bp_2d(x),tspan,X01p,options);

X1pinter=interp1(T1p,X1p,T1,'spline');
deltaX=zeros(size(T1,1),1);
r=zeros(size(T1,1),1);
LCE=zeros(size(T1,1),1);

%loop for calculating the LCE at each time step
for i=1:size(T1,1)
    deltaX(i)= norm(X1(i,:)-X1pinter(i,:));
    r(i)=log(deltaX(i)/deltaX(1));
    LCE(i)= 1/T1(i) *r(i);
    i
end

%plot
figure(3)
plot(T1(2:end),LCE(2:end))
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

