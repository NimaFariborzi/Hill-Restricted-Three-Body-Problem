clc; clear; close all

load('v_u.mat')
load('v_s.mat')

x_u = v_u(1:3);
v_u = v_u(4:6);
x_s = v_s(1:3);
v_s = v_s(4:6);

perturb = 4.2e-1; % Factor for perturbation

% Initial Conditions
% L2 Equilibrium Point
x_eq = [(1/3)^(1/3); 0; 0];
v_eq = [0; 0; 0];
% Add perturbations
v0_u1 = v_eq + perturb*v_u;
v0_u2 = v_eq - perturb*v_u;
v0_s1 = v_eq + perturb*v_s;
v0_s2 = v_eq - perturb*v_s;

% Time to integrate
t_end = 6;

% Integrate
X0_u1 = [x_eq; v0_u1];
X0_u2 = [x_eq; v0_u2];
X0_s1 = [x_eq; v0_s1];
X0_s2 = [x_eq; v0_s2];
tspan = [0 t_end];
options = odeset('AbsTol',1e-12,'RelTol',1e-9);
[~,x_u1] = ode45(@(t,x) HR3BP_Dimless_EOM(t,x,'forward'),tspan,X0_u1,options);
[~,x_u2] = ode45(@(t,x) HR3BP_Dimless_EOM(t,x,'forward'),tspan,X0_u2,options);
[~,x_s1] = ode45(@(t,x) HR3BP_Dimless_EOM(t,x,'backward'),tspan,X0_s1,options);
[~,x_s2] = ode45(@(t,x) HR3BP_Dimless_EOM(t,x,'backward'),tspan,X0_s2,options);

%solving jacobi to implement zero velocity curves
r=sqrt(x_u1(:,1).^2+x_u1(:,2).^2+x_u1(:,3).^2);
J=1/2.*(x_u1(:,4).^2+x_u1(:,5).^2+x_u1(:,6).^2)-1./r-1/2.*(3*x_u1(:,1).^2-x_u1(:,3).^2);

% %plotting zero velocity curves 3D
% f3p=fimplicit3(@(x3fb,y3fb,z3fb)J(1)+1/(sqrt(x3fb^2+y3fb^2+z3fb^2))+3/2*x3fb^2);


%plotting zero velo curves 2D
fp=fimplicit(@(xfb,yfb)J(1)+1/(sqrt(xfb^2+yfb^2))+3/2*xfb^2);
hold on

%extracting coords to fill area of zero velocity
x_upper=fp.XData(fp.YData>0);
x_lower=fp.XData(fp.YData<0);
y_upper=fp.YData(fp.YData>0);
y_lower=fp.YData(fp.YData<0);

patch(x_upper,y_upper,'k')
patch(x_lower,y_lower,'k')

% Plot trajectory in rotating frame
view(2)
hold on
plot(0,0,'ro','MarkerFaceColor','r')
plot_traj(x_u1,'k')
plot_traj(x_s1,[.7 .7 .7])
plot_traj(x_u2,'k')
plot_traj(x_s2,[.7 .7 .7])
plot( x_eq(1),x_eq(2),'kx','LineWidth',1,'MarkerSize',20)
plot(-x_eq(1),x_eq(2),'kx','LineWidth',1,'MarkerSize',20)
text( x_eq(1)-0.02,x_eq(2)+0.04,'L2','FontSize',18)
text(-x_eq(1)+0.02,x_eq(2),'L1','FontSize',18)
xlim([-1 1])
ylim([-1 1])
% zlim([-1 1])
axis equal
xlabel('$x$ (dimensionless)')
ylabel('$y$ (dimensionless)')
zlabel('$z$ (dimensionless)')
title('Hill Restricted 3-Body Problem','Stable and Unstable Manifolds About L2')
legend('zero velocity curves','upper forbidden region','lower forbidden region','Secondary','Unstable Manifolds','Stable Manifolds')
% exportgraphics(gcf,'Manifolds.png','Resolution',300)

%visualzing flow stuck inside zero velocity curve
%im guessing here
x0s=[.3,.2,0];
v0s=[.3,-.04,0];
X0s = [x0s; v0s];
tsend=10;
tspan = [0 tsend];
%solving for traj
[~,xs] = ode45(@(t,x) HR3BP_Dimless_EOM(t,x,'forward'),tspan,X0s,options);
%solving for J
rs=sqrt(xs(:,1).^2+xs(:,2).^2+xs(:,3).^2);
Js=1/2.*(xs(:,4).^2+xs(:,5).^2+xs(:,6).^2)-1./rs-1/2.*(3*xs(:,1).^2-xs(:,3).^2);
%plotting zero velocity curves
figure()
fps=fimplicit(@(xfb,yfb)Js(1)+1/(sqrt(xfb^2+yfb^2))+3/2*xfb^2);
hold on

xs_upper=fps.XData(fps.YData>0);
xs_lower=fps.XData(fps.YData<0);
ys_upper=fps.YData(fps.YData>0);
ys_lower=fps.YData(fps.YData<0);

% patch(xh_upper,yh_upper,'k')
% patch(xh_lower,yh_lower,'k')

% Plot trajectory in rotating frame
view(2)
hold on
plot(0,0,'ro','MarkerFaceColor','r')
plot_traj(xs,'k')
plot( x_eq(1),x_eq(2),'kx','LineWidth',1,'MarkerSize',15)
plot(-x_eq(1),x_eq(2),'kx','LineWidth',1,'MarkerSize',15)
text( x_eq(1)-0.02,x_eq(2)+0.04,'L2','FontSize',10)
text(-x_eq(1)+0.02,x_eq(2),'L1','FontSize',10)
xlim([-1 1])
ylim([-1 1])
% zlim([-1 1])
axis equal
xlabel('$x$ (dimensionless)')
ylabel('$y$ (dimensionless)')
zlabel('$z$ (dimensionless)')
title('Hill Restricted 3-Body Problem','forbidden region covering exits')
legend('zero velocity curves','upper forbidden region','lower forbidden region','Secondary','trajectory')




%visualzing flow hitting the zero velocity curve
x0h=[.5,0,0];
v0h=[.3,0,0];
X0h = [x0h; v0h];
thend=6.5;
tspan = [0 thend];
%solving for traj
[~,xh] = ode45(@(t,x) HR3BP_Dimless_EOM(t,x,'forward'),tspan,X0h,options);
%solving for J
rh=sqrt(xh(:,1).^2+xh(:,2).^2+xh(:,3).^2);
Jh=1/2.*(xh(:,4).^2+xh(:,5).^2+xh(:,6).^2)-1./rh-1/2.*(3*xh(:,1).^2-xh(:,3).^2);
%plotting zero velocity curves
figure()
fph=fimplicit(@(xfb,yfb)Jh(1)+1/(sqrt(xfb^2+yfb^2))+3/2*xfb^2);
hold on

xh_upper=fph.XData(fph.YData>0);
xh_lower=fph.XData(fph.YData<0);
yh_upper=fph.YData(fph.YData>0);
yh_lower=fph.YData(fph.YData<0);

patch(xh_upper,yh_upper,'k')
patch(xh_lower,yh_lower,'k')

% Plot trajectory in rotating frame
view(2)
hold on
plot(0,0,'ro','MarkerFaceColor','r')
plot_traj(xh,'k')
% plot(x0h(1),x0h(2),'bo')
plot( x_eq(1),x_eq(2),'kx','LineWidth',1,'MarkerSize',15)
plot(-x_eq(1),x_eq(2),'kx','LineWidth',1,'MarkerSize',15)
text( x_eq(1)-0.02,x_eq(2)+0.04,'L2','FontSize',10)
text(-x_eq(1)+0.02,x_eq(2),'L1','FontSize',10)
xlim([-1 1])
ylim([-1 1])
% zlim([-1 1])
axis equal
xlabel('$x$ (dimensionless)')
ylabel('$y$ (dimensionless)')
zlabel('$z$ (dimensionless)')
title('Hill Restricted 3-Body Problem','spacecraft reflecting off forbidden region')
legend('zero velocity curves','upper forbidden region','lower forbidden region','Secondary','trajectory')






%solving Traj of 3d orbits with zero velo curves
%im guessing here
x03=[.5,.3,.4];
v03=[0,.4,.2];
X03 = [x03; v03];
t3end=5;
tspan = [0 t3end];
%solving for traj
[~,x3] = ode45(@(t,x) HR3BP_Dimless_EOM(t,x,'forward'),tspan,X0_u1,options);
%solving for J
r3=sqrt(x3(:,1).^2+x3(:,2).^2+x3(:,3).^2);
J3=1/2.*(x3(:,4).^2+x3(:,5).^2+x3(:,6).^2)-1./r3-1/2.*(3*x3(:,1).^2-x3(:,3).^2);
%plotting zero velocity curves 3D
figure()
%I'm trying to make the picture better, hard.... the interval fixes the
%implicit but i need the traj to be fixed as well
% interval=[-1 1 -1 1 -1 1];
f3p=fimplicit3(@(x3fb,y3fb,z3fb)J3(1)+1/(sqrt(x3fb^2+y3fb^2+z3fb^2))+3/2*x3fb^2-1/2*z3fb^2,'EdgeColor','none','FaceAlpha',0.25);
hold on
plot3(0,0,0,'ro','MarkerFaceColor','k')
plot_traj3(x3,'r')
% plot_traj3(x_s1,[.7 .7 .7])
% plot_traj3(x_u2,'b')
% plot_traj3(x_s2,[.7 .7 .7])
% plot3( x_eq(1),x_eq(2),0,'kx','LineWidth',3,'MarkerSize',20)
% plot3(-x_eq(1),x_eq(2),0,'kx','LineWidth',3,'MarkerSize',20)
% text( x_eq(1)-0.02,x_eq(2)+0.04,'L2','FontSize',18)
% text(-x_eq(1)+0.02,x_eq(2),'L1','FontSize',18)
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
axis equal
xlabel('$x$ (dimensionless)')
ylabel('$y$ (dimensionless)')
zlabel('$z$ (dimensionless)')
title('Hill Restricted 3-Body Problem','Zero Velocity Planes')
legend('zero velocity planes','Secondary','trajectory')


function Xdot = HR3BP_Dimless_EOM(~,X,state)
    % Get magnitude of r vector
    r3 = norm(X(1:3))^3;
    % Initialize solution array
    Xdot = [X(4:6); 0; 0; 0];
    % Calculate
    Xdot(4) =  2*X(5) - X(1)/r3 + 3*X(1);
    Xdot(5) = -2*X(4) - X(2)/r3;
    Xdot(6) =         - X(3)/r3 - X(3);
    if strcmp(state,'backward')
        Xdot = Xdot*-1; % Propagate backwards
    end
end

function plot_traj(x,col)
    rx = x(:,1);
    ry = x(:,2);
    plot(rx,ry,'Color',col)
end

function plot_traj3(x,col)
    rx = x(:,1);
    ry = x(:,2);
    rz = x(:,3);
    plot3(rx,ry,rz,'Color',col)
end
