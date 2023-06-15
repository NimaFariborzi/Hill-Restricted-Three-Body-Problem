clc; clear; close all

% Parameters
G = 6.674e-11;
m1 = 5.972e24; % Earth mass, kg
m2 = 7.348e22; % Moon mass, kg
R = 385000e3;  % Mean distance to moon in m
v = 1.022e3;   % Mean moon speed in m/s
omega = v/R;   % Mean moon angular velocity about Earth in rad/s

% Roughly circular orbit about moon %
r_dim = 4000e3;           % Orbit radius in m (4,000 km)
v_dim = sqrt(G*m2/r_dim); % Circular orbit velocity in m/s
r0 = [r_dim; 0; 0];
v0 = [0; v_dim; 0];
t_end = 27.4*3600*24;

% Integrate
x0 = [r0; v0];
tspan = [0 t_end];
options = odeset('AbsTol',1e-14,'RelTol',3e-14);
[t,x] = ode113(@(t,x) HR3BP_EOM(t,x,G*m1,G*m2,R,omega),tspan,x0,options);

% Get results for plotting
rx = x(:,1);
ry = x(:,2);
rz = x(:,3);
vx = x(:,4);
vy = x(:,5);
vz = x(:,6);

% Calculate Jacobi integral
nt = length(t);
J = zeros(1,nt);
for i = 1:nt
    r = norm([rx(i), ry(i), rz(i)]);
    J(i) = 0.5*(vx(i)^2 + vy(i)^2 + vz(i)^2) - G*m2/r - 3/2*G*m1/R^3 * rx(i)^2 + 1/2*G*m1/R^3*rz(i)^2;
end

% Plot trajectory in rotating frame
figure()
view(3)
hold on
sphere3(0,0,0,1738e3,[0.7 0.7 0.7]) % Plot moon
plot3(rx,ry,rz,'k-')
plot3(rx(end),ry(end),rz(end),'ko','MarkerFaceColor','k')
axis equal
xlabel('$x$ (m)')
ylabel('$y$ (m)')
zlabel('$z$ (m)')
title('HR3BP Trajectory','Moon-Centered Rotating Frame')
legend('Moon','Spacecraft Trajectory','Location','northeast')

% Plot Jacobi integral
figure()
plot(t,J/J(1))
xlabel('$t$ (s)')
ylabel('$J/J_0$')
title('(Normalized) Jacobi Integral vs Time')

% Transform coordinates to inertial frame
[x_inert, x_inert_sec] = rotational_to_inertial(x,t,R,omega);

% Plot trajectory in inertial coordinates
figure()
view(3)
hold on
sphere3(0,0,0,6378e3,[0.011719 0.66016 0.98438]) % Plot earth
plot3(x_inert_sec(:,1),x_inert_sec(:,2),x_inert_sec(:,3),'Color',[0.7 0.7 0.7]) % Secondary trajectory
plot3(x_inert(:,1),x_inert(:,2),x_inert(:,3),'k-')
xlabel('$x$ (m)')
ylabel('$y$ (m)')
zlabel('$z$ (m)')
title('HR3BP Trajectory','Earth-Centered Inertial Frame')
axis equal
legend('Earth','Moon Trajectory','Spacecraft Trajectory')

% Interesting Orbit
r0 = [-3e7; 0; 5e7];
v0 = [0; -50; 0];
t_end = 27.4*3600*24;

% Integrate
x0 = [r0; v0];
tspan = [0 t_end];
options = odeset('AbsTol',1e-14,'RelTol',3e-14);
[t,x] = ode113(@(t,x) HR3BP_EOM(t,x,G*m1,G*m2,R,omega),tspan,x0,options);

% Get results for plotting
rx = x(:,1);
ry = x(:,2);
rz = x(:,3);
vx = x(:,4);
vy = x(:,5);
vz = x(:,6);

% Calculate Jacobi integral
nt = length(t);
J = zeros(1,nt);
for i = 1:nt
    r = norm([rx(i), ry(i), rz(i)]);
    J(i) = 0.5*(vx(i)^2 + vy(i)^2 + vz(i)^2) - G*m2/r - 3/2*G*m1/R^3 * rx(i)^2 + 1/2*G*m1/R^3*rz(i)^2;
end

% Plot trajectory in rotating frame
figure()
view(3)
hold on
sphere3(0,0,0,1738e3,[0.7 0.7 0.7]) % Plot moon
plot3(rx,ry,rz,'k-')
plot3(rx(end),ry(end),rz(end),'ko','MarkerFaceColor','k')
axis equal
xlabel('$x$ (m)')
ylabel('$y$ (m)')
zlabel('$z$ (m)')
title('HR3BP Trajectory','Moon-Centered Rotating Frame')
legend('Moon','Spacecraft Trajectory','Location','northeast')

% Plot Jacobi integral
figure()
plot(t,J/J(1))
xlabel('$t$ (s)')
ylabel('$J/J_0$')
title('(Normalized) Jacobi Integral vs Time')

% Transform coordinates to inertial frame
[x_inert, x_inert_sec] = rotational_to_inertial(x,t,R,omega);

% Plot trajectory in inertial coordinates
figure()
view(3)
hold on
sphere3(0,0,0,6378e3,[0.011719 0.66016 0.98438]) % Plot earth
plot3(x_inert_sec(:,1),x_inert_sec(:,2),x_inert_sec(:,3),'Color',[0.7 0.7 0.7]) % Secondary trajectory
plot3(x_inert(:,1),x_inert(:,2),x_inert(:,3),'k-')
xlabel('$x$ (m)')
ylabel('$y$ (m)')
zlabel('$z$ (m)')
title('HR3BP Trajectory','Earth-Centered Inertial Frame')
axis equal
legend('Earth','Moon Trajectory','Spacecraft Trajectory')

% Equations of motion for HR3BP in dimensional coordinates
function Xdot = HR3BP_EOM(~,X,Gm1,Gm2,R,omega)
    % Get magnitude of r vector
    r3 = norm(X(1:3))^3;
    R3 = R^3;
    % Initialize solution array
    Xdot = [X(4:6); 0; 0; 0];
    % Calculate
    Xdot(4) =  2*omega*X(5) + 3*Gm1/R3*X(1) - Gm2/r3*X(1);
    Xdot(5) = -2*omega*X(4)                 - Gm2/r3*X(2);
    Xdot(6) =                 - Gm1/R3*X(3) - Gm2/r3*X(3);
end

% Function to draw sphere with given radius and position
function sphere3(x,y,z,r,col)
    [X,Y,Z] = sphere(50);
    X2 = X * r;
    Y2 = Y * r;
    Z2 = Z * r;
    surf(X2+x,Y2+y,Z2+z,'EdgeColor','none','FaceColor',col);
end