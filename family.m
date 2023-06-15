clc; clear; close all

% Initialize IC arrays
xi    = {};
J     = {};
t_end = {};

xi{end+1}    = 0.56291;
J{end+1}     = 4;
t_end{end+1} = 1.63594 * 2;

xi{end+1}    = -0.56291;
J{end+1}     = 4;
t_end{end+1} = 1.63594 * 2;

xi{end+1}    = 0.58201;
J{end+1}     = 4.1;
t_end{end+1} = 1.57627 * 2;

xi{end+1}    = -0.58201;
J{end+1}     = 4.1;
t_end{end+1} = 1.57627 * 2;

xi{end+1}    = 0.60040;
J{end+1}     = 4.2;
t_end{end+1} = 1.45398 * 2;

xi{end+1}    = -0.60040;
J{end+1}     = 4.2;
t_end{end+1} = 1.45398 * 2;

xi{end+1}    = 0.55029;
J{end+1}     = 4.3;
t_end{end+1} = 0.96514 * 2;

xi{end+1}    = -0.55029;
J{end+1}     = 4.3;
t_end{end+1} = 0.96514 * 2;

xi{end+1}    = 0.43684;
J{end+1}     = 4.4;
t_end{end+1} = 0.71494 * 2;

xi{end+1}    = -0.43684;
J{end+1}     = 4.4;
t_end{end+1} = 0.71494 * 2;

xi{end+1}    = .28350;
J{end+1}     = 4.49999;
t_end{end+1} = 0.61294 * 2;

% For 2D orbits, convert (xi, J) to initial position, velocity
n = length(xi);
r0    = cell(1,n);
v0    = cell(1,n);
x0    = cell(1,n);
tspan = cell(1,n);
for i = 1:length(xi)
    r0{i} = [xi{i}; 0; 0];
    eta_dot = sqrt(3*xi{i}^2 + 2/abs(xi{i}) - J{i});
    if xi{i} < 0
        eta_dot = eta_dot * -1;
    end
    v0{i} = [0; eta_dot; 0];
    x0{i} = [r0{i}; v0{i}];
    tspan{i} = [0 t_end{i}];
end

% Integrate
t = cell(1,n);
x = cell(1,n);
options = odeset('AbsTol',1e-12,'RelTol',1e-12);
for i = 1:n
    [t_out,x_out] = ode113(@(t,x) HR3BP_Dimless_EOM(t,x,'forward'),tspan{i},x0{i},options);
    t{i} = t_out;
    x{i} = x_out;
end

% Plot trajectory in rotating frame
figure()
view(2)
hold on
cellfun(@(x) PlotTraj(x,'k'), x)
axis equal
xlabel('$\rho_x$')
ylabel('$\rho_y$')
title('HR3BP Trajectories','Family of Stable Periodic Orbits')
pause(0.1)
% exportgraphics(gcf,'Family.png','Resolution',300)

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

% Function to draw sphere with given radius and position
function sphere3(x,y,z,r,col)
    [X,Y,Z] = sphere(50);
    X2 = X * r;
    Y2 = Y * r;
    Z2 = Z * r;
    surf(X2+x,Y2+y,Z2+z,'EdgeColor','none','FaceColor',col);
end

function PlotTraj(x,col)
    plot3(x(:,1),x(:,2),x(:,3),'Color',col)
end