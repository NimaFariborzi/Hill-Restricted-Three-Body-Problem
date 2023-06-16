function [X_inertial, X_inertial_sec] = rotational_to_inertial(X,t,R,omega)
%ROTATIONAL_TO_INERTIAL(X,t,R,omega) transforms trajectories in a rotating frame
%about the secondary body to inertial coordinates about the primary body.
%   - X is the dimensional trajectory in a rotating frame centered on the
%       secondary body.
%   - t is the time vector.
%   - R is the distance between the primary and secondary body.
%   - omega is the angular velocity of the secondary body about the primary.
%Coordinates are transformed to an inertial frame assuming the secondary body
%starts on the x-axis to the right of the primary, and orbits counterclockwise
%in the x-y plane.
    nt = length(t);
    V = R*omega; % Velocity of secondary about primary
    X_inertial     = zeros(size(X));
    X_inertial_sec = zeros(size(X));
    for i = 1:nt
        % Get rotation angle from x-axis
        theta = omega*t(i);
        % Calculate position and velocity of secondary in inertial coordinates
        X_sec = [R* cos(theta), R* sin(theta), 0 ...
                 V*-sin(theta), V* cos(theta), 0];
        % Calculate spacecraft position and velocity about the secondary,
        % rotated by angle theta
        Q = [cos(theta), -sin(theta);
             sin(theta),  cos(theta)]; % Rotation matrix
        xy_rot = Q*X(i,1:2)';
        uv_rot = Q*X(i,4:5)';
        % Position: Secondary position plus rotated position relative to
        % secondary
        X_inertial(i,1:3) = X_sec(1:3) + [xy_rot', X(i,3)];
        % Velocity: Secondary position plus rotated velocity relative to
        % secondary
        X_inertial(i,4:6) = X_sec(4:6) + [uv_rot', X(i,6)];
        % Secondary position and velocity
        X_inertial_sec(i,:) = X_sec;
    end
end