function [x, y, v, strikes] = mover(x, y, v, npart, L, h, mpv_exc, mpv_unp, tau)

% mover - Function to move particles by free flight.
%         Also handles collisions with walls.
%
% Inputs
%       x           x-positions of the particles
%       y           y-positions of the particles
%       v       Velocities of the particles
%       npart       Number of particles in the system
%       L           System length
%       h           System height
%       mpv_exc     Most probable velocity off the wall in excited plasma
%       mpv_unp     Most probably velocity off the wall in unperturbed plasma
%       tau         Timestep
%
% Outputs
%       x, y, v     Updated positions and velocities
%       strikes     Number of particles striking each wall
%       delv        Change of y-velocity at each wall

% Move all particles pretending walls are absent
% Split according to temp
x_old = x;          % Remember original x-positions
x(:) = x_old(:) + v(:,1) * tau;
y_old = y;          % Remember original y-positions
y(:) = y_old(:) + v(:,2) * tau;

% Loop over all particles
strikes = [0 0];
xwall = [0 L];
ywall = [0 h];
direction = [1 -1];         % Direction of particle leaving wall
stdev_exc = mpv_exc/sqrt(2);
stdev_unp = mpv_unp/sqrt(2);

for i = 1:npart
    
    % Test if particle strikes either wall
    if( x(i) <= 0 )
        flag = 1;           % Particle strikes left wall
    elseif( x(i) >= L )
        flag = 2;           % Particle strikes right wall
    else
        flag = 0;           % Particle strikes neither wall
    end
    
    % If particle strikes a wall, reset its position and velocity.
    % Record velocity change.
    if( flag > 0 )
        strikes(flag) = strikes(flag) + 1;
        vyInitial = v(i,2);
        % Reset velocity components as biased Maxwellian,
        % Exponential dist. in x; Gaussian in y and z
        % Split according to local temp
        if( flag == 1 )
            v(i, 1) = direction(flag) * abs(sqrt(-log(abs(1-rand(1))))) * mpv_exc;
            v(i, 2) = stdev_exc * randn(1);
            v(1, 3) = stdev_exc * randn(1);
        else
            v(i, 1) = direction(flag) * abs(sqrt(-log(abs(1-rand(1))))) * mpv_unp;
            v(i, 2) = stdev_unp * randn(1);
            v(1, 3) = stdev_unp * randn(1);
        end
        
        % Time of flight after leaving wall
        dtr = tau * (x(i) - xwall(flag))/(x(i)-x_old(i));
        % Reset position after leaving wall
        x(i) = xwall(flag) + v(i,1)*dtr;
    end
    
    % Now for y-direction
    
    % Test if particle strikes either wall
    if( y(i) <= 0 )
        flag = 1;           % Particle strikes bottom wall
    elseif( y(i) >= h )
        flag = 2;           % Particle strikes top wall
    else
        flag = 0;           % Particle strikes neither wall
    end
    
    % If particle strikes a wall, reset its position and velocity. 
    % Record velocity change.
    if( flag > 0 )
        strikes(flag) = strikes(flag) + 1;
        % Reset velocity components as biased Maxwellian.
        % Exponential dist. in y; Gaussian in x and z
        % Split according to local temp
        if( x(i) < L/2 )
            v(i, 1) = stdev_exc * randn(1);
            v(i, 2) = direction(flag) * abs(sqrt(-log(abs(1-randn(1))))) * mpv_exc;
            v(i, 3) = stdev_exc * randn(1);
        else
            v(i, 1) = stdev_unp * randn(1);
            v(i, 2) = direction(flag) * abs(sqrt(-log(abs(1-randn(1))))) * mpv_unp;
            v(i, 3) = stdev_unp * randn(1);
        end
       
        % Time of flight after leaving wall
        dtr = tau * (y(i) - ywall(flag))/(y(i) - y_old(i));
        % Reset position after leaving wall
        y(i) = ywall(flag) + v(i, 2)*dtr;
    end
end

return;
        
        