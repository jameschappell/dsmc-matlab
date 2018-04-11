% dscmne - program to simulate a dilute gas using DSMC algorithm
clear all;
help dsmcne;    % clear memory and print header

% Initialise constants (particle mass, diameter, etc.)
boltz = 1.3806e-23;         % Boltzmann's constant (J/K)
mass = 1.44e-25;            % Mass of Rb atom (kg)
diam = 6.06e-10;            % Effective diameter of Rb atom (m) from http://en.wikipedia.org/wiki/Rubidium
T = 473;                    % Initial temperature (K)
density = 7.0e14/(0.01)^3;  % Number density of Rb (m^-3)
L = 2.0;                    % System size (m)
h = 0.04;                   % System height (m)

Volume = h*L*L;             % Volume of the system (m^3)
npart = input('Enter number of simulation particles: ');
eff_num = density * Volume / npart;
fprintf('Each simulation particle represents %g atoms\n', eff_num);
mfp = 1/(sqrt(2) * pi * diam^2 * density);      % mean free path
fprintf('Mean free path is %g mm\n', mfp * 1000.0);

fprintf('System length %.1f m is %g mean free paths\n', L, L / mfp);
mpv = sqrt(2 * boltz * T / mass);       % most probable initial velocity
fprintf('Most probable initial velocity = %g m/s\n', mpv);

% Assign random positions and velocities to particles
rand('state', 1);           % initialise random number generators
randn('state', 1);

% Creating initial gas density distribution

L_plasma = 0.5 * L;

dl_bin = 2e-2;
lbin = dl_bin/2:dl_bin:L;   % bins for x-histogram

x = L_plasma * rand(npart,1);   % assign random x-positions
y = h * rand(npart,1);          % assign random y-positions

% Assign thermal velocities using Gaussian random numbers
v = sqrt(boltz * T / mass) * randn(npart,3);

% Initialise variables used for evaluating collisions
dl = 0.5e-3;                % cells are rectangular: dl in x and y
ncell = (L/dl) * (h/dl);    % number of cells
tau = 0.2 * dl / mpv;       % set timestep tau
fprintf('Number of cells = (%.1f)x(%.1f) = %.1f\n', L/dl, h/dl, ncell);
fprintf('%.1f particles per cell\n', npart/ncell);
fprintf('Timestep tau = %g microseconds\n', tau * 1e6);

vrmax = 3 * mpv * ones(ncell, 1);   % Estimated max relative speed in a cell
selxtra = zeros(ncell, 1);          % Used by collision routine "collider"
coeff = 0.5 * eff_num * pi * diam^2 * tau / (Volume / ncell);

% Declare structure for lists used in sorting
sortData = struct('ncell', ncell, ...
                  'npart', npart, ...
                  'cell_n', zeros(ncell, 1), ...
                  'index', zeros(ncell, 1), ...
                  'Xref', zeros(npart, 1));
              
% Loop for desired number of timesteps
colSum = 0;
strikeSum = [0 0];
t_max = input('Enter total period to simulate (ms): ');
nstep = round(t_max * 1e-3 / tau);
fprintf('Number of steps = %.1f\n', nstep);

dl_bin = 5e-3;
xbins = dl_bin/2:dl_bin:L;
ybins = dl_bin/2:dl_bin:h;

string = 'Initial distribution';

fig = figure;
values = hist3([y x], {ybins xbins});
imagesc(xbins, ybins, values)
colorbar
axis equal
caxis_v = caxis;
title(string);
xlabel('x (m)');
ylabel('y (m)');
png = sprintf('results/%06dmks.png', 0);
saveas(fig, png);
pause(1);

close(fig);
fig = figure();
hist(x, lbin);
title(string);
xlabel('x (m)');
ylabel('Density');
png = sprintf('results/h%06dmks.png', 0);
saveas(fig, png);
pause(1);

tic;
for istep = 1:nstep;
    
    % Move all particles
    [x, y, v, strikes] = mover(x, y, v, npart, L, h, mpv, tau);
    strikeSum = strikeSum + strikes;
    
    % Sort the particles into cells
    sortData = sorter(x, y, L, dl, sortData);
    
    % Evaluate collisions among the particles
    [v, vrmax, selxtra, col] = ...
        collider(v, vrmax, tau, selxtra, coeff, sortData);
    colSum = colSum + col;
    
    % periodically display progress
    if( rem(istep, 100) < 1)
        elapsed_time = toc;
        
        fprintf('Finished %g of %g steps, Collisions = %g\n', ...
            istep, nstep, colSum);
        
        string = sprintf(...
            'Done %g of %g steps; %g collisions; t = %gs, %.2f hours left', ...
            istep, nstep, colSUm, istep*tau, elapsed_time*(nstep/istep - 1)/3600.0);
            
        close(fig);
        fig = figure();
        values = hist3([y x], {ybins, xbins});
        imagesc(xbins, ybins, values)
        colorbar
        axis equal
        caxis(caxis_v)
        title(string);
        xlabel('x (m)');
        ylabel('y (m)');
        png = sprintf('results/%06dmks.png', round(istep * tau * 1e6));
        saveas(fig, png);
        pause(0.1);
        
        close(fig);
        fig = figure();
        hist(x, lbin);
        title(string);
        xlabel('x (m)');
        ylabel('Density');
        png = sprintf('results/h%06dmks.png', round(istep * tau * 1e6));
        saveas(fig, png);
        pause(0.1);
        
    end
end

fprintf('Press Enter to exit\n');
pause;

        
     