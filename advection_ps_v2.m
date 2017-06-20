%% Advection Equation - Pseudospectral
% JAM 4/13/17.
%
% Solves the advection equation
%
% $$ \partial_t f + v\partial_x f = 0, $$
%
% with constant $v$ over a periodic interval $[-L/2,+L/2]$ using a pseudospectral
% method with RK4 timestepping. 
%

%% Parameters
% Specify the simulation parameters
clearvars; % clear all variables
L = 2*pi;  % Length of x domain [-L/2,+L/2]
Nx = 128;  % Number of x computational gridpoints (of form 2^n)
dx = L/Nx; % grid spacing
v = 1.0;   % velocity
T = 2*2*pi;  % final time
cfl = 0.5; % Courant number v dt/dx
dt = cfl*dx/v;  % timestep
Nt = round(T/dt); % approx number of timesteps
                  % simulation final time is Nt*dt
                  
%% Preliminary Calculations
% Set up the grid and mode numbers
dto2 = dt/2.0;
xf = linspace(-L/2,+L/2,Nx+1);  % full discretized x domain
x = xf(2:Nx+1); % computational x grid (periodic)
%k = (2*pi/L)*[0:(Nx/2-1) Nx/2 (-Nx/2+1):-1]; % wavenumbers
k = (2*pi/L)*[0:(Nx/2-1) 0 (-Nx/2+1):-1]; % wavenumbers (Trefethen)

%% Initial Distribution
% Must be periodic over the interval
%k0=10.0; f0 = cos(k0*x);
%f0 = exp(-x.*x/0.5);
f0 = exp(-x).*sin(x)
%f0 = (abs(x)<1);
%
% Some plotting setup
plotrange = [-L/2,L/2,min(f0)-0.1,max(f0)+0.1];

%% Time stepping
% Step forward in time and show the evolution of the distribution
f = f0;
for n = 1:Nt
    % 2nd order R-K
    % The argument of ifft should be symmetric to within roundoff;
    %   force it to be and return a real result
    % vdfdx1 = v*ifft(i*k.*fft(f),'symmetric');   
    % f2 = f - deltato2*vdfdx1;
    % vdfdx2 = v*ifft(i*k.*fft(f2),'symmetric');
    % f = f - deltat*vdfdx2;

    % 4th order R-K
    vdfdx1 = v*ifft(1i*k.*fft(f),'symmetric');
    f2 = f - dto2*vdfdx1;
    vdfdx2 = v*ifft(1i*k.*fft(f2),'symmetric');
    f3 = f - dto2*vdfdx2;
    vdfdx3 = v*ifft(1i*k.*fft(f3),'symmetric');
    f4 = f - dt*vdfdx3;
    vdfdx4 = v*ifft(1i*k.*fft(f4),'symmetric');
    f = f - dt*(vdfdx1/6 + vdfdx2/3 + vdfdx3/3 + vdfdx4/6);
    
    plot(x,f0,'-xb',x,f,'-r','LineWidth',1.0); 
    axis(plotrange); %pbaspect([2 1 1]) 
    pause(0.025)
end

%% Final Plot
% Add titles to final plot generated in the loop
ax = gca; 
ax.FontSize = 12;
title('Pseudospectral','interpreter','latex','FontSize',18);
xlabel('$x$','interpreter','latex','FontSize',18);
ylabel({'$f(x)$'},'interpreter','latex','FontSize',18);



