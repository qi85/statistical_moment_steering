%%
% Generate set of parameters called by scripts and functions
%%

options.ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);

%% Simulation
load("curve.mat")
% Containts params_ODE

N_nodes = 15;
Tf = curve.Tf;
N_edges = N_nodes - 1;
tk = linspace(0,Tf,N_nodes);

%% Dynamics & Control
nx = 4;
nu = 2;
B = [zeros(2);eye(2)];

%% Initial Conditions
mu0 = curve.X0;
P0 = diag([0.05,0.05,0.05,0.05].^2);

%% Constraints
constraint.muf = curve.Xf;
constraint.r_3sigma = 0.25;

%% General Helpers
H_r = [eye(2),zeros(2)];
H_v = [zeros(2),eye(2)];
n_chi2 = @(epsilon, n) sqrt(chi2inv(1 - epsilon, n));
Es = @(i) (nx*(i-1)+1):(nx*(i)); % indices of ith sigma point in aggregate sigma point vector 