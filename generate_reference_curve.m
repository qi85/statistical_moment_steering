%%
% Generates a reference trajectory for the nonlinear oscillator.
%%

clear;
close all;

initialization;

options.fsolve = optimset('TolFun',1e-8,'TolX',1e-8);
options.ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);

params_ODE.alpha = 0.25; % spring constant
params_ODE.beta = 0.75; % nonlinear
params_ODE.kappa = 0.5; % coupling

x0 = 2;
ydot0_guess = 5;
Tf_guess = 12;

X0_optm = fsolve(@(X_optm) nonlinearEOM_reflection(X_optm, x0, params_ODE), [ydot0_guess;Tf_guess], options.fsolve);

%% Plot Trajectory
X0 = [x0;0;0;X0_optm(1)];
Tf = X0_optm(end);
[~,X_sol] = ode45(@(t,X) nonlinearEOM(t,X,params_ODE), [0,Tf], X0, options.ode);
X_sol = X_sol.';
Xf = X_sol(:,end);

figure;
plot(X_sol(1,:),X_sol(2,:),'b-')
hold on
plot(X0(1),X0(2),'g*')
hold on
plot(Xf(1),Xf(2),'r*')
grid on
axis equal

%% Save Data
curve.X0 = X0;
curve.Xf = Xf;
curve.Tf = Tf;
save("curve.mat",'curve','params_ODE')

%% Shooting Function
function Psi = nonlinearEOM_reflection(X_optm, x0, params_ODE)
    options.ode = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12);

    ydot0 = X_optm(1);
    Tf = X_optm(end);

    X0 = [x0;0;0;ydot0];

    [~,X_hist] = ode45(@(t,X) nonlinearEOM(t,X,params_ODE), [0,Tf], X0, options.ode);
    Xf = X_hist(end,:).';

    %% F = 0 and DF
    Psi = Xf(2:3);
end