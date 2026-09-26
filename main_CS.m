%%
% Solves Linearized Covariance Steering using 
% F. Liu, G. Rapakoulias, and P. Tsiotras, “Optimal covariance steering
% for discrete-time linear stochastic systems,” IEEE Transactions on
% Automatic Control, vol. 70, no. 4, 2025.
%%

clear;
close all;

initialization;
params_simulation;

%% Generate Reference
Ak = zeros(nx,nx,N_edges);
Bk = zeros(nx,nu,N_edges);
ck = zeros(nx,N_edges);
Rk = eye(nu);

refTraj.tk = tk;
STM0 = reshape(eye(nx),nx^2,1);
c0 = zeros(nx,1);

refTraj.muk_ref = zeros(nx,N_nodes);
refTraj.muk_ref(:,1) = mu0;
for k = 1:N_edges
    sol = ode45(@(t,X) nonlinearEOM(t,X,params_ODE,1), [refTraj.tk(k),refTraj.tk(k+1)], [refTraj.muk_ref(:,k);STM0;c0], options.ode);
    zend = deval(sol, refTraj.tk(k+1));

    Ak(:,:,k) = reshape(zend(nx+1:nx+nx^2), nx, nx);
    Bk(:,:,k) = Ak(:,:,k) * B;
    ck(:,k) = Ak(:,:,k) * zend(nx+nx^2+1:nx+nx^2+nx);

    refTraj.muk_ref(:,k+1) = zend(1:nx);
end


%% Solve Convex Problem
Pf = (constraint.r_3sigma/3)^2*eye(nx/2);

cvx_begin
    variable muk(nx,N_nodes)
    variable Pk(nx,nx,N_nodes) symmetric

    % Convexification Variables
    variable Uk(nu,nx,N_edges)
    variable Yk(nu,nu,N_edges) symmetric

    variable ub_k(nu,N_edges)

    subject to

        % Mean Constraint
        muk(:,1) == mu0;
        for k = 1:N_edges
            muk(:,k+1) == Ak(:,:,k)*muk(:,k) + Bk(:,:,k)*ub_k(:,k) + ck(:,k);
        end

        % Covariance
        Pk(:,:,1) == P0;
        for k = 1:N_edges
            Pk(:,:,k+1) == Ak(:,:,k)*Pk(:,:,k)*Ak(:,:,k).' + Bk(:,:,k)*Uk(:,:,k)*Ak(:,:,k).' + Ak(:,:,k)*Uk(:,:,k).'*Bk(:,:,k).' + Bk(:,:,k)*Yk(:,:,k)*Bk(:,:,k).';

            [Pk(:,:,k),Uk(:,:,k).';
             Uk(:,:,k),Yk(:,:,k)] == semidefinite(nx+nu);
        end

        % Final Constraints
        Pf - H_r*Pk(:,:,end)*H_r.' == semidefinite(nx/2);
        muk(:,end) == constraint.muf;

        OBJ = 0;
        for k = 1:N_edges
            OBJ = OBJ + quad_form(ub_k(:,k),Rk) ...
                      + trace(Rk*Yk(:,:,k));
        end

    minimize OBJ
cvx_end

muk    = full(muk);
Pk = full(Pk);
Uk    = full(Uk);
Yk    = full(Yk);
ub_k  = full(ub_k);

K_k = zeros(nu,nx,N_edges);
for k = 1:N_edges
    K_k(:,:,k) = Uk(:,:,k)*Pk(:,:,k)^(-1);
end

refTraj.Pk_ref = Pk;
refTraj.ub_k = ub_k; 
refTraj.K_k = K_k;

% Get Continuous Mean
muk_hist = zeros(nx,N_nodes);
muk_hist(:,1) = mu0;
for k = 1:N_edges
    mu_plus = muk_hist(:,k) + B*refTraj.ub_k(:,k);

    if k == 1
        sol = ode45(@(t,X) nonlinearEOM(t,X,params_ODE), [refTraj.tk(k),refTraj.tk(k+1)], mu_plus, options.ode);
    else
        sol = odextend(sol, @(t,X) nonlinearEOM(t,X,params_ODE), [refTraj.tk(k),refTraj.tk(k+1)], mu_plus, options.ode);
    end
    muk_hist(:,k+1) = sol.y(1:nx,end);
end

t_hist = sol.x;
mu_hist = deval(sol,t_hist);
refTraj.t_hist = t_hist;
refTraj.mu_hist = mu_hist;

% Dummy Variables so no plotting issues later
refTraj.ns = 1;
refTraj.xk_ref = nan(nx*refTraj.ns,N_nodes);

save("refTraj_CS.mat","refTraj")

%% Plotting
figure;
plot(refTraj.mu_hist(1,:),refTraj.mu_hist(2,:),'b-')
grid on
axis equal
