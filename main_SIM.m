function [X_sol, Xk_hist, uk_hist] = main_SIM(refTraj)
    %%
    % Runs single iteration of Monte Carlo from main_MC.m
    % Inputs: refTraj
    % Outputs: ODE solution (X_sol), Xk history (Xk_hist), control history (uk_hist)
    %%

    params_simulation;

    % Initial Gaussian
    Xk = mu0 + chol(P0,'lower')*randn(nx,1);

    % Start Sim
    Xk_hist = zeros(nx, N_nodes);
    Xk_hist(:,1) = Xk;
    
    uk_hist = zeros(nu, N_edges);
    for k = 1:N_edges
        % Control
        uk = refTraj.ub_k(:,k) + refTraj.K_k(:,:,k)*(Xk - refTraj.muk_ref(:,k));
        Xk_plus = Xk + B*uk;

        uk_hist(:,k) = uk;
        if k == 1
            sol = ode45(@(t,X) nonlinearEOM(t,X,params_ODE), [refTraj.tk(k),refTraj.tk(k+1)], Xk_plus, options.ode);
        else
            sol = odextend(sol,@(t,X) nonlinearEOM(t,X,params_ODE), [refTraj.tk(k),refTraj.tk(k+1)], Xk_plus, options.ode);
        end

        Xk = sol.y(1:nx,end);
        Xk_hist(:,k+1) = Xk;
    end

    X_sol = sol;
end