function [J_nl, infeas, xk_nl, J_struct, g_all, h_all] = cvx_nonlinear(opt, refTraj)
    params_simulation;
    ns = refTraj.ns;

    %% Nonlinear Cost
    J_nl = 0;

    %% Objective Function
    f0_nl = objective_fcn(opt.objective, refTraj, "nonlinear");

    J_nl = J_nl + f0_nl;
    J_struct.f0_nl = f0_nl;

    %% Dynamics
    eta_dyn = zeros(nx*ns,N_edges);
    xk_nl = zeros(nx*ns,N_nodes);
    xk_nl(:,1) = reshape(refTraj.xi0,nx*ns,1);
    for k = 1:N_edges
        for i = 1:ns
            s_ind = Es(i);

            xk = refTraj.xk_ref(s_ind,k);
            xkp1 = refTraj.xk_ref(s_ind,k+1);

            % Control
            uk = refTraj.ub_k(:,k) + refTraj.K_k(:,:,k)*refTraj.zk_ref(s_ind,k);
            xi_plus = xk + B*uk;

            % Nonlinear 
            [~,X_hist] = ode45(@(t,X) nonlinearEOM(t,X,params_ODE), [refTraj.tk(k),refTraj.tk(k+1)], xi_plus, options.ode);
            xk_nl(s_ind,k+1) = X_hist(end,:).';

            eta_dyn(s_ind,k) = xkp1 - xk_nl(s_ind,k+1);
        end
    end
    g_dyn_nl = reshape(eta_dyn,nx*ns*N_edges,1); 
    penalty_dyn_nl = penalty_fcn(opt.w, opt.lambda_dyn, g_dyn_nl,[],[]);

    J_nl = J_nl + penalty_dyn_nl;
    J_struct.penalty_dyn_nl = penalty_dyn_nl;
    g_all.g_dyn_nl = g_dyn_nl;

    infeas = norm(eta_dyn(:),inf);

    %% Other Constraints
    if opt.constraints.skew
        xik_ref = reshape(refTraj.xk_ref(:,end), nx ,ns);
        [~,~,gammaf,~] = momentProcessing(refTraj.wi, xik_ref);

        gamma_skewfinal_nl = H_r*gammaf;
        h_skewfinal_nl = [ gamma_skewfinal_nl - opt.constraints.skew_eps;
                          -gamma_skewfinal_nl - opt.constraints.skew_eps];
        penalty_skewfinal = penalty_fcn(opt.w,[], [], opt.mu_skewfinal, h_skewfinal_nl);
        J_nl = J_nl + penalty_skewfinal;
    
        J_struct.penalty_skewfinal = penalty_skewfinal;
        h_all.h_skewfinal_nl = h_skewfinal_nl;
    else
        h_all.h_skewfinal_nl = 0;
    end
end
