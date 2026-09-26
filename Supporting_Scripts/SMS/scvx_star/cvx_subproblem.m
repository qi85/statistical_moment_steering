function [cvxTraj, L_subproblem, L_struct, cvx_status] = cvx_subproblem(opt, refTraj, cvx_Jacobians)
    params_simulation;
    ns = refTraj.ns;

    Az = cvx_Jacobians.Az;
    Amu = cvx_Jacobians.Amu;
    Agamma_k = cvx_Jacobians.Agamma_k;
    Akurt_k = cvx_Jacobians.Akurt_k;

    Ak = cvx_Jacobians.Ak;
    Bk = cvx_Jacobians.Bk;
    ck = cvx_Jacobians.ck;

    n_skew = 2*size(H_r,1); % Position, by positive and negative skewness constraint
    
    cvx_begin quiet
        %% Optimization Variables
        variable dx_k(nx*ns, N_nodes) 
        variable dub_k(nu, N_edges) 
        variable dK_k(nu, nx, N_edges) 
    
        % Moments
        variable dz_k(nx*ns, N_nodes)
        variable muk(nx, N_nodes)
        
        % Slack Variables
        variable eta_dyn(nx*ns,N_edges)
        variable eta_skewfinal(n_skew)
        variable eta_kurtfinal(1)

        if opt.objective == "cvar"
            variable CVaR_tau(1)
            variable uik_hist(nu,ns,N_edges)
        end
    
        subject to
            %% Dynamical Constraint
            for k = 1:N_edges
                for i = 1:ns
                    s_ind = Es(i);
    
                    xk = refTraj.xk_ref(s_ind,k) + dx_k(s_ind,k);
                    xkp1 = refTraj.xk_ref(s_ind,k+1) + dx_k(s_ind,k+1);
                    
                    uk = (refTraj.ub_k(:,k) + dub_k(:,k)) ...
                           + refTraj.K_k(:,:,k)*refTraj.zk_ref(s_ind,k) ...
                           + dK_k(:,:,k)*refTraj.zk_ref(s_ind,k) + refTraj.K_k(:,:,k)*dz_k(s_ind,k);
    
                    eta_dyn(s_ind,k) == xkp1 - (Ak{i,k}*xk + Bk{i,k}*uk + ck{i,k});

                    if opt.objective == "cvar"
                        uik_hist(:,i,k) == uk;
                    end
                end
            end
    
		    %% Trust Region
		    norm(dx_k(:),inf)  <= opt.r_trust;
            norm(dub_k(:),inf) <= opt.r_trust;
            norm(dK_k(:),inf)  <= opt.r_trust;
	    
            %% Initial Distrubtion Contraint
            refTraj.xk_ref(:,1) + dx_k(:,1) == reshape(refTraj.xi0,nx*ns,1);
    
            %% Convex Moments
            for k = 1:N_nodes
                dz_k(:,k) == Az*dx_k(:,k); 
                muk(:,k) == refTraj.muk_ref(:,k) + Amu*dx_k(:,k);
            end
    
            %% Additional Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Mean
            constraint.muf == muk(:,end);
    
            % Covariance
            for k = N_nodes
                norm(H_r*z2Psqrt(refTraj.wi,  reshape(refTraj.zk_ref(:,k) + dz_k(:,k),nx,ns)),   2) <= constraint.r_3sigma/3;
            end
    
            % Final Skewness
            if opt.constraints.skew
                [ H_r*(refTraj.gammak_ref(:,end) + Agamma_k(:,:,end)*dz_k(:,end)) - opt.constraints.skew_eps;
                 -H_r*(refTraj.gammak_ref(:,end) + Agamma_k(:,:,end)*dz_k(:,end)) - opt.constraints.skew_eps] <= eta_skewfinal;
                eta_skewfinal >= 0;
            end
    
            %% Objective Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cvx_vars.dub_k = dub_k;
            cvx_vars.dK_k = dK_k;
            cvx_vars.dz_k = dz_k;
            if opt.objective == "cvar"
                cvx_vars.CVaR_tau = CVaR_tau;
                cvx_vars.uik_hist = uik_hist;
            end
            f0_cvx = objective_fcn(opt.objective, refTraj, "cvx", cvx_vars);

            %% Penalty Functions
            % Dynamics
            g_dyn_cvx = reshape(eta_dyn,nx*ns*N_edges,1);
            penalty_dyn_cvx = penalty_fcn(opt.w, opt.lambda_dyn, g_dyn_cvx,[],[]);
    
            % Skew Final
            if opt.constraints.skew 
                h_skewfinal_cvx = eta_skewfinal;
                penalty_skewfinal_cvx = penalty_fcn(opt.w, [], [], opt.mu_skewfinal, h_skewfinal_cvx);
            else
                penalty_skewfinal_cvx = 0;
            end
    
            L_subproblem = f0_cvx + penalty_dyn_cvx + penalty_skewfinal_cvx;

            minimize L_subproblem
    cvx_end
    % Change Data Type Just in Case
    dx_k = full(dx_k);
    dz_k = full(dz_k);
    dub_k = full(dub_k);
    dK_k = full(dK_k);
    
    % Full State of CVX Solution
    cvxTraj = refTraj;
    cvxTraj.xk_ref = refTraj.xk_ref + dx_k;
    cvxTraj.zk_ref = refTraj.zk_ref + dz_k;
    cvxTraj.ub_k = refTraj.ub_k + dub_k;
    cvxTraj.K_k = refTraj.K_k + dK_k;

    if opt.objective == "cvar"
        cvxTraj.CVaR_tau = CVaR_tau;
    end
    
    L_struct.f0_cvx = f0_cvx;
    L_struct.penalty_dyn_cvx = penalty_dyn_cvx;
    L_struct.penalty_skewfinal_cvx = penalty_skewfinal_cvx; 
    
    if 0
        %% Debug Plot
        figure;
        for k = 1:N_nodes
            xik = reshape(refTraj.xk_ref(:,k),nx,ns);
            plot(xik(1,:),xik(2,:),'kx','LineStyle','none')
            hold on
    
            xik = reshape(cvxTraj.xk_ref(:,k),nx,ns);
            plot(xik(1,:),xik(2,:),'rx','LineStyle','none')
        end
        grid on
        axis equal
    end

end
