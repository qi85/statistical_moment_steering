function f0 = objective_fcn(objective, refTraj, objective_mode, cvx_vars)
    params_simulation;
    ns = refTraj.ns;
    Rk = eye(nu);

    if nargin < 4
        cvx_vars = struct();
    end

    switch objective
        case "cvar"
            f0 = objective_cvar(refTraj, objective_mode, cvx_vars, N_edges, ns);

        case "dv99"
            f0 = objective_dv99(refTraj, objective_mode, cvx_vars, N_edges, nu, n_chi2);

        case "quadratic"
            f0 = objective_quadratic(refTraj, objective_mode, cvx_vars, N_edges, Rk);

        otherwise
            error("Unknown objective function: %s", objective);
    end
end

function f0 = objective_cvar(refTraj, objective_mode, cvx_vars, N_edges, ns)
    alpha = 0.99;
    ui_norm_total = zeros(1,ns);
    for k = 1:N_edges
        if objective_mode == "cvx"
            uik = cvx_vars.uik_hist(:,:,k);
            ui_norm_total = ui_norm_total + norms(uik,2,1);
        else
            zik = reshape(refTraj.zk_ref(:,k),[],ns);
            uik = refTraj.ub_k(:,k) + refTraj.K_k(:,:,k)*zik;
            ui_norm_total = ui_norm_total + vecnorm(uik,2,1);
        end
    end

    if objective_mode == "cvx"
        CVaR_tau = cvx_vars.CVaR_tau;
        ui_norm_zeros = [ui_norm_total-CVaR_tau;zeros(1,ns)];
        f0 = CVaR_tau + 1/(1-alpha)*sum(refTraj.wi.*max(ui_norm_zeros,[],1));
    else
        f0 = nonlinearCVaR(ui_norm_total, refTraj.wi, alpha);
    end
end

function f0 = nonlinearCVaR(losses, wi, alpha)
    losses = losses(:);
    wi = wi(:);
    ns = length(wi);

    % CVar = inf_tau (...)
    cvx_begin quiet
        variable CVaR_tau(1)
        variable tail_excess(ns)

        minimize(CVaR_tau + 1/(1-alpha)*(wi.'*tail_excess))

        subject to
            tail_excess >= losses - CVaR_tau;
            tail_excess >= 0;
    cvx_end

    if ~contains(cvx_status, "Solved")
        error("Nonlinear CVaR tau optimization failed: %s", cvx_status);
    end

    CVaR_tau = full(CVaR_tau);
    f0 = CVaR_tau + 1/(1-alpha)*sum(wi.*max(losses-CVaR_tau,0));
end

function f0 = objective_dv99(refTraj, objective_mode, cvx_vars, N_edges, nu, n_chi2)
    f0 = 0;
    for k = 1:N_edges
        if objective_mode == "cvx"
            ub_k = refTraj.ub_k(:,k) + cvx_vars.dub_k(:,k);
            K_k = refTraj.K_k(:,:,k) + cvx_vars.dK_k(:,:,k);
        else
            ub_k = refTraj.ub_k(:,k);
            K_k = refTraj.K_k(:,:,k);
        end

        f0 = f0 + norm(ub_k,2) ...
             + n_chi2(0.01, nu) * norm(K_k*chol(refTraj.Pk_ref(:,:,k),'lower'),2);
    end
end

function f0 = objective_quadratic(refTraj, objective_mode, cvx_vars, N_edges, Rk)
    f0 = 0;
    for k = 1:N_edges
        if objective_mode == "cvx"
            ub_k = refTraj.ub_k(:,k) + cvx_vars.dub_k(:,k);
            K_k = refTraj.K_k(:,:,k) + cvx_vars.dK_k(:,:,k);
        else
            ub_k = refTraj.ub_k(:,k);
            K_k = refTraj.K_k(:,:,k);
        end

        control_cov_norm = norm(chol(Rk)*K_k*chol(refTraj.Pk_ref(:,:,k),'lower'),'fro');
        if objective_mode == "cvx"
            ub_cost = quad_form(ub_k,Rk);
            control_cov_cost = pow_pos(control_cov_norm,2);
        else
            ub_cost = ub_k.'*Rk*ub_k;
            control_cov_cost = control_cov_norm^2;
        end

        f0 = f0 + ub_cost + control_cov_cost;
    end
end
