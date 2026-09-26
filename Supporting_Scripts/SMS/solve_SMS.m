function [refTraj, opt_hist] = solve_SMS(refTraj, opt, opt_params, flag_plot)
    %%
    % Solve Statistical Moment Steering using SCvx*
    % K. Oguri, “Successive convexification with feasibility guarantee via
    % augmented lagrangian for non-convex optimal control problems,”
    % IEEE Conference on Decision and Control, 2023.
    %%

    params_simulation;
    ns = refTraj.ns;
    
    % Params History
    opt_hist.opt_params = opt_params;
    opt_hist.iter = [];
    opt_hist.deltaJ = [];
    opt_hist.infeas = [];
    opt_hist.iter_accept = [];
    
    %% Initialize SCVX Optimization Variables
    opt.status.iter = 0;
    opt.status.iter_accept = "Y";
    
    opt.deltaJ = inf;
    opt.deltaL = inf;
    opt.infeas = inf;
    opt.delta = inf;
    opt.r_trust = opt_params.r_init;
    opt.w = opt_params.w_init;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Lagrange Multipliers
        opt.lambda_dyn = zeros(nx*ns*N_edges,1); 
        opt.mu_skewfinal = zeros(2*size(H_r,1),1); % position skewness, both +/- sides
        opt.mu_kurtfinal = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Invariant Matrices
    [Az, Amu] = linearize_z(refTraj.wi, nx, ns); 
    cvx_Jacobians.Az = Az;
    cvx_Jacobians.Amu = Amu;

    %% Run SCVX*
    tic;
    while abs(opt.deltaJ) > opt_params.tol_opt || opt.infeas > opt_params.tol_feas || opt.status.iter_accept == "N"
        opt.status.iter = opt.status.iter + 1;
    
        %% Linearization
        if opt.status.iter_accept == "Y"
            % Update moments with new sigma points
            refTraj = updateMoments(refTraj);

            % Linearize Moments
            [Agamma_k,Akurt_k] = linearize_moments(refTraj); 
            cvx_Jacobians.Agamma_k = Agamma_k;
            cvx_Jacobians.Akurt_k = Akurt_k;

            % Linearize Sigma Point Dynamics
            [Ak,Bk,ck] = linearize_dyn(refTraj); 
            cvx_Jacobians.Ak = Ak;
            cvx_Jacobians.Bk = Bk;
            cvx_Jacobians.ck = ck;
        end
        
        %% Convex Subproblem
        [cvxTraj, L_subproblem, L_struct, cvx_status] = cvx_subproblem(opt, refTraj, cvx_Jacobians);
        if cvx_status ~= "Solved"
            fprintf("\n CVX Warning: "+cvx_status+"\n")
            break
        end

        %% Nonlinear Problem and Infesibility
        % Ref solution
        [J_ref,~,~,J_struct_ref,g_ref,h_ref] = cvx_nonlinear(opt, refTraj);
    
        % New solution
        [J_star, infeas, xk_nl, J_struct_new, g_all, h_all] = cvx_nonlinear(opt, cvxTraj);
    
        opt.deltaJ = J_ref - J_star;
        opt.deltaL = J_ref - L_subproblem;
        opt.infeas = infeas;
    
        % if opt.deltaL < 0
        % 	dbstop; % something is wrong - stop for debug
        % end

        if flag_plot
            %% CVX vs Ref vs Nonlinear
            fig = figure(1000);
            clf(fig)
            for k = 1:N_nodes
                xik_ref = reshape(refTraj.xk_ref(:,k), nx ,ns);
                plot(xik_ref(1,:),xik_ref(2,:),'kx')
                hold on
    
                xik_cvx = reshape(cvxTraj.xk_ref(:,k), nx ,ns);
                plot(xik_cvx(1,:),xik_cvx(2,:),'bx')
                hold on
    
                xik_nl = reshape(xk_nl(:,k), nx ,ns);
                plot(xik_nl(1,:),xik_nl(2,:),'rx')
                hold on
            end
            axis equal
            grid on
            xlabel('x')
            ylabel('y')

            legend('Reference','CVX','Nonlinear')

            %% Dynamics Violation
            eta_dyn_nl = reshape(g_all.g_dyn_nl,nx*ns,N_edges);
            eta_norm = zeros(ns, N_edges);
            for k = 1:N_edges
                eta_k = reshape(eta_dyn_nl(:,k), nx, ns);
                eta_norm(:,k) = vecnorm(eta_k, 2, 1);
            end

            fig2 = figure(1001);
            clf(fig2)
           
            colormap(fig2, cool)
            cmap = cool(256);
            cmin = min(eta_norm,[],'all');
            cmax = max(eta_norm,[],'all');
            
            xik_nl = reshape(xk_nl(:,1), nx ,ns); 
            plot(xik_nl(1,:),xik_nl(2,:),'k.')
            hold on
            for k = 1:N_edges
                xik_nl = reshape(xk_nl(:,k+1), nx ,ns); % Plot the next one
  
                for i = 1:ns
                    idx = round(1 + 255*(eta_norm(i,k) - cmin)/(cmax - cmin + eps));
                    plot(xik_nl(1,i), xik_nl(2,i), '.','Color',cmap(idx,:))
                end

                hold on
            end
            axis equal
            grid on
            xlabel('x')
            ylabel('y')
            colorbar
            clim([cmin cmax])
            ylabel(colorbar, 'Dynamical Violation')
        end
    
        %% Step Acceptence Criteria
        if abs(opt.deltaL) == 0
            opt.rho = 1;
        else
            opt.rho = opt.deltaJ/opt.deltaL;
        end

        %% Step Acceptence
        if ((1 - opt_params.eta0) <= opt.rho) && (opt.rho <= (1 + opt_params.eta0)) 
            opt.status.iter_accept = "Y";
            refTraj = cvxTraj; % Accept Solution 
    
            % Lagrange Multiplier Update
            opt = updateLM(opt, opt_params, g_all, h_all);
        else
            opt.status.iter_accept = "N";
        end
    
        %% Trust Region Update
        opt = updateTrustRegion(opt, opt_params);

        %% Save History
        opt_hist.iter = [opt_hist.iter;opt.status.iter];
        opt_hist.deltaJ = [opt_hist.deltaJ;opt.deltaJ];
        opt_hist.infeas = [opt_hist.infeas;opt.infeas];
        opt_hist.iter_accept = [opt_hist.iter_accept;opt.status.iter_accept];
    
        %% Output Message
        if mod(opt.status.iter,20) == 1
            fprintf("\n\n")
            fprintf("Iter   Accept?         rho        deltaL        deltaJ         Infeas         Trust        w    \n")
            fprintf("================================================================================================\n")
        end
        status_output = sprintf("%4d%8s  %12.4e%14.4e%14.4e%15.4e%14.3e%3s%10.2e\n",opt.status.iter, ...
                                                                       opt.status.iter_accept, ...
                                                                       opt.rho, ...
                                                                       opt.deltaL, ...
                                                                       opt.deltaJ, ...
                                                                       opt.infeas, ...
                                                                       opt.r_trust, opt.status.trust_update, ...
                                                                       opt.w);
        fprintf(status_output)
    end
    toc;
end
