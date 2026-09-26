function refTraj = propagateSigmaPoints_nonlinear(refTraj)
    params_simulation;
    ns = refTraj.ns;

    xik = zeros(nx,ns,N_nodes);
    xik(:,:,1) = refTraj.xi0;
    for i = 1:ns
        for k = 1:N_edges
            uk = refTraj.ub_k(:,k) + refTraj.K_k(:,:,k)*(xik(:,i,k) - refTraj.muk_ref(:,k));
            xi_plus = xik(:,i,k) + B*uk;

            if k == 1
                sol = ode45(@(t,X) nonlinearEOM(t,X,params_ODE), [refTraj.tk(k),refTraj.tk(k+1)], xi_plus, options.ode);
            else
                sol = odextend(sol, @(t,X) nonlinearEOM(t,X,params_ODE), [refTraj.tk(k),refTraj.tk(k+1)], xi_plus, options.ode);
            end
            xik(:,i,k+1) = sol.y(1:nx,end);
        end

        % Save Mean Trajectory
        if i == 1
            t_hist = sol.x;
            x_hist = zeros(nx*ns,length(t_hist));
        end
        x_hist((i-1)*nx+1:i*nx,:) = deval(sol,t_hist);
    end

    % Save Info
    for k = 1:N_nodes
        refTraj.xk_ref(:,k) = reshape(xik(:,:,k),nx*ns,1);
    end
    refTraj = updateMoments(refTraj);

    [~, Amu] = linearize_z(refTraj.wi, nx, ns); 
    mu_hist = Amu*x_hist;
    refTraj.t_hist = t_hist;
    refTraj.mu_hist = mu_hist;
end