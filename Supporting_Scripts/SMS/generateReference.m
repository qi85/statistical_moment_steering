function refTraj = generateReference()
    %%
    % Generates initial reference solution for statistical moment steering
    % Uses "Scaled Linear Covariance" technique from 
    % D. C. Qi, K. Oguri, P. Singla, and M. R. Akella, “Non-gaussian
    % distribution steering in nonlinear dynamics with conjugate unscented
    % transformation,” arXiv preprint, 2025.
    %%

    params_simulation;

    refTraj.tk = tk;

    %% Uncontrolled Lin Cov Scaled
    refTraj.ub_k = zeros(nu, N_edges);
    refTraj.K_k = zeros(nu, nx, N_edges);

    % Nominal Trajectory
    sol_nominal = ode45(@(t,X) nonlinearEOM(t,X,params_ODE), tk, mu0, options.ode);
    Xk_nominal = deval(sol_nominal, tk);

    Pk_linear = zeros(nx,nx,N_nodes);
    Pk_linear(:,:,1) = P0;

    STM0 = reshape(eye(nx),nx^2,1);
    c0 = zeros(nx,1);

    % Scaled Linear Covariance
    for k = 1:N_edges
        [~,X_hist] = ode45(@(t,X) nonlinearEOM(t,X,params_ODE,1), [refTraj.tk(k),refTraj.tk(k+1)], [Xk_nominal(:,k);STM0;c0], options.ode);
        Ak = reshape(X_hist(end,nx+1:nx^2 + nx),nx,nx);

        Pk_linear_next = Ak*Pk_linear(:,:,k)*Ak.';

        %[~,D0] = eig(Pk_linear(:,:,k));
        [V,~] = eig(Pk_linear_next);

        Pk_linear_normalized = V*P0*V.';

        Pk_linear(:,:,k+1) = Pk_linear_normalized;
    end

    %% Initial Reference Sigma Points
    [xi0,wi,ns] = sigmaPoints_CUT4G(mu0, P0);
    refTraj.wi = wi;
    refTraj.ns = ns;
    refTraj.xi0 = xi0; % Initial Distribution

    xk_ref = zeros(nx*ns,1);
    for k = 1:N_nodes
        [xik,~,~] = sigmaPoints_CUT4G(Xk_nominal(:,k), Pk_linear(:,:,k));

        xk_ref(:,k) = reshape(xik,nx*ns,1);
    end
    refTraj.xk_ref = xk_ref;

    if 0
        %% Debug Plot
        X_nominal = deval(sol_nominal, sol_nominal.x);

        figure;
        plot(X_nominal(1,:),X_nominal(2,:),'b-')
        hold on
        grid on
        axis equal
        for k = 1:N_nodes
            xik = reshape(refTraj.xk_ref(:,k),nx,ns);
            plot(xik(1,:),xik(2,:),'kx','LineStyle','none')
        end
    end
end