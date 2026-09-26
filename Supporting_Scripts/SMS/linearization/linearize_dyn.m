function [Ak_cell,Bk_cell,ck_cell] = linearize_dyn(refTraj)
    params_simulation;
    ns = refTraj.ns;

    STM0 = reshape(eye(nx),nx^2,1);
    c0 = zeros(nx,1);

    Ak_cell = cell(ns,N_edges);
    Bk_cell = cell(ns,N_edges);
    ck_cell = cell(ns,N_edges);
    for k = 1:N_edges
        xik = reshape(refTraj.xk_ref(:,k),nx,ns);
        zik = reshape(refTraj.zk_ref(:,k),nx,ns);
        for i = 1:ns
            xik_plus = xik(:,i) + B*(refTraj.ub_k(:,k) + refTraj.K_k(:,:,k)*zik(:,i));
            sol = ode45(@(t,X) nonlinearEOM(t,X,params_ODE,1), [refTraj.tk(k),refTraj.tk(k+1)], [xik_plus;STM0;c0], options.ode);
            zend = deval(sol, refTraj.tk(k+1));

            Ak = reshape(zend(nx+1:nx+nx^2), nx, nx);
            Bk = Ak * B;
            ck = Ak * zend(nx+nx^2+1:nx+nx^2+nx);

            Ak_cell{i,k} = Ak; 
            Bk_cell{i,k} = Bk;
            ck_cell{i,k} = ck;
        end
    end
end