function [Agamma_k,Akurt_k] = linearize_moments(refTraj)
    params_simulation;
    wi = refTraj.wi;
    ns = refTraj.ns;

    Agamma_k = zeros(nx, nx*ns, N_nodes);
    Akurt_k = zeros(nx, nx*ns, N_nodes);
    for k = 1:N_nodes  
        zik = reshape(refTraj.zk_ref(:,k),nx,ns);

        % Skewness
        Agamma_k(:,:,k) = linearize_Cm(wi, zik, 3);

        % Kurtosis
        Akurt_k(:,:,k) = linearize_Cm(wi, zik, 4);
    end
end