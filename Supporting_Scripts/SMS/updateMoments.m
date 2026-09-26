function refTraj = updateMoments(refTraj)
    [nx_ns,N_nodes] = size(refTraj.xk_ref);

    ns = refTraj.ns;
    wi = refTraj.wi;
    nx = nx_ns/ns;

    zk_ref = zeros(nx*ns,N_nodes);
    muk_ref = zeros(nx,N_nodes);
    Pk_ref = zeros(nx,nx,N_nodes);
    gammak_ref = zeros(nx,N_nodes);
    kurtk_ref = zeros(nx,N_nodes);
    for k = 1:N_nodes
        xik_ref = reshape(refTraj.xk_ref(:,k), nx ,ns);

        % Reference Moments
        [muk,Pk,gammak,kurtk] = momentProcessing(wi, xik_ref);

        muk_ref(:,k) = muk;
        Pk_ref(:,:,k) = Pk;
        gammak_ref(:,k) = gammak;
        kurtk_ref(:,k) = kurtk;

        zik_ref = xik_ref - muk;
        zk_ref(:,k) = reshape(zik_ref,nx*ns,1);
    end

    refTraj.zk_ref = zk_ref;
    refTraj.muk_ref = muk_ref;
    refTraj.Pk_ref = Pk_ref;
    refTraj.gammak_ref = gammak_ref;
    refTraj.kurtk_ref = kurtk_ref;
end