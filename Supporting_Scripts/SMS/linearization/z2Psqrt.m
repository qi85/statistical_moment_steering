function P_sqrt = z2Psqrt(wi,zik)
    % Have to do this for cvx 
    [nx,ns] = size(zik);
    w_matrix = sqrt(wi).*ones(nx,ns);

    P_sqrt = w_matrix.*zik;
end