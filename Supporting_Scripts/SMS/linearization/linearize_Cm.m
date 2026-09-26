function Am = linearize_Cm(wi, zi_ref, m)
    wi = vpa(wi);
    zi_ref = vpa(zi_ref);

    wi = vpa(wi);
    zi_ref = vpa(zi_ref);

    %% Reference Sigma Points
    [nx,ns] = size(zi_ref);

    %% m-th standardized moment
    Ez2 = sum(wi.*(zi_ref.^2),2);
    Ezm = sum(wi.*(zi_ref.^m),2);

    % More efficient way
    Ez_i = repmat(eye(nx),1,ns);
    Ezm_minius_i = repmat(eye(nx),1,ns);
    Ez_i(Ez_i == 1) = wi.*zi_ref;
    Ezm_minius_i(Ezm_minius_i == 1) = wi.*(zi_ref.^(m-1));

    Am = m*Ezm_minius_i.*Ez2.^(-m/2) - m * Ezm.* (Ez2.^(-m/2 - 1)).* Ez_i;
    Am = double(Am);
end