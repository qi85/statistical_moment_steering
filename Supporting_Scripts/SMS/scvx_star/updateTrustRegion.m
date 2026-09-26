function opt = updateTrustRegion(opt, opt_params)
    %% Trust Region Update
    if ((1 - opt_params.eta2) <= opt.rho) && (opt.rho <= (1 + opt_params.eta2))
        opt.r_trust = min(opt_params.r_max, opt.r_trust*opt_params.alpha2);

        opt.status.trust_update = "(^)";
    elseif ((1 - opt_params.eta1) <= opt.rho) && (opt.rho <= (1 + opt_params.eta1))
        opt.r_trust = opt.r_trust;

        opt.status.trust_update = "(-)";
    else
        opt.r_trust = max(opt_params.r_min, opt.r_trust/opt_params.alpha1);

        opt.status.trust_update = "(v)";
    end
end