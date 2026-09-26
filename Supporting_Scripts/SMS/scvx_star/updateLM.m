function opt = updateLM(opt, opt_params, g_all, h_all)
        %% Multiplier and Stationary Tol Update
        if abs(opt.deltaJ) < opt.delta

            %% Lagrange Multiplier
            opt.lambda_dyn = opt.lambda_dyn + opt.w*g_all.g_dyn_nl;

            if opt.constraints.skew
                opt.mu_skewfinal = opt.mu_skewfinal + opt.w*h_all.h_skewfinal_nl;
                opt.mu_skewfinal(opt.mu_skewfinal<0) = 0;
            end
            % if opt.constraints.kurt == "ON"
            %     opt.mu_kurtfinal = opt.mu_kurtfinal + opt.w*h_all.h_kurtfinal_nl;
            %     opt.mu_kurtfinal(opt.mu_kurtfinal<0) = 0;
            % end

            %% w_{k+1}
            opt.w = opt_params.beta*opt.w;
            if opt.w > opt_params.w_max
                opt.w = opt_params.w_max;
            end

            %% delta_{k+1}
            if isinf(opt.delta)
                opt.delta = abs(opt.deltaJ);
            else
                opt.delta = opt.delta*opt_params.gamma;
            end
        end
end