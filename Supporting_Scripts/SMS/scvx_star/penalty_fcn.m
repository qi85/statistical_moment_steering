function [penalty] = penalty_fcn(w,lambda,g,mu,h)
    % h must be non-negative
    if isnumeric(h)
	    h(h<0) = 0;
    end
    
    % equality
    if isempty(lambda)
	    penalty_eq = 0;
    else
        penalty_eq = dot(lambda, g) + w/2 * dot(g,g);
    end
    
    % inequality
    if isempty(mu)
	    penalty_ineq = 0;
    else
        penalty_ineq = dot(mu, h) + w/2 * dot(h,h);
    end
    
    % overall
    penalty = penalty_eq + penalty_ineq;
end



