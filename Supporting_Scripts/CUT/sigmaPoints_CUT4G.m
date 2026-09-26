function [xi,wi,N] = sigmaPoints_CUT4G(mu,P)
    %% 
    % Cubature rule from
    % N. Adurthi, P. Singla, and T. Singh, “Conjugate unscented transfor-
    % mation: Applications to estimation and control,” Journal of Dynamic
    % Systems, Measurement, and Control, vol. 140, no. 3, pp. 1–22, 2018
    %
    % Inputs: Mean (mu), Covariance Matrix (P)
    % Outputs: Cubature points (xi), weights (wi), and number of points (N)
    %%
   
    % n: Dimension
    n = length(mu);
    
    %% 4-th Order CUT - Gaussian
    % Distances
    r1 = sqrt((n+2)/2);
    r2 = sqrt((n+2)/(n-2));
    
    % Weights (w0 = 0)
    w1 = 4/(n+2)^2;
    w2 = (n-2)^2/(2^n*(n+2)^2);
    
    %% Calculate Sigma Points
    % Principle Axis
    sigma = zeros(n,2*n);
    sigma(:,1:n) = eye(n);
    sigma(:,n+1:end) = -eye(n);
    
    % Conjugate Axis    
    c_n = (dec2bin(0:2^n-1) - '0').';
    c_n(c_n == 0) = -1;
    
    %% Transform to Fit Given Mean and Covariance
    xi = [r1*sigma,r2*c_n]; % Sigma Points
    xi = chol(P,'lower')*xi + mu;

    wi = [repmat(w1,1,2*n),repmat(w2,1,2^n)]; % Weights
    N = length(wi);
end 