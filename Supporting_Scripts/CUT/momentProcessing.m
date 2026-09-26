function [mu,P,gamma,kurt] = momentProcessing(wi, xi)
    %% 
    % Processes Cubature Points From 
    % N. Adurthi, P. Singla, and T. Singh, “Conjugate unscented transfor-
    % mation: Applications to estimation and control,” Journal of Dynamic
    % Systems, Measurement, and Control, vol. 140, no. 3, pp. 1–22, 2018
    %
    % Inputs: Cubature points (xi), weights (wi)
    % Outputs: Mean (mu), Covariance (P), Skewness (gamma), kurtosis (kurt)
    %%

    [n,~] = size(xi);

    % Increase Precision with VPA
    wi = vpa(wi);
    xi = vpa(xi);

    % Diagonal Moments
    Ex = sym(zeros(n,4));
    for m = 1:4
        Ex(:,m) = sum(wi.*(xi.^m),2);
    end

    % Mean
    mu = sum(wi.*xi,2);    

    % Covariance 
    Exx = zeros(n);
    for i = 1:length(wi)
        Exx = Exx + wi(i)*(xi(:,i)*xi(:,i).');
    end
    P = Exx - mu*mu.';
    sigma = sqrt(diag(P));

    % Skewness (Diagonal Only)
    %gamma = (Ex(:,3) - 3*mu.*sigma.^2 - mu.^3)./(sigma.^3);
    gamma = (Ex(:,3) - 3*mu.*Ex(:,2) + 2*mu.^3)./(sigma.^3);

    % Kurtosis (Diagonal Only)
    kurt = (Ex(:,4) - 4*Ex(:,3).*mu + 6*Ex(:,2).*mu.^2 - 3*mu.^4)./(sigma.^4);

    mu = double(mu);
    P = double(P);
    gamma = double(gamma);
    kurt = double(kurt);
end