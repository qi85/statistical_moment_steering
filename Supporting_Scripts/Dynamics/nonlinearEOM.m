function Xdot = nonlinearEOM(t,X,params,flag_linearized)
    % Default behavior if flag not provided or empty
    if nargin < 4 || isempty(flag_linearized)
        flag_linearized = 0;
    end

    % Unpack Parameters and State Variables
    alpha = params.alpha;
    beta = params.beta;
    kappa = params.kappa;

    x = X(1);
    y = X(2);
    xdot = X(3);
    ydot = X(4);
    
    % EOM
    xddot = -alpha*x-beta*x^3-kappa*(x-y);
    yddot = -alpha*y-beta*y^3-kappa*(y-x);

    Xdot = [xdot;ydot;xddot;yddot];

    if flag_linearized
        % Linearize System
        n = 4;

        STM = reshape(X(n+1:n+n^2),n,n);

        A = [zeros(2),  eye(2);
             -alpha-3*beta*x^2-kappa, kappa, 0,0;
             kappa, -alpha-3*beta*y^2-kappa, 0,0];

        STMdot = A*STM;

        c = [xdot;ydot;xddot;yddot] - A*[x;y;xdot;ydot];
        cdot = STM\c;

        Xdot = [Xdot;reshape(STMdot,n^2,1);cdot];
    end
end