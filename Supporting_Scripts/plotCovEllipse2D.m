function [X,Y] = plotCovEllipse2D(P,nSigma)
    % Plot 3D Covariance Ellipse with n-Sigma

    % Position Only (2 x 2)
    [eigVec,eigVal] = eig(P); 

    % Rotation
    vec1 = eigVec(:,1);
    vec2 = eigVec(:,2);
    
    NB = [vec1,vec2]; 

    % Standard Deviations in Ellipsoidal Frame
    nStd = nSigma*sqrt(diag(eigVal));

    theta = linspace(0,2*pi,1000);
    X_E = nStd(1)*cos(theta);
    Y_E = nStd(2)*sin(theta);

    r_E = [X_E;Y_E];

    r = NB*r_E;
    X = r(1,:);
    Y = r(2,:);
end