function [Az, Amu] = linearize_z(wi, nx, ns)
    %% Selector Matrices
    Ei = @(i,nx,ns) [zeros(nx,nx*(i-1)),eye(nx),zeros(nx,nx*(ns-i))]; 
    Ibar = repmat(eye(nx),ns,1);

    %% Mean
    Amu = zeros(nx, nx*ns);
    for i = 1:ns
        Amu = Amu + wi(i)*Ei(i,nx,ns);
    end
    Az = eye(nx*ns) - Ibar*Amu;
end