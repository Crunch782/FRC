% Computes the Terminal Conditions and the Objective Functional

function [uT, vT, pT, cT, J] = Terminal(uT,vT,pT,cT,s,x,y,nx,ny,dx,dy)

    % First, compute the Objective Functional
    ck = zeros(nx,ny);
    gammak = zeros(nx,ny);
    J = 0;
    [ck, K, L] = offt2(cT,x,y,nx,ny,dx,dy);
    kappa2 = K.^2 + L.^2;
    kappa2 = kappa2.^s;
    gammak = ck;
    ck = abs(ck).^2;
    ck = ck./kappa2;
    ck(nx./2+1,ny./2+1) = 0;
    J = sum(sum(ck));
    J = J./2;
    
    % Now compute the adjoint terminal scalar
    gammak = gammak./kappa2;
    gammak(nx./2+1,ny./2+1) = 0;
    
    [cT, ~, ~] = oifft2(gammak,K,L,x,y);
    cT = real(cT);
    
    % and the other adjoint variables
    uT = zeros(nx,ny);
    vT = zeros(nx,ny);
    pT = zeros(nx,ny);

end