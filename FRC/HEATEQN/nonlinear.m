% Computes the Nonlinear advection term

function nl = nonlinear(u)

    param;
    
    nl = -((u + u(ip,jc)).^2 - (u + u(im,jc)).^2)./(4*dx);

end