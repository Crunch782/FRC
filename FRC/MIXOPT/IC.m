% Sets the initial condition for the KH test

function u0 = IC(xc,ym,U0,Pj,Ly,Rj,Ax,lx)
    
    % Define grid and u1,u2
    [X,Y] = meshgrid(xc,ym);
    X = X'; 
    Y = Y';
    
    u1 = (U0./2)*(1 + tanh((Pj./2)*(1 - abs(Ly./2 - Y)./Rj)));
    u2 = Ax*sin(2*pi*X./lx);
    
    u0 = u1.*(1+u2);
    

end