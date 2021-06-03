% Generates a white noise IC

function [u0, v0] = ICRand(nx,ny,xc,ym,xm,yc,dx,dy,lxmin,lxmax,lymin,lymax)

    Mu = ones(nx,ny);
    Mv = ones(nx,ny);
    
    % Initial condition
    u0 = Mu.*noise2(nx,ny,xc,ym,dx,dy,lxmin,lxmax,lymin,lymax);
    v0 = Mv.*noise2(nx,ny,xm,yc,dx,dy,lxmin,lxmax,lymin,lymax);

end