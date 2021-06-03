% Generates a white noise IC

function [u0, v0] = ICRand(nx,ny,dx,dy,lxmin,lxmax,lymin,lymax)

    Mu = ones(nx,ny);
    Mv = ones(nx,ny);
    
    % Initial condition
    u0 = Mu.*noise2(nx,ny,dx,dy,lxmin,lxmax,lymin,lymax);
    v0 = Mv.*noise2(nx,ny,dx,dy,lxmin,lxmax,lymin,lymax);
    
    % Sym
    %u0 = (Yu < 0).*u0;
    %v0 = (Yv < 0).*v0;

end