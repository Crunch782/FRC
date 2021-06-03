% Ordered fft2
% f: inverse fourier transform of g
% x: vector of space along k
% y: vector of space along h


function [f, X, Y] = oifft2(g,k,h,xx,yy)

    [nh, nk] = size(k);
    dk = diff(k(1,:)); dk=dk(1,1);
    dh = diff(h(:,1)); dh=dh(1,1);

    [X, Y] = meshgrid(xx,yy);
    X=X';
    Y=Y';

    g = fftshift(g);
    f = (nk*nh)*ifft2(g);

end