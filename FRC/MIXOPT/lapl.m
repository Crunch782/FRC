% Calculates the laplacian

function lp = lapl(u,im,ic,ip,jm,jc,jp,dx,dy)

    lp = (u(ip, jc)-2*u+u(im,jc))./(dx.^2) + (u(ic,jp)-2*u+u(ic,jm))./(dy.^2);

end