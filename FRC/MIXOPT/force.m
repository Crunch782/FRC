% Computes the Forcing term for the Adjoint Eqns

function [fx, fy] = force(c,cd,ic,jc,im,jm,dx,dy)

    fx = c.*(cd-cd(im,jc))./dx;
    fy = c.*(cd-cd(ic,jm))./dy;

end