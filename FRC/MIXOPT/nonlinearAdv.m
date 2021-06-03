% Calculates the Nonlinear term in the momentum eqn

function hc = nonlinearAdv(u,v,c,dx,dy,im,ic,ip,jm,jc,jp)

    % Calculate Hc
    hc = (0.5./dx)*((u(ip, jc).*(c(ip,jc) + c)) - u.*(c(im,jc) + c));
    hc = hc + (0.5./dy)*((v(ic, jp).*(c(ic,jp) + c)) - v.*(c(ic,jm) + c));
    hc = hc*(-1);
    
end