% Calculates the Nonlinear term in the momentum eqn

function [hu, hv] = nonlinearADJ(u,v,ud,vd,dx,dy,im,ic,ip,jm,jc,jp)

    % Calculate Hu
    hu = (0.25./dx)*((ud+ud(ip,jc)).*(u+u(ip,jc)) - (ud+ud(im,jc)).*(u+u(im,jc))); % duudx
    hu = hu + (0.25./dy)*((u+u(ic,jp)).*(vd(ic,jp)+vd(im,jp)) - (u+u(ic,jm)).*(vd+vd(im,jc)));
    hu = (-1)*hu;
    
    % Calculate Hv
    hv = (0.25./dy)*((vd+vd(ic,jp)).*(v+v(ic,jp)) - (vd+vd(ic,jm)).*(v+v(ic,jm))); % duudx
    hv = hv + (0.25./dx)*((ud(ip,jc)+ud(ip,jm)).*(v+v(ip,jc)) - (ud+ud(ic,jm)).*(v+v(im,jc)));
    hv = (-1)*hv;
    
end