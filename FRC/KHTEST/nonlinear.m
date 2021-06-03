% Calculates the Nonlinear term in the momentum eqn

function [hu, hv] = nonlinear(u,v,dx,dy,im,ic,ip,jm,jc,jp)

    % Calculate Hu
    hu = (0.25./dx)*((u+u(ip,jc)).^2 - (u+u(im,jc)).^2); % duudx
    hu = hu + (0.25./dy)*((u+u(ic,jp)).*(v(ic,jp)+v(im,jp)) - (u+u(ic,jm)).*(v+v(im,jc)));
    hu = (-1)*hu;
    
    % Calculate Hv
    hv = (0.25./dy)*((v+v(ic,jp)).^2 - (v+v(ic,jm)).^2); % duudx
    hv = hv + (0.25./dx)*((u(ip,jc)+u(ip,jm)).*(v+v(ip,jc)) - (u+u(ic,jm)).*(v+v(im,jc)));
    hv = (-1)*hv;
    
end