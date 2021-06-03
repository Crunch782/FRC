% Calculates the Laplacian of function on the Grid

function lapu = NSE_calc_lap(u)

    param;
    
    lapu = (u(ip, jc)-2*u+u(im,jc))./(dx.^2) + (u(ic,jp)-2*u+u(ic,jm))./(dy.^2);

end
