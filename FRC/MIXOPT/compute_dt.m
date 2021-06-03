% Compute the timestep in accordance with the CFL condition

function dt = compute_dt(u,v,dx,dy,cfl)

    umax = max(max(abs(u)))./dx;
    vmax = max(max(abs(v)))./dy;
    
    dt = cfl./(umax + vmax);

end