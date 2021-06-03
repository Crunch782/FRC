% Compute dt with CFL condition

function dt = compute_dt()

    param;
    dt = (cfl./2)./(1./(dx.^2) + 1./(dy.^2));

end