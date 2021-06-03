% Computes the RHS for the ADI method

function rhs = RHS(u, f, dt)

    param;
    rhs = zeros(nxm, nym);
    
    rhs = NSE_calc_lap(u);
    rhs = rhs + f;
    rhs = rhs*dt;
    
end