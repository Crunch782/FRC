% Main function for Evaluating the Objective Function with DAL


function [J, dJ]=J_eval(Qg)

    param;

    % Convert to fields
    u = reshape(Qg(1:nxm*nym),nxm,nym);
    v = reshape(Qg(nxm*nym+1:end),nxm,nym);
    J = 0;
    dJ = [];

    % Set the parameters
    param;

    % Set the Passive Scalar
    c = ICScalar(xm,ym,dx,dy,Vol);

    % Now initialize Pressure
    p = ICPress(u,v,dx,dy,im,ic,ip,jm,jc,jp,al,cl,alpl,xs2l);

    % Check Normalization
    fprintf('\nInitial KE = %.2f\n',0.5*(L2Norm(u,1).^2 + L2Norm(v,1).^2));

    % Main Loop to integrate from t=0 to t=T
    [u, v, p, c] = NSE_Solver(u,v,p,c);

    % Compute the Terminal Conditions for the Adjoint Eqns
    [u, v, p, c, J] = Terminal(u,v,p,c,s,xm,ym,nxm,nym,dx,dy);

    % Main Loop to integrate from t=T to t=0
    [u, v, ~, ~] = ADJNSE_Solver(u,v,p,c);
    
    % Compute Grad
    dJ = [u(:) ; v(:)];
    
end


