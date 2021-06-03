% Function for Solving NSE

function [u, v, p, c] = NSE_SolverDiag(u,v,p,c,cd,vl,vu)

    param;

    % Initialize local variables
    px = zeros(nxm,nym);
    py = zeros(nxm,nym);
    rhsu = zeros(nxm,nym);
    rhsv = zeros(nxm,nym);
    rhsc = zeros(nxm,nym);
    rhscd = zeros(nxm,nym);
    du = zeros(nxm,nym);
    dv = zeros(nxm,nym);
    du1 = zeros(nxm,nym);
    dv1 = zeros(nxm,nym);
    dc = zeros(nxm,nym);
    dc1 = zeros(nxm,nym);
    dcd = zeros(nxm,nym);
    dcd1 = zeros(nxm,nym);
    Q = zeros(nxm,nym);
    Qb = zeros(nxm,nym);
    phi = zeros(nxm,nym);
    phib = zeros(nxm,nym);
    
    % time
    t = 0;
    iters = 0;
    
    % Data Arrays   
    nff = 1;
    t_array = zeros(1,nf+1);
    V_array = zeros(1,nf+1);
    M_array = zeros(1,nf+1);
    
    fprintf('\n Writing Data at time %.3f ',t);
    [V,M] = fdata(c,cd,u,v,xc,yc,xm,ym,nxm,nym,dx,dy,im,jc,ic,jm,t,tag,vl,vu,s,nff);
    t_array(nff) = 0;
    V_array(nff) = 1;
    M_array(nff) = 1;
    nff = nff + 1;
       
    while(iters < nd)
       
        % Compute the RHS for the ADI method
        
        [au1, av1] = nonlinear(u,v,dx,dy,im,ic,ip,jm,jc,jp);
        [au2, av2] = nonlinear(u+(dt./2)*au1,v+(dt./2)*av1,dx,dy,im,ic,ip,jm,jc,jp);
        [au3, av3] = nonlinear(u+(dt./2)*au2,v+(dt./2)*av2,dx,dy,im,ic,ip,jm,jc,jp);
        [au4, av4] = nonlinear(u+dt*au3,v+dt*av3,dx,dy,im,ic,ip,jm,jc,jp);
        au = (1./6)*(au1+2*au2+2*au3+au4);
        av = (1./6)*(av1+2*av2+2*av3+av4);
        
        rhsu = dt*(-px + au + (1./Re)*lapl(u,im,ic,ip,jm,jc,jp,dx,dy));
        rhsv = dt*(-py + av + (1./Re)*lapl(v,im,ic,ip,jm,jc,jp,dx,dy));
        
        % Preform the ADI method to solve the momentum eqns for u*,v*
        du1 = ADI_step(ax,cx,alpx,xs2x,rhsu');
        du = ADI_step(ay,cy,alpy,xs2y,du1');
        u = u + du;
        
        dv1 = ADI_step(ax,cx,alpx,xs2x,rhsv');
        dv = ADI_step(ay,cy,alpy,xs2y,dv1');
        v = v + dv;
        
        % Solve the Poisson Eqn for the corrector
        Q = (1./dt)*( (u(ip,jc)-u)./dx  +  (v(ic,jp)-v)./dy);
        Qb = fft(Q);
        Qb(1,1) = 0;
        phib = POI_step(al,cl,alpl,xs2l,Qb);
        phi = ifft(phib);
        phi = real(phi);
        
        % Correct the velocity field
        u = u - (dt./dx)*(phi-phi(im,jc));
        v = v - (dt./dy)*(phi-phi(ic,jm));
        
        % Compute the new pressure
        p = p + phi - (dt./(2*Re))*lapl(phi,im,ic,ip,jm,jc,jp,dx,dy);
        
        % Compute the new pressure gradient
        px = (p - p(im,jc))./dx;
        py = (p - p(ic,jm))./dy;
        
        % Data for Integral Eqn
        
        ac1 = nonlinearAdv(u,v,c,dx,dy,im,ic,ip,jm,jc,jp);
        ac2 = nonlinearAdv(u,v,c+(dt./2)*ac1,dx,dy,im,ic,ip,jm,jc,jp);
        ac3 = nonlinearAdv(u,v,c+(dt./2)*ac2,dx,dy,im,ic,ip,jm,jc,jp);
        ac4 = nonlinearAdv(u,v,c+dt*ac3,dx,dy,im,ic,ip,jm,jc,jp);
        ac = (1./6)*(ac1+2*ac2+2*ac3+ac4);
        rhsc = dt*(ac + (1./Pe)*lapl(c,im,ic,ip,jm,jc,jp,dx,dy));
        rhscd = dt*((1./Pe)*lapl(cd,im,ic,ip,jm,jc,jp,dx,dy));
              
        % Preform the ADI method to solve the momentum eqns for u*,v*
        dc1 = ADI_step(ax,cx,alpx,xs2x,rhsc');
        dc = ADI_step(ay,cy,alpy,xs2y,dc1');
        c = c + dc;

        dcd1 = ADI_step(ax,cx,alpx,xs2x,rhscd');
        dcd = ADI_step(ay,cy,alpy,xs2y,dcd1');
        cd = cd + dcd;
                          
        % Update time and iterations
        t = t + dt;        
        iters = iters+1;

        % Create and store the Data
        if mod(iters, nf) == 0
            fprintf('\n Writing Data at time %.3f ',t);
            [V,M] = fdata(c,cd,u,v,xc,yc,xm,ym,nxm,nym,dx,dy,im,jc,ic,jm,t,tag,vl,vu,s,nff);
            t_array(nff) = t;
            V_array(nff) = V;
            M_array(nff) = M;
            nff = nff + 1;
        end
        
    end
    
    %Save arrays
    save(strcat('Data/',tag,'/V_data.dat'),'V_array','-ASCII');
    save(strcat('Data/',tag,'/M_data.dat'),'M_array','-ASCII');
    save(strcat('Data/',tag,'/t_data.dat'),'t_array','-ASCII');
    
end