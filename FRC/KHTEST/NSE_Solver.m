% Function for Solving NSE

function [u, v, p, c] = NSE_Solver(u,v,p,c,dt,maxits,ax,cx,alpx,xs2x,ay,cy,alpy,xs2y,al,cl,alpl,xs2l)

    param;

    % Initialize local variables
    px = zeros(nxm,nym);
    py = zeros(nxm,nym);
    rhsu = zeros(nxm,nym);
    rhsv = zeros(nxm,nym);
    rhsc = zeros(nxm,nym);
    du = zeros(nxm,nym);
    dv = zeros(nxm,nym);
    du1 = zeros(nxm,nym);
    dv1 = zeros(nxm,nym);
    dc = zeros(nxm,nym);
    dc1 = zeros(nxm,nym);
    Q = zeros(nxm,nym);
    Qb = zeros(nxm,nym);
    phi = zeros(nxm,nym);
    phib = zeros(nxm,nym);
    
    % time
    t = 0;
    
    % Data
    K = [];
    tt = [];
    I = [];
    Kc = [];
    Ic = [];
    
    % File for Flow Data
    %filename = 'RK4_KHIC_HR.txt';
    %fopen(filename,'w+');
    
    % Calculate the initial convection term
    %[hu, hv] = nonlinear(u,v,dx,dy,im,ic,ip,jm,jc,jp);
    %hc = nonlinearAdv(u,v,c,dx,dy,im,ic,ip,jm,jc,jp);
    
    axis tight manual
    AX = gca;
    AX.NextPlot = 'replaceChildren';
    F(maxits) = struct('cdata',[],'colormap',[]);
    
    iters = 0;
    while(iters < maxits)
       
        % Compute the RHS for the ADI method
        %{
        rhsu = dt*(-px - (1./2)*hu + (1./Re)*lapl(u,im,ic,ip,jm,jc,jp,dx,dy));
        rhsv = dt*(-py - (1./2)*hv + (1./Re)*lapl(v,im,ic,ip,jm,jc,jp,dx,dy));
        
        [hu, hv] = nonlinear(u,v,dx,dy,im,ic,ip,jm,jc,jp);
        rhsu = rhsu + dt*(3./2)*hu;
        rhsv = rhsv + dt*(3./2)*hv;
        %}
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
        KE = 0.5*(L2Norm(u,1).^2 + L2Norm(v,1).^2);
        I1 = sum(sum(lapl(u,im,ic,ip,jm,jc,jp,dx,dy).*u + lapl(v,im,ic,ip,jm,jc,jp,dx,dy).*v));
        I1 = I1*(dx*dy)./(Re*Vol);
        I2 = (1./Vol)*sum(sum(au.*u + av.*v))*(dx*dy);
        I(iters+1) = I1+I2;
        K(iters+1) = KE;
        tt(iters+1) = t;
        
        % Now we solve the equation for the passive scalar
        %{
        rhsc = dt*((1./2)*hc + (1./Pe)*lapl(c,im,ic,ip,jm,jc,jp,dx,dy));
        hc = nonlinearAdv(u,v,c,dx,dy,im,ic,ip,jm,jc,jp);
        rhsc = rhsc + dt*(3./2)*hc;
        %}
        ac1 = nonlinearAdv(u,v,c,dx,dy,im,ic,ip,jm,jc,jp);
        ac2 = nonlinearAdv(u,v,c+(dt./2)*ac1,dx,dy,im,ic,ip,jm,jc,jp);
        ac3 = nonlinearAdv(u,v,c+(dt./2)*ac2,dx,dy,im,ic,ip,jm,jc,jp);
        ac4 = nonlinearAdv(u,v,c+dt*ac3,dx,dy,im,ic,ip,jm,jc,jp);
        ac = (1./6)*(ac1+2*ac2+2*ac3+ac4);
        rhsc = dt*(ac + (1./Pe)*lapl(c,im,ic,ip,jm,jc,jp,dx,dy));
              
        % Preform the ADI method to solve the momentum eqns for u*,v*
        dc1 = ADI_step(ax,cx,alpx,xs2x,rhsc');
        dc = ADI_step(ay,cy,alpy,xs2y,dc1');
        c = c + dc;
        
        % Data for Integral Eqn
        KE = 0.5*L2Norm(c,1).^2;
        I1 = sum(sum(lapl(c,im,ic,ip,jm,jc,jp,dx,dy).*c));
        I1 = I1*(dx*dy)./(Pe*Vol);
        I2 = (1./Vol)*sum(sum(ac.*c))*(dx*dy);
        Ic(iters+1) = I1+I2;
        Kc(iters+1) = KE;
        
                        
        % Update time and iterations
        t = t + dt;        
        iters = iters+1;
        if mod(iters, 50) == 0
            fprintf('\n t == %.3f ',t);
            D = (u(ip,jc)-u)./dx + (v(ic,jp)-v)./dy;
            div = L2Norm(D,0);
            fprintf('____ divergence == %.14f',div);
            %{
            KU = L2Norm(u,1).^2;
            KV = L2Norm(v,1).^2;
            KC = L2Norm(c,1).^2;
            D = (u(ip,jc)-u)./dx + (v(ic,jp)-v)./dy;
            div = L2Norm(D,0);
            flowdata = [t KU KV KC div];
            save(filename,'flowdata','-ASCII','-append');
            %}
            
        end
        
    end
    
    
    dkcdt = diff(Kc)./dt;
    dkdt = diff(K)./dt;
    tt = tt + dt./2;
    Icnew = zeros(1, length(Ic)-1);
    for i = 1:length(Ic)-1
        Icnew(i) = (Ic(i)+Ic(i+1))./2;
    end
    
    Inew = zeros(1, length(I)-1);
    for i = 1:length(I)-1
        Inew(i) = (I(i)+I(i+1))./2;
    end
    
    tt(end) = [];
    
    errc = (dt./(T-dt))*sum(sum(abs(dkcdt - Icnew).^2));
    errc = sqrt(errc)
    
    err = (dt./(T-dt))*sum(sum(abs(dkdt - Inew).^2));
    err = sqrt(err)
    
    plot(tt,dkcdt,'-r')
    hold on
    plot(tt,Icnew,'-b')
    plot(tt,dkdt, '-g')
    plot(tt,Inew, '-y')
    xlim([dt./2 T-dt./2])
    ylim([min(min(Ic)) max(max(dkdt))])
    xlabel('t')
    legend('dkcdt','Ic','dkdt','I')
    
    
end