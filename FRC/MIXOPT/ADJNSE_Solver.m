% Function for Solving NSE

function [u, v, p, c] = ADJNSE_Solver(u,v,p,c)

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
    
    % Initialize Storage for Direct variables
    un = zeros(nxm,nym);
    vn = zeros(nxm,nym);
    cn = zeros(nxm,nym);
    un1 = zeros(nxm,nym);
    vn1 = zeros(nxm,nym);
    cn1 = zeros(nxm,nym);
    unh = zeros(nxm,nym);
    vnh = zeros(nxm,nym);
    cnh = zeros(nxm,nym);
    
    % File for Flow Data
    %{
    filename = 'adjoint.txt';
    fopen(filename,'w+');
    %}
    
    % time
    t = T;
    
    % Data
    %{
    K = [];
    tt = [];
    I = [];
    Kc = [];
    Ic = [];
    %}
    
    %Save Initial fields
    iters = nt;
    
    while(iters > 0)
        
        % Load the Checkpoints
        % First step
        un = load(strcat('Checkpoints/',tag,'/u_',num2str(iters),'.dat'),'-ASCII');
        vn = load(strcat('Checkpoints/',tag,'/v_',num2str(iters),'.dat'),'-ASCII');
        cn = load(strcat('Checkpoints/',tag,'/c_',num2str(iters),'.dat'),'-ASCII');
        
        % Second step
        un1 = load(strcat('Checkpoints/',tag,'/u_',num2str(iters-1),'.dat'),'-ASCII');
        vn1 = load(strcat('Checkpoints/',tag,'/v_',num2str(iters-1),'.dat'),'-ASCII');
        cn1 = load(strcat('Checkpoints/',tag,'/c_',num2str(iters-1),'.dat'),'-ASCII');
        
        % Intermediate
        unh = (un+un1)./2;
        vnh = (vn+vn1)./2;
        cnh = (cn+cn1)./2;
        
        % Solve the Diffusion Eqn
        ac1 = nonlinearAdv(un,vn,c,dx,dy,im,ic,ip,jm,jc,jp);
        ac2 = nonlinearAdv(unh,vnh,c+(dt./2)*ac1,dx,dy,im,ic,ip,jm,jc,jp);
        ac3 = nonlinearAdv(unh,vnh,c+(dt./2)*ac2,dx,dy,im,ic,ip,jm,jc,jp);
        ac4 = nonlinearAdv(un1,vn1,c+dt*ac3,dx,dy,im,ic,ip,jm,jc,jp);
        ac = (1./6)*(ac1+2*ac2+2*ac3+ac4);
        rhsc = dt*(-ac + (1./Pe)*lapl(c,im,ic,ip,jm,jc,jp,dx,dy));
              
        % Preform the ADI method to solve the momentum eqns for u*,v*
        dc1 = ADI_step(ax,cx,alpx,xs2x,rhsc');
        dc = ADI_step(ay,cy,alpy,xs2y,dc1');
        c = c + dc;
        
        % Data for Integral Eqn
        %{
        KE = 0.5*L2Norm(c,1).^2;
        I1 = sum(sum(lapl(c,im,ic,ip,jm,jc,jp,dx,dy).*c));
        I1 = I1*(dx*dy)./(Pe*Vol);
        I2 = (-1./Vol)*sum(sum(ac.*c))*(dx*dy);
        Ic(iters+1) = I1+I2;
        Kc(iters+1) = KE;
        %}
        
        % Compute Advection
        [au1, av1] = nonlinearADJ(u,v,un,vn,dx,dy,im,ic,ip,jm,jc,jp);
        [au2, av2] = nonlinearADJ(u+(dt./2)*au1,v+(dt./2)*av1,unh,vnh,dx,dy,im,ic,ip,jm,jc,jp);
        [au3, av3] = nonlinearADJ(u+(dt./2)*au2,v+(dt./2)*av2,unh,vnh,dx,dy,im,ic,ip,jm,jc,jp);
        [au4, av4] = nonlinearADJ(u+dt*au3,v+dt*av3,un1,vn1,dx,dy,im,ic,ip,jm,jc,jp);
        au = (1./6)*(au1+2*au2+2*au3+au4);
        av = (1./6)*(av1+2*av2+2*av3+av4);
        
        % Compute the Forcing Term
        [cforu1, cforv1] = force(c,cn,ic,jc,im,jm,dx,dy);
        [cforu2, cforv2] = force(c,cnh,ic,jc,im,jm,dx,dy);
        [cforu3, cforv3] = force(c,cnh,ic,jc,im,jm,dx,dy);
        [cforu4, cforv4] = force(c,cn1,ic,jc,im,jm,dx,dy);
        cforu = 1/6*dt*(cforu1+2*cforu2+2*cforu3+cforu4);
        cforv = 1/6*dt*(cforv1+2*cforv2+2*cforv3+cforv4);
        
        % Compute the RHS
        rhsu = dt*(px - au - cforu + (1./Re)*lapl(u,im,ic,ip,jm,jc,jp,dx,dy));
        rhsv = dt*(py - av - cforv + (1./Re)*lapl(v,im,ic,ip,jm,jc,jp,dx,dy));
        
        % Preform the ADI method to solve the momentum eqns for u*,v*
        du1 = ADI_step(ax,cx,alpx,xs2x,rhsu');
        du = ADI_step(ay,cy,alpy,xs2y,du1');
        u = u + du;
        
        dv1 = ADI_step(ax,cx,alpx,xs2x,rhsv');
        dv = ADI_step(ay,cy,alpy,xs2y,dv1');
        v = v + dv;
        
        % Solve the Poisson Eqn for the corrector
        Q = -(1./dt)*( (u(ip,jc)-u)./dx  +  (v(ic,jp)-v)./dy);
        Qb = fft(Q);
        Qb(1,1) = 0;
        phib = POI_step(al,cl,alpl,xs2l,Qb);
        phi = ifft(phib);
        phi = real(phi);
        
        % Correct the velocity field
        u = u + (dt./dx)*(phi-phi(im,jc));
        v = v + (dt./dy)*(phi-phi(ic,jm));
        
        % Compute the new pressure
        p = p + phi - (dt./(2*Re))*lapl(phi,im,ic,ip,jm,jc,jp,dx,dy);
        
        % Compute the new pressure gradient
        px = (p - p(im,jc))./dx;
        py = (p - p(ic,jm))./dy;
        
        % Data for Integral Eqn
        %{
        KE = 0.5*(L2Norm(u,1).^2 + L2Norm(v,1).^2);
        I1 = sum(sum(lapl(u,im,ic,ip,jm,jc,jp,dx,dy).*u + lapl(v,im,ic,ip,jm,jc,jp,dx,dy).*v));
        I1 = I1*(dx*dy)./(Re*Vol);
        I2 = (-1./Vol)*sum(sum(au.*u + av.*v))*(dx*dy);
        I3 = (1./Vol)*sum(sum(cforu.*u+cforv.*v))*dx*dy;
        I(iters+1) = I1+I2-I3;
        K(iters+1) = KE;
        tt(iters+1) = t;
        %}
        
        % Update time and iterations
        t = t - dt;        
        iters = iters-1;
        
        %{
        if mod(iters, 50) == 0         
            fprintf('\n t == %.3f',t);
            KU = L2Norm(u,1).^2;
            KV = L2Norm(v,1).^2;
            KC = L2Norm(c,1).^2;
            D = (u(ip,jc)-u)./dx + (v(ic,jp)-v)./dy;
            div = L2Norm(D,0);
            flowdata = [t KU KV KC div];
            save(filename,'flowdata','-ASCII','-append');
        end
        %}
        
   
    end
    
    %{
    dkcdt = diff(Kc)./dt;
    dkdt = diff(K)./dt;
    
    Icnew = zeros(1, length(Ic)-1);
    for i = 1:length(Ic)-1
        Icnew(i) = (Ic(i)+Ic(i+1))./2;
    end
    
    Inew = zeros(1, length(I)-1);
    for i = 1:length(I)-1
        Inew(i) = (I(i)+I(i+1))./2;
    end
    
    ttnew = zeros(1, length(tt)-1);
    for i = 1:length(I)-1
        ttnew(i) = (tt(i)+tt(i+1))./2;
    end
    
    % since dtau = -dt
    dkdt = dkdt*(-1);
    dkcdt = dkcdt*(-1);
    
    errc = (dt./(T-2*dt))*sum(sum(abs(dkcdt - Icnew).^2));
    errc = sqrt(errc)
    
    err = (dt./(T-2*dt))*sum(sum(abs(dkdt - Inew).^2));
    err = sqrt(err)
    
    save(strcat('dkc.dat'),'dkcdt','-ASCII');
    save(strcat('Ic.dat'),'Icnew','-ASCII');
    save(strcat('dk.dat'),'dkdt','-ASCII');
    save(strcat('I.dat'),'Inew','-ASCII');
    save(strcat('tt.dat'),'tt','-ASCII');
    
    
    dkdt = real(log(dkdt));
    Inew = real(log(Inew));
    plot(ttnew,dkcdt,'-r')
    hold on
    plot(ttnew,Icnew,'--b')
    plot(ttnew,dkdt, '-g')
    plot(ttnew,Inew, '--y')
    xlim([2*dt T-dt./2])
    ylim([-30 1])
    xlabel('t')
    legend('dkcdt','Ic','dkdt','I')
    hold off
    %}
    
    
    
end
