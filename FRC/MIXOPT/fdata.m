% Computes and saves the diagnostic data

function [V,M] = fdata(c,cd,u,v,xc,yc,xm,ym,nx,ny,dx,dy,im,jc,ic,jm,t,dirname,vl,vu,s,nff,T,Re)

    % Compute the Variance
    V1 = L2Norm(c,0).^2;
    V2 = L2Norm(cd,0).^2;
    V = V1./V2;
    
    %Compute the Mix-Norm
    [ck, K, L] = offt2(c,xm,ym,nx,ny,dx,dy);
    kappa2 = K.^2 + L.^2;
    kappa2 = kappa2.^s;
    ck = abs(ck).^2;
    ck = ck./kappa2;
    ck(nx./2+1,ny./2+1) = 0;
    M1 = sum(sum(ck));
    
    [cdk, ~, ~] = offt2(cd,xm,ym,nx,ny,dx,dy);
    cdk = abs(cdk).^2;
    cdk = cdk./kappa2;
    cdk(nx./2+1,ny./2+1) = 0;
    M2 = sum(sum(cdk));
    
    M = M1./M2;
    
    % Now compute the vorticity + scalar plots
    w = vort(u,v,dx,dy,im,jc,ic,jm);
    PlotFields(c,xm,ym,w,xc,yc,tag,t,vl,vu,nff,s,T,Re);
    

end