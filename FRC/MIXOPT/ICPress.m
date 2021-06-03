% Computes the inital pressure

function p0 = ICPress(u,v,dx,dy,im,ic,ip,jm,jc,jp,al,cl,alpl,xs2l)

    [hu, hv] = nonlinear(u,v,dx,dy,im,ic,ip,jm,jc,jp);
    divUV = (hu(ip,jc)-hu)./dx + (hv(ic,jp)-hv)./dy;
    
    Qb = fft(divUV);
    Qb(1,1) = 0;
    phib = POI_step(al,cl,alpl,xs2l,Qb);
    p0 = ifft(phib);
    p0 = real(p0);

end