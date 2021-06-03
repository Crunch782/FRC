% Calculates the vorticity field

function w = vort(u,v,dx,dy,im,jc,ic,jm)

    w = (v-v(im,jc))./dx - (u - u(ic,jm))./dy;  

end