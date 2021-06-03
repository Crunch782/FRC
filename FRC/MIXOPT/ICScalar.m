% Sets the IC for the Passive Scalar

function c0 = ICScalar(x,y,dx,dy,Vol)

    [X,Y] = meshgrid(x,y);
    X=X';
    Y=Y';
    
    c0 = tanh(6*(X-pi./2))-tanh(6*(X-(3./2)*pi)) - 1;
    
    meanc0 = sum(sum(c0))*(dx*dy)./Vol;
    while meanc0 > 1e-14
        c0 = c0 - meanc0;
        meanc0 = sum(sum(c0))*(dx*dy)./Vol        
    end

end