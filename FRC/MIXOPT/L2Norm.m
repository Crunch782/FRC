% Computes the L2Norm 

function I = L2Norm(u, ave)

    param;
    
    I = sum(sum(u.^2))*(dx*dy);
    
    if ave == 1
        I = I./Vol;
    end
        
    I = sqrt(I);
        
end