% Performs the  adi method

function f = ADI_step(aa, ac, alpha, xs2, f)

    param;
    
    [m,n] = size(f);
    
    for k = 2:n
        f(:,k) = f(:,k)-aa(k)*f(:,k-1); 
    end
    
    f(:,n) = f(:,n)*alpha(n);
    
    for k=n-1:-1:1
        f(:,k) = f(:,k)*alpha(k)-ac(k)*f(:,k+1); 
    end
    
    xs1 = zeros(m,1);
    xs1 = (f(:,1) + f(:,n))./(1 + xs2(:,1) + xs2(:, n));
    
    for k=1:n
        f(:,k) = f(:,k) - xs1*xs2(k); 
    end

end