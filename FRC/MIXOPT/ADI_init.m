% Inializes the ADI method

function [aa, ac, alpha, xs2] = ADI_init(aa, ab, ac)

    n = length(aa);
    
    % Modify to get b*
    ab(1) = ab(1)-aa(1);
    ab(n) = ab(n)-ac(n);
    
    % Initialize X2
    xs2 = zeros(1,n);
    xs2(1) = aa(1);
    xs2(n) = ac(n);
    
    % Compute coeffs for M*X1 eqn
    alpha = zeros(1,n);
    alpha(1) = 1./ab(1);
    for k=2:n
        alpha(k) =1./(ab(k)-(ac(k-1)*alpha(k-1))*aa(k)); 
        aa(k) = aa(k)*alpha(k-1);
        ac(k-1) = ac(k-1)*alpha(k-1);
    end
    
    ab(1) = xs2(1);
    for k=2:n
        ab(k) = xs2(k) - aa(k)*ab(k-1);
    end
    
    % Solve for X2
    xs2(n) = ab(n)*alpha(n);
    for k=n-1:-1:1
        xs2(k) = ab(k)*alpha(k)-ac(k)*xs2(k+1);
    end

end