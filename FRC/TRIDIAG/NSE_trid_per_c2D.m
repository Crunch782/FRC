% Solves a Periodic Tridiagnonal Matrix System Using the Thomas Algorithm

function fi = NSE_trid_per_c2D(aa, ab, ac, fi)

    [m, n] = size(ab);

    % Solve the periodic trid diag
    
    v      = zeros(m,n);
    v(:,1)   = aa(:,1);
    v(:,end) = ac(:,end);
    
    ab(:,1)   = ab(:,1)-aa(:,1);
    ab(:,end) = ab(:,n)-ac(:,n);
    
    X1 = NSE_trid_c2D(aa, ab, ac, fi);
    X2 = NSE_trid_c2D(aa, ab, ac, v);
    
    % Compute X*
    Xs = (X1(:,1)+X1(:,end))./(1+X2(:,1)+X2(:,end));
    
    % Final Solution
    fi = X1-Xs(:).*X2;
    

end

function fi = NSE_trid_c2D(aa, ab, ac, fi)

    [m,n] = size(ab);
    X = zeros(m,n);
    
    % Tridiagonal solver (no periodicity)
    
    %Define beta,gamma 
    beta  = zeros(m,n);
    gamma = zeros(m,n);
    
    beta(:,1)  = ab(:,1); 
    gamma(:,1) = fi(:,1)./beta(:,1);
    
    for k = 2:n
        beta(:,k)  = ab(:,k) - (ac(:,k-1)./beta(:,k-1)).*aa(:,k);
        gamma(:,k) = (fi(:,k) - aa(:,k).*gamma(:,k-1))./beta(:,k);
    end
    
    % Solve for X
    X(:,end)     = gamma(:,end);
    for k = n-1:-1:1
        X(:,k) = gamma(:,k) - (ac(:,k)./beta(:,k)).*X(:,k+1);
    end
    
    fi = X;

end