% Computes the source/sink term f

function f = sourcef(X,Y, type)

    param;
    if type == 'L'
        f = (a.^2 + b.^2)*exactsol(X,Y);
    elseif type == 'NL'
        f = (a.^2 + b.^2 - 2*a*sin(a*X).*sin(b*Y)).*exactsol(X,Y);
    end
end