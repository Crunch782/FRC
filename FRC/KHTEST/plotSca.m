% Plots a pcolor plot of the scalar field

function h = plotSca(C,x,y)

    %Grid
    [X,Y] = meshgrid(x,y);
    X = X';
    Y = Y';
    
    h = pcolor(X,Y,C);
    set(h, 'EdgeColor', 'none');
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
    xlabel('xm')
    ylabel('ym')
    title('Scalar Field')
    colorbar
    
end