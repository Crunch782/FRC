% Plots a pcolor plot of the vorticity field

function h = plotVort(W,x,y)

    %Grid
    [X,Y] = meshgrid(x,y);
    X = X';
    Y = Y';
    
    h = pcolor(X,Y,W);
    set(h, 'EdgeColor', 'none');
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
    xlabel('xc')
    ylabel('yc')
    title('Vorticity Field')
    colorbar
    
end