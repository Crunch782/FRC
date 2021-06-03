% Plots a pcolor plot of the vorticity field

function h = plotVort(W,x,y,vl,vu)

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
    %caxis([vl vu])
    pbaspect([1 1 1])
    
end