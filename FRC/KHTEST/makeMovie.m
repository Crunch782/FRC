% Creates a movie of the flow evolution

function m = makeMovie(W,x,y,vid,t)

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
    title(num2str(t))
    colorbar
    drawnow
    image = getframe(gcf);
    size(image.cdata)
    writeVideo(vid,image);

end