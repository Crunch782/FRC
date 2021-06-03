% Creates isocontour plot of the exact and numerical solutions
% Also plots the convergence

function p = plotVisuals(un, ue, X,Y, e, nn, ncont)

    param;
    %figure;
    hold off
    umin = min(min(un));
    umax = max(max(un));
    ucont=umin+[0:ncont-1]*(umax-umin)/(ncont-1);
    
    % Plot the isocontours on the grid (xc,ym)
    
    %subplot(1,2,1)
    contour(X,Y, un, ucont, 'b')
    hold on
    contour(X,Y, ue, ucont, '--r')
    xlabel('xc')
    ylabel('ym')
    title('Isocontours of u_{num} and u_{ex}')
    xlim([xc(1) xc(end)])
    ylim([ym(1) ym(end)])
    legend('u_{num}','u_{ex}')
    
    %{
    subplot(1,2,2)
    semilogy(nn,e)
    title('Convergence')
    xlabel('Iterations')
    ylabel('||u_{n+1} - u_{n}||_2')
    xlim([1 500])
    ylim([0 1])
    %}
    drawnow;
    zoom on;
    
    

end

