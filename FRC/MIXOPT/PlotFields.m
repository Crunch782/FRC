% Plots both the Vorticity and Scalar Fields

function pf = PlotFields(c,xm,ym,w,xc,yc,tag,t,vl,vu,nff,s,T,Re)

    pf = 1;
    figure;
    
    subplot(1,2,1)
    plotVort(w,xc,yc,vl,vu);
    
    subplot(1,2,2)
    plotSca(c,xm,ym);
    tl = append('t = ',num2str(t));
    nm = num2str(nff);
    sgtitle(tl)
    ttag = append('\bf (s, T, Re) = (',num2str(s),', ',num2str(T),', ',num2str(Re),')');
    text(-6, -4, ttag)
    
    fpath = char(strcat('Data/',tag,'/Figures'));
    
    saveas(gcf, fullfile(fpath, nm), 'jpg'); 

end