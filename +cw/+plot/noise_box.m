function noise_box(is_noise)
    xl = xlim;
    yl = ylim;
        
    if is_noise
        line(xlim,[yl(1) yl(1)],'LineWidth',10,'Color','r')
        line(xlim,[yl(2) yl(2)],'LineWidth',10,'Color','r')
        line([xl(1) xl(1)], ylim,'LineWidth',10,'Color','r')
        line([xl(2) xl(2)], ylim,'LineWidth',10,'Color','r')
    else
        line(xlim,[yl(1) yl(1)],'LineWidth',10,'Color','g')
        line(xlim,[yl(2) yl(2)],'LineWidth',10,'Color','g')
        line([xl(1) xl(1)], ylim,'LineWidth',10,'Color','g')
        line([xl(2) xl(2)], ylim,'LineWidth',10,'Color','g')
    end
end