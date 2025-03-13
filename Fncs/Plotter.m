%% figure settings
function Plotter(titlename,xname,yname,xlim,ylim,legendArrary,flag)

xlabel(xname);
ylabel(yname);
title(titlename);
grid on;
box on;
% axis square;
set(gca, 'FontName', 'Arial', 'FontSize', 12);
set(gcf,'Position', [0, 0, 480, 400]);
set(gca, 'LineWidth', 1.25);
set(gca, 'XLim',xlim,'YLim',ylim);
% ylim tight;
if flag.LegendON_OFF
    lgd = legend(legendArrary{:},'Location','best');
    set(lgd,'FontName','Arial','FontSize',12);
    if flag.Legendflage
        set(lgd, 'Color', 'none'); % 设置图例框的颜色为'none'
        set(lgd, 'Box', 'off');
    end
end
end
% lgd = legend('Experimental','Simulated','Experimental R_p','Location','best');
% set(lgd, 'Color', 'none'); % 设置图例框的颜色为'none'
% set(lgd, 'Box', 'off');
% yticks([-5, -4, -3, -2,-1, 0, 1,2,3,4,5]);
