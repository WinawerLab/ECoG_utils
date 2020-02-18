function ecog_plotPRFs(channels, results, stimulus, coloropt)  

% PRFs
res_in_pix = size(stimulus{1},1);
% define pix center

for el = 1:nChans
    
    subplot(plotDim1,plotDim2,el); hold on
    plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        

    sd = results.rfsize(el);
    %[xx, yy] = meshgrid(linspace(-1,1,res));
    %[th, r] = cart2pol(xx, yy);
    p = results.params(1,:,el);
    im = makegaussian2d(250,p(1)+75,p(2)+75,p(3)/sqrt(p(5)),p(3)/sqrt(p(5))); 

    % plot pRF
    imagesc(im);
    if coloropt == 1
        colormap(parula)
    else
        colormap(1-gray)
    end

    % draw stimulus
    h1 = k_drawellipse(125,125,0,50,50); % circle indicating stimulus extent
    set(h1,'Color',[0 0 0],'LineWidth',1, 'LineStyle', ':');
    h2 = straightline(125,'h','k:');     % line indicating horizontal meridian
    h3 = straightline(125,'v','k:');     % line indicating vertical meridian
    if coloropt == 1
        set(h1,'Color',[1 1 1]);
        set(h2,'Color',[1 1 1]);
        set(h3,'Color',[1 1 1]);
    end
    % plot pRF center and sd  
    h1 = k_drawellipse(p(2)+75,p(1)+75,0,sd,sd);      % 
    h2 = k_drawellipse(p(2)+75,p(1)+75,0,2*sd,2*sd);  % 
    set(h1,'Color', [0 0 0],'LineWidth',2,'LineStyle', '-');
    set(h2,'Color', [0 0 0],'LineWidth',2,'LineStyle', '-');
    h3 = scatter(p(2)+75,p(1)+75,'wo','filled');
    if coloropt == 1
        set(h3,'CData',[1 0 0]);
    end

    axis square;
    set(gca, 'XTick', [1 25:25:250], 'XTickLabel', round(([1 25:25:250]*0.16)-20));
    set(gca, 'YTick', [1 25:25:250], 'YTickLabel', round(([1 25:25:250]*0.16)-20));
    xlabel('X-position (deg)');
    ylabel('Y-position (deg)');
    title(plotTitle);
    set(gca, 'XLim', [25 225],'YLim',[25 225])
end
saveas(gcf, fullfile(plotSaveDir,'modelfits', figureName), 'png'); close;

end