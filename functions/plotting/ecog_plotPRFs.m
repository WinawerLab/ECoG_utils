function ecog_plotPRFs(results, stimulus, channels, degPerPix, colorOpt)  

if ~exist('degPerPix', 'var') || isempty(degPerPix)
    % The stimulus is 100 pixels (in both height and weight), and this corresponds to
    % 16.6 degrees of visual angle:
    degPerPix= 16.6/100;
end

if ~exist('colorOpt', 'var') || isempty(colorOpt)
    colorOpt = 0; % 0 = gray, 1 = parula 
end
    
if ~iscell(stimulus), stimulus = {stimulus}; end

% PRFs
stimRes = size(stimulus{1},1);
% draw the prfs at 3x the screen size (adjust axes later)
drawRes = stimRes * 3;
% define center pix 
centerPix = drawRes/2;

figure; hold on
nChans = size(results.ang,1);

plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
%tiledlayout(plotDim1, plotDim2, 'TileSpacing','compact'); % Matlab 2019B
for el = 1:nChans
    %nexttile
    subplot(plotDim1,plotDim2,el); hold on
    %plotTitle = sprintf('%s %s %s R2: %0.1f ecc: %0.1f ang: %0.1f', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el}, ...
    %    results.R2(el), results.ecc(el)*degPerPix, results.ang(el));        
    plotTitle = sprintf('%s R2 %0.1f ecc %0.1f ang %d', channels.name{el}, results.R2(el), results.ecc(el)*degPerPix, round(results.ang(el)));        
    
    sd = results.rfsize(el);
    %[xx, yy] = meshgrid(linspace(-1,1,res));
    %[th, r] = cart2pol(xx, yy);
    p = results.params(1,:,el);
    im = makegaussian2d(drawRes,p(1)+stimRes,p(2)+stimRes,p(3)/sqrt(p(5)),p(3)/sqrt(p(5))); 

    % plot pRF
    imagesc(im);
    if colorOpt == 1
        colormap(parula)
    else
        colormap(1-gray)
    end

    % draw stimulus
    h1 = k_drawellipse(centerPix,centerPix,0,drawRes/6,drawRes/6); % circle indicating stimulus extent
    set(h1,'Color',[0 0 0],'LineWidth',1, 'LineStyle', ':');
    h2 = straightline(centerPix,'h','k:');     % line indicating horizontal meridian
    h3 = straightline(centerPix,'v','k:');     % line indicating vertical meridian
    if colorOpt == 1
        set(h1,'Color',[1 1 1]);
        set(h2,'Color',[1 1 1]);
        set(h3,'Color',[1 1 1]);
    end
    % plot pRF center and sd  
    h1 = k_drawellipse(p(2)+stimRes,p(1)+stimRes,0,sd,sd);      % 
    h2 = k_drawellipse(p(2)+stimRes,p(1)+stimRes,0,2*sd,2*sd);  % 
    set(h1,'Color', [0 0 0],'LineWidth',2,'LineStyle', '-');
    set(h2,'Color', [0 0 0],'LineWidth',2,'LineStyle', '-');
    h3 = scatter(p(2)+stimRes/2,p(1)+stimRes/2,'wo','filled');
    if colorOpt == 1
        set(h3,'CData',[1 0 0]);
    end
    
    % format axes
    title(plotTitle);
    axis square;
    set(gca, 'XTick', [0:drawRes/6:drawRes], 'XTickLabel', round(([0:drawRes/6:drawRes]*degPerPix)-centerPix*degPerPix));
    set(gca, 'YTick', [0:drawRes/6:drawRes], 'YTickLabel', round(([0:drawRes/6:drawRes]*degPerPix)-centerPix*degPerPix));
    %if el == 1; xlabel('X-position (deg)'); ylabel('Y-position (deg)');end
    
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
    set(gca, 'XLim', [stimRes/2 drawRes-stimRes/2],'YLim',[stimRes/2 drawRes-stimRes/2])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    set(gca, 'FontSize', 14)
end

end