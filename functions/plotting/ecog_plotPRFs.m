function ecog_plotPRFs(results, stimulus, channels, chan_ind, degPerPix, colorOpt)  

if ~exist('chan_ind', 'var') || isempty(chan_ind)
    chan_ind = 1:height(channels);
end

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

f_ind = checkForHDgrid(channels);


for f = 1:length(f_ind)
    
    figure; hold on
    chan_ind = f_ind{f};
    nChans = length(chan_ind);

    plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
    %tiledlayout(plotDim1, plotDim2, 'TileSpacing','compact'); % Matlab 2019B
    for ii = 1:nChans

        %nexttile
        el = chan_ind(ii);

        subplot(plotDim1,plotDim2,ii); 
        %plotTitle = sprintf('%s R2 %0.1f ecc %0.1f ang %d', channels.name{el}, results.xval(el,:), results.ecc(el)*degPerPix, round(results.ang(el)));        
        %plotTitle = sprintf('%s R2=%0.1f e=%0.1f a=%d', channels.name{el}, results.R2(el,:), results.ecc(el)*degPerPix, round(results.ang(el)));        
        plotTitle = sprintf('%s R2 %0.1f', channels.name{el}, results.R2(el,:));        
        sd = results.rfsize(el);
        %[xx, yy] = meshgrid(linspace(-1,1,res));
        %[th, r] = cart2pol(xx, yy);
        p = results.params(1,:,el);
        im = makegaussian2d(drawRes,p(1)+stimRes,p(2)+stimRes,p(3)/sqrt(p(5)),p(3)/sqrt(p(5))); 

        % plot pRF
        imagesc(im);hold on
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
        h3 = scatter(p(2)+stimRes,p(1)+stimRes,'ko','filled');
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
        setsubplotaxes();
    end
end

end