dataPth     = '/Volumes/server/Projects/BAIR/Analyses/visual/sub-chaam/analyzePRF_nyu/';

load(fullfile(dataPth, 'chaam_analyzePRF_results_ecog.mat'));

nChan = length(chanNames);

%% Visualize the location of each voxel's pRF

% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 16.6 degrees of visual angle:
cfactor = 16.6/100;

%% Plot each electrode in separate subplots FULL TIMESERIES FITS

figure; hold on;

for cc = 1:nChan
        
    subplot(3,2,cc); hold on
    set(gcf,'Units','points','Position',[100 100 400 400]);

    xpos = results.ecc(cc) * cos(results.ang(cc)/180*pi) * cfactor;
    ypos = results.ecc(cc) * sin(results.ang(cc)/180*pi) * cfactor;
    ang = results.ang(cc)/180*pi;
    sd = results.rfsize(cc) * cfactor;
    h = k_drawellipse(xpos,ypos,ang,sd,sd);  % 
    set(h,'Color','k','LineWidth',2);
    set(scatter(xpos,ypos,'r.'),'CData',[0 0 0]);
    rsq = results.R2(cc);

    % plot edits
    h = k_drawellipse(0,0,0,8.3,8.3); % circle indicating stimulus extent
    set(h,'Color',[0.5 0.5 0.5],'LineWidth',1, 'LineStyle', ':');
    axis([-20 20 -20 20]);
    straightline(0,'h','k:');       % line indicating horizontal meridian
    straightline(0,'v','k:');       % line indicating vertical meridian
    axis square;
    set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
    xlabel('X-position (deg)');
    ylabel('Y-position (deg)');

    title(sprintf('%s R2 = %s', chanNames{cc}, num2str(round(rsq,2))));
end
set(gcf, 'Position', [1000 500 1000 1000])

%% All in one plot FULL

figure; hold on;
set(gcf,'Units','points','Position',[100 100 400 400]);
cmap = hsv(nChan);
for cc = 1:nChan 
      xpos = results.ecc(cc) * cos(results.ang(cc)/180*pi) * cfactor;
      ypos = results.ecc(cc) * sin(results.ang(cc)/180*pi) * cfactor;
      ang = results.ang(cc)/180*pi;
      sd = results.rfsize(cc) * cfactor;
      h = k_drawellipse(xpos,ypos,ang,sd,sd);  % 
      h.Annotation.LegendInformation.IconDisplayStyle = 'off';
      set(h,'Color',cmap(cc,:),'LineWidth',2);
      set(scatter(xpos,ypos,'ro', 'filled'),'CData',cmap(cc,:));
end
h = k_drawellipse(0,0,0,8.3,8.3); % circle indicating stimulus extent
set(h,'Color',[0.5 0.5 0.5],'LineWidth',1, 'LineStyle', ':');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
axis([-20 20 -20 20]);
straightline(0,'h','k-');       % line indicating horizontal meridian
straightline(0,'v','k-');       % line indicating vertical meridian
axis square;
set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
xlabel('X-position (deg)');
ylabel('Y-position (deg)');
legend(chanNames);
set(gca, 'FontSize', 18);
set(gcf, 'Position', [1000 500 1000 1000])

%% Visualize the location of each voxel's pRF 100 BOOTSTRAPS OF TIMESERIES

clear xpos ypos ang sd rsq

% Plot each electrode in separate subplots:
figure; hold on;

for cc = 1:nChan
        
    subplot(3,2,cc); hold on
    set(gcf,'Units','points','Position',[100 100 400 400]);

    for ii = 1:size(results_boot,2)

        theseresults = results_boot(ii); 
        
        xpos(cc,ii) = theseresults.ecc(cc) * cos(theseresults.ang(cc)/180*pi) * cfactor;
        ypos(cc,ii)= theseresults.ecc(cc) * sin(theseresults.ang(cc)/180*pi) * cfactor;
        ang(cc,ii) = theseresults.ang(cc)/180*pi;
        sd(cc,ii) = theseresults.rfsize(cc) * cfactor;
        h = k_drawellipse(xpos(cc,ii),ypos(cc,ii),ang(cc,ii),sd(cc,ii),sd(cc,ii));  
        set(h,'Color',[0.5 0.5 0.5],'LineStyle', '-','LineWidth',1);
        scatter(xpos(cc,ii),ypos(cc,ii),10,[0.5 0.5 0.5], 'o', 'filled');
        rsq(cc,ii) = theseresults.R2(cc);
    end
    
    title(sprintf('%s median R2 = %s', chanNames{cc}, num2str(round(median(rsq(cc,:),2),1))));
    h = k_drawellipse(median(xpos(cc,:),2),median(ypos(cc,:),2),median(ang(cc,:),2),median(sd(cc,:),2),median(sd(cc,:),2)); 
    set(h,'Color','k','LineStyle', '-','LineWidth',3);
    h = scatter(median(xpos(cc,:),2),median(ypos(cc,:),2),'ko', 'filled');
    set(h, 'MarkerEdgeColor', 'k', 'SizeData',10);
    
    % plot edits
    h = k_drawellipse(0,0,0,8.3,8.3); % circle indicating stimulus extent
    set(h,'Color',[0.5 0.5 0.5],'LineWidth',1, 'LineStyle', ':');
    axis([-20 20 -20 20]);
    straightline(0,'h','k:');       % line indicating horizontal meridian
    straightline(0,'v','k:');       % line indicating vertical meridian
    axis square;
    set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
    xlabel('X-position (deg)');
    ylabel('Y-position (deg)');

end
set(gcf, 'Position', [1000 500 1000 1000])

%% All in one plot BOOTSTRAPS

figure; hold on;
set(gcf,'Units','points','Position',[100 100 400 400]);
cmap = hsv(nChan);
for cc = 1:nChan 
      h = k_drawellipse(median(xpos(cc,:),2),median(ypos(cc,:),2),median(ang(cc,:),2),median(sd(cc,:),2),median(sd(cc,:),2));
      h.Annotation.LegendInformation.IconDisplayStyle = 'off';
      set(h,'Color',cmap(cc,:),'LineWidth',2);
      set(scatter(median(xpos(cc,:),2),median(ypos(cc,:),2),'ro', 'filled'),'CData',cmap(cc,:));
      set(h,'Color',cmap(cc,:));
end
h = k_drawellipse(0,0,0,8.3,8.3); % circle indicating stimulus extent
set(h,'Color',[0.5 0.5 0.5],'LineWidth',1, 'LineStyle', ':');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
axis([-20 20 -20 20]);
straightline(0,'h','k-');       % line indicating horizontal meridian
straightline(0,'v','k-');       % line indicating vertical meridian
axis square;
set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
xlabel('X-position (deg)');
ylabel('Y-position (deg)');
legend(chanNames);
set(gca, 'FontSize', 18);
set(gcf, 'Position', [1000 500 1000 1000])
