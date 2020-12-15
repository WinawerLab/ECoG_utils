function ecog_plotPRFtsfits(data, stimulus, results, channels, chan_ind)

if ~exist('chan_ind', 'var') || isempty(chan_ind)
    chan_ind = 1:height(channels);
end

% Plot PRF timecourses data and fits
% Using example code from Kendrick Kays website
if ~iscell(data), data = {data}; end
if ~iscell(stimulus), stimulus = {stimulus}; end

% Define some variables
nPixels = size(stimulus{1}, 1);
res = [nPixels nPixels];                % row x column resolution of the stimuli
resmx = nPixels;                        % maximum resolution (along any dimension)
hrf = results.options.hrf;              % HRF that was used in the model
degs = results.options.maxpolydeg;      % vector of maximum polynomial degrees used in the model

% Pre-compute cache for faster execution
[~,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% Prepare the stimuli for use in the model
stimulusPP = {};
for cc=1:length(stimulus)
  stimulusPP{cc} = squish(stimulus{cc},2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP{cc} = [stimulusPP{cc} cc*ones(size(stimulusPP{cc},1),1)];  % this adds a dummy column to indicate run breaks
end

% Define the model function.  This function takes parameters and stimuli as input and
% returns a predicted time-series as output.  Specifically, the variable <pp> is a vector
% of parameter values (1 x 5) and the variable <dd> is a matrix with the stimuli (frames x pixels).
% Although it looks complex, what the function does is pretty straightforward: construct a
% 2D Gaussian, crop it to <res>, compute the dot-product between the stimuli and the
% Gaussian, raise the result to an exponent, and then convolve the result with the HRF,
% taking care to not bleed over run boundaries.
modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));

% Construct projection matrices that fit and remove the polynomials.
% Note that a separate projection matrix is constructed for each run.
polymatrix = {};
for cc=1:length(degs)
  polymatrix{cc} = projectionmatrix(constructpolynomialmatrix(size(data{cc},2),0:degs(cc)));
end

f_ind = checkForHDgrid(channels);

for f = 1:length(f_ind)
    
    figure; hold on
    chan_ind = f_ind{f};
    nChans = length(chan_ind);

    plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);


    for ii  = 1:nChans

        el = chan_ind(ii);

        subplot(plotDim1,plotDim2,ii); hold on
        plotTitle = sprintf('%s %s %s R2 %0.1f', channels.name{el},channels.benson14_varea{el}, channels.wang15_mplbl{el},results.R2(el));        
        % For each run, collect the data and the model fit.  We project out polynomials
        % from both the data and the model fit.  This deals with the problem of
        % slow trends in the data.
        datats = {};
        modelts = {};
        for cc=1:length(data)
          datats{cc} =  polymatrix{cc}*data{cc}(el,:)';
          modelts{cc} = polymatrix{cc}*modelfun(results.params(1,:,el),stimulusPP{cc});
        end

        set(gcf,'Units','points');
        plot(cat(1,datats{:}),'k-', 'LineWidth', 2);
        plot(cat(1,modelts{:}),'r-','LineWidth', 2);
        if el == 1, legend('Data', 'Model prediction'); end
        %if el == 1, xlabel('PRF stimulus'); ylabel('Broadband response');end
        set(gca, 'FontSize', 14)
        set(gca, 'XTickLabel', []);
        title(plotTitle);
        axis tight
        setsubplotaxes();
    end
end
end