function ecog_plotPRFtsfits(data,results,tr, opt,channels)


% Define some variables
res = [100 100];                    % row x column resolution of the stimuli
resmx = 100;                        % maximum resolution (along any dimension)
hrf = results.options.hrf;          % HRF that was used in the model
degs = results.options.maxpolydeg;  % vector of maximum polynomial degrees used in the model

% Pre-compute cache for faster execution
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);

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

for el = 1:nChans
    subplot(plotDim1,plotDim2,el); hold on
    plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        

    vx = el;
    % For each run, collect the data and the model fit.  We project out polynomials
    % from both the data and the model fit.  This deals with the problem of
    % slow trends in the data.
    datats = {};
    modelts = {};
    for cc=1:length(data)
      datats{cc} =  polymatrix{cc}*data{cc}(vx,:)';
      modelts{cc} = polymatrix{cc}*modelfun(results.params(1,:,vx),stimulusPP{cc});
    end

    set(gcf,'Units','points');
    plot(cat(1,datats{:}),'k-', 'LineWidth', 2);
    plot(cat(1,modelts{:}),'r-','LineWidth', 2);
    xlabel('PRF stimulus','FontSize', 28);
    if el == 1, legend('Data', 'Model prediction'); end
    set(gca, 'FontSize', 18)
    title(plotTitle);
    axis tight
end

end