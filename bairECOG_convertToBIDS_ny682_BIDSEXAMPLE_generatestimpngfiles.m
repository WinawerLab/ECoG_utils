
%% VISUALIZE STIMULI
for ii = 1:size(stimulus.images,3)-1
    subplot(6,6,ii);hold on
    imshow(stimulus.images(:,:,ii));
end

%% WRITE OUT STIMULI
readDir = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/682/stimdata/';
writeDir = '/Users/winawerlab/matlab/git/bids-examples/ieeg_visual_multimodal/stimuli';

taskList = {'spatialobject', 'spatialpattern', 'temporalpattern'};

%% RUN 1

for kk = 1:length(taskList)
    load(fullfile(readDir, sprintf('sub-ny682_ses-nyuecog01_task-%s_run-1.mat', taskList{kk})));

    stimcount = 0;
    for ii = 1:size(stimulus.im_cell,2)
        im = stimulus.im_cell{ii};
        for jj = 1:size(im,3)
            %I = im(:,:,jj); % WRITE EMPTY STIMULI
            I = zeros(size(im(:,:,jj)));
            stimcount = stimcount + 1;
            if jj < 10
                filename = sprintf('stim_%d_%s_0%d', stimulus.cat(stimcount),stimulus.categories{ii}, jj);
            else
                filename = sprintf('stim_%d_%s_%d', stimulus.cat(stimcount),stimulus.categories{ii}, jj);
            end
            filename = [filename '.png'];    
            imwrite(I,fullfile(writeDir,filename),'png');
        end
    end
end

%% RUN 2

for kk = 1:length(taskList)
    load(fullfile(readDir, sprintf('sub-ny682_ses-nyuecog01_task-%s_run-2.mat', taskList{kk})));
    stimcount = 0;
    for ii = 1:size(stimulus.im_cell,2)
        im = stimulus.im_cell{ii};
        ADD = size(im,3);
        for jj = 1:size(im,3)
            %I = im(:,:,jj); % WRITE EMPTY STIMULI
            I = zeros(size(im(:,:,jj)));
            stimcount = stimcount + 1;
            if jj + ADD < 10
                filename = sprintf('stim_%d_%s_0%d', stimulus.cat(stimcount),stimulus.categories{ii}, jj+ADD);
            else
                filename = sprintf('stim_%d_%s_%d', stimulus.cat(stimcount),stimulus.categories{ii}, jj+ADD);
            end
            filename = [filename '.png'];    
            imwrite(I,fullfile(writeDir,filename),'png');
        end
    end
end
