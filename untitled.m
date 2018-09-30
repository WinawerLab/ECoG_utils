fileName = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/648/NY648_Winawer.edf';

data = ft_read_data(fileName);
hdr = ft_read_header(fileName);

ft_write_data('testtest.edf', data, 'header', hdr, 'dataformat', 'edf');

data2 = ft_read_data('testtest.edf');


data = edf2fieldtrip(fileName);

cfg            = [];
cfg.dataset    = fileName;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

%%
%%
%fileName = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/648/NY648_Winawer.edf';

cd ~/Desktop
fileName = 'examplerecording.edf';
data_in = ft_read_data(fileName); 
hdr_in = ft_read_header(fileName);

ft_write_data('exampledata.edf', data_in, 'header', hdr_in, 'dataformat', 'edf');
data_out = ft_read_data('exampledata.edf');

ft_write_data('exampledata', data_in, 'header', hdr_in, 'dataformat', 'brainvision_eeg');
data_out = ft_read_data('exampledata.eeg');



