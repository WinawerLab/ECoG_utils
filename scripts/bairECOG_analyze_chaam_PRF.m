
% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGCAR';

% Dataset specs
subject     = 'chaam'; 
session     = 'umcuecog';
task        = 'prf';
description = 'reref';

[data, channels, events] = bidsEcogGetPreprocData(dataPth, subject, session, task, [], description);