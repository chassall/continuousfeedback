% This script adds files and information to a BIDS folder that has been 
% created by Brain Visions's BV2BIDS command line tool
% 
% Other m-files required: 
% bids-matlab toolbox (https://github.com/bids-standard/bids-matlab)

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% specify BIDS folder
bidsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_Gnomes_Hassall/data';

% specify source folder (contain the raw behavioural files
% that will be bidsified below)
sourceFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_Gnomes_Hassall/data_source';

%% read dataset_description.json, created by BV2BIDS, and modify
% https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
datasetFile = fullfile(bidsFolder, 'dataset_description.json');
dataset = bids.util.jsondecode(datasetFile);
dataset.DatasetType = 'raw';
% dataset.ReferencesAndLinks = {''};
options.indent = '  '; % Adds white space, easier to read
bids.util.jsonwrite(datasetFile, dataset, options);

%% write participant tsv files, not created by BV2BIDS
% https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html

% load participant file
participantFile = fullfile(bidsFolder, 'participants.tsv');
participants.participant_id = {'sub-01'; 'sub-02'; 'sub-03'; 'sub-04'; 'sub-05'; 'sub-06'; 'sub-07'; 'sub-08'; 'sub-09'; 'sub-10'; 'sub-11'; 'sub-12'; 'sub-13'; 'sub-14'; 'sub-15'; 'sub-16'; 'sub-17'; 'sub-18'; 'sub-19'; 'sub-20'; 'sub-21'};
participants.participant = {'01'; '02';'03'; '04';'05';'06';'07';'08';'09'; '10'; '11'; '12'; '13'; '14'; '15'; '16'; '17'; '18'; '19'; '20'; '21'};
participants.datetime = {'19-May-2021 11:56:37'; '21-May-2021 10:24:48'; '26-May-2021 09:36:23';... 
'27-May-2021 10:34:13'; '28-May-2021 11:00:31'; '03-Jun-2021 11:12:24';...
'10-Jun-2021 10:43:34'; '14-Jun-2021 12:27:05'; '16-Jun-2021 09:56:57';...
'16-Jun-2021 12:57:28'; '10-Nov-2021 15:30:09'; '17-Nov-2021 14:58:10'; ...
'19-Nov-2021 14:31:18'; '26-Nov-2021 14:34:36'; '03-Dec-2021 13:41:07';...
'08-Dec-2021 14:35:15'; '06-Jan-2022 15:11:54'; '02-Feb-2022 10:23:30';...
'10-Feb-2022 11:20:46'; '16-Feb-2022 10:44:16'; '25-Feb-2022 13:03:09'};

participants.age = [23; 22; 23; 25; 27; 29; 26; 41; 33; 25; 24; 25; 24; 22; 28; 27; 21; 23; 24; 26; 24];
participants.sex = {'F';'M'; 'M'; 'F'; 'F'; 'F'; 'F'; 'M'; 'F'; 'M'; 'F'; 'F'; 'F'; 'F'; 'F'; 'F'; 'F'; 'F'; 'F'; 'M'; 'F'};

participants.handedness = {'R'; 'L'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'L'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'};

bids.util.tsvwrite(participantFile, participants);

%% write participant json file, not created by BV2BIDS
% https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
participantJSONFile = fullfile(bidsFolder, 'participants.json');
pInfoDesc.participant.Description = 'participant number as input by the tester';
pInfoDesc.datetime.Description = 'date and time at start of task';
pInfoDesc.age.Description = 'self-reported age of participant';
pInfoDesc.age.Units = 'years';
pInfoDesc.sex.Description = 'self-reported sex of participant';
pInfoDesc.sex.Levels.M = 'male';
pInfoDesc.sex.Levels.F = 'female';
pInfoDesc.handedness.Description = 'self-reported handedness of participant';
pInfoDesc.handedness.Levels.L = 'left-handed';
pInfoDesc.handedness.Levels.R = 'right-handed';
pInfoDesc.handedness.Levels.LR = 'ambidextrous';
options.indent = '  '; % Adds white space, easier to read
bids.util.jsonwrite(participantJSONFile,pInfoDesc,options);

%% read eeg json for each participant, created by BV2BIDS, and modify
% https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html#electroencephalography
whichPs = 1:21;
for p = 1:length(whichPs)
    pString = ['sub-' num2str(whichPs(p),'%0.02i')];
    eegJSONFile = fullfile(bidsFolder, pString,'eeg',[pString '_task-gnomes_eeg.json']);
    eeg = bids.util.jsondecode(eegJSONFile);
    eeg.InstitutionName = 'University of Oxford';
    eeg.InstitutionAddress = 'Warneford Hospital, Oxford, OX3 7JX';
    eeg.InstitutionalDepartmentName = 'Department of Psychiatry';
    eeg.ManufacturersModelName = 'actiCHamp Plus';
    options.indent = '  '; % Adds white space, easier to read
    bids.util.jsonwrite(eegJSONFile, eeg, options);
end

%% write events json for each participant, not created by BV2BIDS
whichPs = 1:21;
for p = 1:length(whichPs)
    pString = ['sub-' num2str(whichPs(p),'%0.02i')];
    eventsJSONFile = fullfile(bidsFolder, pString,'eeg',[pString '_task-gnomes_events.json']);

    eInfoDesc.trial_type.Description = 'BrainVision event number and type (string)';
    eInfoDesc.value.Description = 'Event value (string)';
    eInfoDesc.value.Levels.S1 = 'Condition 1 fixation cross';
    eInfoDesc.value.Levels.S2 = 'Condition 2 fixation cross';
    eInfoDesc.value.Levels.S3 = 'Condition 3 fixation cross';
    eInfoDesc.value.Levels.S4 = 'Condition 4 fixation cross';
    eInfoDesc.value.Levels.S5 = 'Condition 5 fixation cross';
    eInfoDesc.value.Levels.S6 = 'Condition 6 fixation cross';
    eInfoDesc.value.Levels.S11 = 'Condition 1 cue';
    eInfoDesc.value.Levels.S12 = 'Condition 2 cue';
    eInfoDesc.value.Levels.S13 = 'Condition 3 cue';
    eInfoDesc.value.Levels.S14 = 'Condition 4 cue';
    eInfoDesc.value.Levels.S15 = 'Condition 5 cue';
    eInfoDesc.value.Levels.S16 = 'Condition 6 cue';
    eInfoDesc.value.Levels.S21 = 'Condition 1 bar';
    eInfoDesc.value.Levels.S22 = 'Condition 2 bar';
    eInfoDesc.value.Levels.S23 = 'Condition 3 bar';
    eInfoDesc.value.Levels.S24 = 'Condition 4 bar';
    eInfoDesc.value.Levels.S25 = 'Condition 5 bar';
    eInfoDesc.value.Levels.S26 = 'Condition 6 bar';
    eInfoDesc.value.Levels.S31 = 'Condition 1 response';
    eInfoDesc.value.Levels.S32 = 'Condition 2 response';
    eInfoDesc.value.Levels.S33 = 'Condition 3 response';
    eInfoDesc.value.Levels.S34 = 'Condition 4 response';
    eInfoDesc.value.Levels.S35 = 'Condition 5 response';
    eInfoDesc.value.Levels.S36 = 'Condition 6 response';
    eInfoDesc.value.Levels.S41 = 'Condition 1 animation start';
    eInfoDesc.value.Levels.S42 = 'Condition 2 animation start';
    eInfoDesc.value.Levels.S43 = 'Condition 3 animation start';
    eInfoDesc.value.Levels.S44 = 'Condition 4 animation start';
    eInfoDesc.value.Levels.S45 = 'Condition 5 animation start';
    eInfoDesc.value.Levels.S46 = 'Condition 6 animation start';
    eInfoDesc.value.Levels.S51 = 'Condition 1 animation end';
    eInfoDesc.value.Levels.S52 = 'Condition 2 animation end';
    eInfoDesc.value.Levels.S53 = 'Condition 3 animation end';
    eInfoDesc.value.Levels.S54 = 'Condition 4 animation end';
    eInfoDesc.value.Levels.S55 = 'Condition 5 animation end';
    eInfoDesc.value.Levels.S56 = 'Condition 6 animation end';
    eInfoDesc.onset.Description = 'Event onset';
    eInfoDesc.onset.Units = 'milisecond';
    eInfoDesc.duration.Description = 'Event duration';
    eInfoDesc.duration.Units = 'milisecond';
    eInfoDesc.channel.Description = 'Channel number';
    eInfoDesc.channel.Levels.x0 = 'All channels';
    options.indent = '  '; % Adds white space, easier to read
    bids.util.jsonwrite(eventsJSONFile,eInfoDesc,options);
end

%% load behavioural data and save as tsv files
whichPs = 1:21;
for p = 1:length(whichPs)
    pString = num2str(whichPs(p),'%0.02i');
    pBIDSString = ['sub-' pString];
    thisDir = dir(fullfile(sourceFolder,pBIDSString,'beh','*.txt')); % This file has a header
    thisData = readtable(fullfile(sourceFolder,pBIDSString,'beh',thisDir.name));

    thisData = renamevars(thisData,thisData.Properties.VariableNames,{'trial','condition','image','target','fixationTime','rtStart','rtEnd','rt','responseX','responseY','responseProp','pointsThisRound','pointsTotal','preFeedbackTime'});

    % Make a struct out of the behavioural data
    beh.trial = thisData.trial;
    beh.condition = thisData.condition;
    beh.image = thisData.image;
    beh.target = thisData.target;
    beh.fixationTime = thisData.fixationTime;
    beh.rtStart = thisData.rtStart;
    beh.rtEnd = thisData.rtEnd;
    beh.rt = thisData.rt;
    beh.responseX = thisData.responseX;
    beh.responseY = thisData.responseY;
    beh.responseProp = thisData.responseProp;
    beh.pointsThisRound = thisData.pointsThisRound;
    beh.pointsTotal = thisData.pointsTotal;
    beh.preFeedbackTime = thisData.preFeedbackTime;
    
    behFolder = fullfile(bidsFolder,pBIDSString,'beh');
    if ~exist(behFolder, 'dir')
        mkdir(behFolder);
    end
    behFile = [pBIDSString '_task-gnomes_beh.tsv'];
    bids.util.tsvwrite(fullfile(behFolder,behFile), beh);

end

%% write beh json for each participant
whichPs = 1:21;
for p = 1:length(whichPs)
    pString = ['sub-' num2str(whichPs(p),'%0.02i')];
    behJSONFile = fullfile(bidsFolder, pString,'beh',[pString '_task-gnomes_beh.json']);

    bInfoDesc.trial.Description = 'Trial number (integer)';	
    bInfoDesc.condition.Description	 = 'Condition number (integer) - see README for details';
    bInfoDesc.image.Description	 = 'Image number (integer) - see /code/task/images for images';
    bInfoDesc.target.Description = 'Target level - the animated bar will rise to this height (proportion)';
    bInfoDesc.fixationTime.Description = 'Fixation cross duration (seconds)';
    bInfoDesc.rtStart.Description = 'Start of response period (seconds relative to trial start)';
    bInfoDesc.rtEnd.Description = 'Response time (seconds relative to trial start)';
    bInfoDesc.rt.Description = 'Response time (seconds)';
    bInfoDesc.responseX.Description = 'Horizontal response location (pixels)';
    bInfoDesc.responseY.Description = 'Vertical response location (pixels)';
    bInfoDesc.responseProp.Description = 'Response (proportion of bar height)';
    bInfoDesc.pointsThisRound.Description = 'Points for this trial (integer)';
    bInfoDesc.pointsTotal.Description = 'Total points (integer)';
    bInfoDesc.preFeedbackTime.Description = 'Delay before feedback onset (seconds)';

    options.indent = '  '; % Adds white space, easier to read
    delete(behJSONFile);
    bids.util.jsonwrite(behJSONFile,bInfoDesc,options);
end