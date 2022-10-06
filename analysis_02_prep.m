% Preprocess EEG for the continuous feedback processing project
% 
% Other m-files required: 
% EEGLAB toolbox https://github.com/sccn/eeglab
% /private/ccn_prep.m
% /private/ccnlabactichamp.locs
% /private/ccn_check.m
% /private/find_artifacts.m
% /private/make_erp.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Participant strings
ps = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};

% Set data folders - set as needed
if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_Gnomes_Hassall\data';
else
    dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_Gnomes_Hassall/data';
end

% ICA settings - which triggers? Window size?
icaTriggers = {'S  1','S  2','S  3','S  4','S  5','S  6'}; % Fixation cross
icaWindow = [-0.2  2.8];

% Artifact check settings - which triggers? Window size? Baseline?
checkTriggers = icaTriggers;
checkWindow = [-0.2 0.6];
baseline = [-200 0];

toRemove = {};
allBadChannels = {};

filters = [0.1 30];
badChannels = {};
reference = {'TP9','TP10'};
appendString = '';

artifactSettings.maxMin = 100;
artifactSettings.level = 100;
artifactSettings.step = 40;
artifactSettings.lowest = 0.1;

for p = 1:length(ps)

    rng(2021); % Set for consistency
    
    subName = ['sub-' ps{p}];
    taskName = 'gnomes';
    
    ccn_prep(dataFolder,subName,taskName,icaTriggers,icaWindow,filters,reference,badChannels,appendString);
    allBadChannels{p} = ccn_check(dataFolder,subName,taskName,icaTriggers,checkWindow,baseline,artifactSettings);
end