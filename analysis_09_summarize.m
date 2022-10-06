% Combine results across participants, and save in results folder
% 
% Other m-files required: 
% EEGLAB toolbox https://github.com/sccn/eeglab
% Unfold toolbox: https://github.com/unfoldtoolbox/unfold
% /private/num2bv.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

close all; clear all;

if ispc
    projectFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_Gnomes_Hassall';
else
    projectFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_Gnomes_Hassall';
end
dataFolder = [projectFolder '/data'];
resultsFolder = [projectFolder '/analysis/results'];

% Participant numbers
ps = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};

%% GLM, collapsed across condition, animation start included, no
% regularization, all trials included
% Output from analysis_07_glm_collapsed.m
incAnimStart = 1;
useReg = 0;
excludeEarly = 0; 
allX = {};
allArtifactProp = [];
allBeta = [];
for p = 1:length(ps)
    saveFile = ['sub-' ps{p} '_task-gnomes_glmcollapsed_' num2str(incAnimStart) '_' num2str(useReg) '_' num2str(excludeEarly) '.mat'];
    saveFolder = [dataFolder '/derivatives/glmres/sub-' ps{p}];
    load(fullfile(saveFolder,saveFile),'chanlocs','srate','betaI','X','beta','artifactProp');
    allX{p} = X;
    allArtifactProp(p) = artifactProp;
    allBeta(p,:,:) = beta;
end

save(fullfile(resultsFolder,['deconv_collapsed_' num2str(incAnimStart) '_' num2str(useReg)  '_' num2str(excludeEarly) '.mat'])); 

%% GLM, collapsed across condition, animation start included, WITH
% regularization, all trials included
% Output from analysis_07_glm_collapsed.m
incAnimStart = 1;
useReg = 1;
excludeEarly = 0; 
allX = {};
allArtifactProp = [];
allBeta = [];
allCVErrors = [];
for p = 1:length(ps)
    saveFile = ['sub-' ps{p} '_task-gnomes_glmcollapsed_' num2str(incAnimStart) '_' num2str(useReg) '_' num2str(excludeEarly) '.mat'];
    saveFolder = [dataFolder '/derivatives/glmres/sub-' ps{p}];
    load(fullfile(saveFolder,saveFile),'chanlocs','srate','betaI','X','beta','artifactProp','cvErrors');
    allX{p} = X;
    allArtifactProp(p) = artifactProp;
    allBeta(p,:,:) = beta;
end
save(fullfile(resultsFolder,'cv_results_colapsed.mat'),'allCVErrors');
save(fullfile(resultsFolder,['deconv_collapsed_' num2str(incAnimStart) '_' num2str(useReg)  '_' num2str(excludeEarly) '.mat'])); 


%% GLM, split by condition, animation start included, no regularization
% Output from analysis_07_glm_deconv.m
incAnimStart = 1;
useReg = 0;
allX = {};
allArtifactProp = [];
allBeta = [];
for p = 1:length(ps)
    saveFile = ['sub-' ps{p} '_task-gnomes_glm_' num2str(incAnimStart) '_' num2str(useReg) '.mat'];
    saveFolder = [dataFolder '/derivatives/glmres/sub-' ps{p}];
    load(fullfile(saveFolder,saveFile),'chanlocs','srate','betaI','X','beta','artifactProp','cvErrors');
    allX{p} = X;
    allArtifactProp(p) = artifactProp;
    allBeta(p,:,:) = beta;
end
save(fullfile(resultsFolder,['deconv_' num2str(incAnimStart) '_' num2str(useReg)  '.mat'])); 

%% GLM, split by condition, animation start included, WITH regularization
% Output from analysis_07_glm_deconv.m
incAnimStart = 1;
useReg = 1;
allX = {};
allArtifactProp = [];
allBeta = [];
allCVErrors = [];
for p = 1:length(ps)
    saveFile = ['sub-' ps{p} '_task-gnomes_glm_' num2str(incAnimStart) '_' num2str(useReg) '.mat'];
    saveFolder = [dataFolder '/derivatives/glmres/sub-' ps{p}];
    load(fullfile(saveFolder,saveFile),'chanlocs','srate','betaI','X','beta','artifactProp','cvErrors');
    allX{p} = X;
    allArtifactProp(p) = artifactProp;
    allBeta(p,:,:) = beta;
    allCVErrors(p,:) = cvErrors;
end
save(fullfile(resultsFolder,'cv_results.mat'),'allCVErrors');
save(fullfile(resultsFolder,['deconv_' num2str(incAnimStart) '_' num2str(useReg)  '.mat'])); 


%% GLM, collapsed across conditions, inc. average reward regressor
% Output from analysis_07_glm_collapsedaverew.m
incAnimStart = 1;
useReg = 0;
allX = {};
allArtifactProp = [];
allBeta = [];
for p = 1:length(ps)
    saveFile = ['sub-' ps{p} '_task-gnomes_glmcollapsedaverew_' num2str(incAnimStart) '_' num2str(useReg) '.mat'];
    saveFolder = [dataFolder '/derivatives/glmres/sub-' ps{p}];
    load(fullfile(saveFolder,saveFile),'chanlocs','srate','betaI','X','beta','artifactProp','cvErrors');
    allX{p} = X;
    allArtifactProp(p) = artifactProp;
    allBeta(p,:,:) = beta;
end
save(fullfile(resultsFolder,['deconv_collapsedaverew_' num2str(incAnimStart) '_' num2str(useReg) '.mat'])); 

%% GLM, collapsed across conditions, inc. average reward regressor
% Output from analysis_07_glm_collapsedaverew.m
incAnimStart = 1;
useReg = 1;
allX = {};
allArtifactProp = [];
allBeta = [];
for p = 1:length(ps)
    saveFile = ['sub-' ps{p} '_task-gnomes_glmcollapsedaverew_' num2str(incAnimStart) '_' num2str(useReg) '.mat'];
    saveFolder = [dataFolder '/derivatives/glmres/sub-' ps{p}];
    load(fullfile(saveFolder,saveFile),'chanlocs','srate','betaI','X','beta','artifactProp','cvErrors');
    allX{p} = X;
    allArtifactProp(p) = artifactProp;
    allBeta(p,:,:) = beta;
end
save(fullfile(resultsFolder,['deconv_collapsedaverew_' num2str(incAnimStart) '_' num2str(useReg) '.mat'])); 
