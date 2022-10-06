% Correlational analysis for the continuous feedback processing project
% 
% Other m-files required: 
% EEGLAB toolbox https://github.com/sccn/eeglab
% Unfold toolbox: https://github.com/unfoldtoolbox/unfold

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

close all; clear all; init_unfold();

incBarHeight = 0; % Include regressor for bar height?
split = 0; % Split rising/falling reward regressors
excludeFall = 0; % Exclude falling part of reward

if ispc
    projectFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_Gnomes_Hassall';
else
    projectFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_Gnomes_Hassall';
end
dataFolder = [projectFolder '/data'];
resultsFolder = [projectFolder '/analysis/results'];

% Participant numbers
ps = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};

allCoeff = [];
allCoeffPE = [];
allCoeffCollapsed = [];
allPower = [];

whichAnalysis = 'eeg';
allArtifacts = [];

% Loop through participants
for p = 1:length(ps)
    
    % Load EEG
    load(fullfile([dataFolder '/derivatives/eegprep/sub-' ps{p}],['/sub-' ps{p} '_task-gnomes_eegprep.mat']),'EEG');
    bEEG = EEG;
    
    switch whichAnalysis
        case 'eeg'
            whichBand = [];
        case 'delta'
           whichBand = delta;
           bEEG.data = whichBand;
        case 'theta'
            whichBand = theta;
            bEEG.data = whichBand;
        case 'alpha'
            whichBand = alpha;
            bEEG.data = whichBand;
        case 'beta'
            whichBand = beta;
            bEEG.data = whichBand;
        case 'gamma'
            whichBand = gamma;
            bEEG.data = whichBand;
        case 'beta1'
            whichBand = beta1;
            bEEG.data = whichBand;
        case 'beta2'
            whichBand = beta2;
            bEEG.data = whichBand;
    end
    
    % Load regressors
    load( fullfile([dataFolder '/derivatives/behmod/sub-' ps{p}],['/sub-' ps{p} '_task-gnomes_reg.mat']) ,'barHeight','barHeightSplit','instReward','instPE','instRewardSplit','instRewardRiseFall','expectancy','instRewardSplit','expectancy','instRewardRiseFallSplit','barHeightRiseFallSplit')

    % x = [instReward'];
    if incBarHeight
        x = [barHeightSplit' instRewardSplit'];
        xCollapsed = [barHeight instReward];
    else
        if split
            x = [squeeze(instRewardRiseFallSplit(:,1,:))' squeeze(instRewardRiseFallSplit(:,2,:))'];
        else
            if excludeFall
                x = squeeze(instRewardRiseFallSplit(:,1,:))'; % Rising part only
            else
                 
                 x = instRewardSplit'; % Rising and falling parts
                 xPE = instPE';
            end

            xCollapsed = instReward';
        end
    end

    % Artifacts - use Unfold as ERPLab's basicrap includes overlapping
    % segments
    winrej = uf_continuousArtifactDetect(EEG,'amplitudeThreshold',150,'windowsize',2000,'stepsize',100,'combineSegments',[]);

    % Remove bad samples from model
    toRemove = [];
    isBad = zeros(size(x,1),1);
    for i = 1:size(winrej,1)
        toRemove = [toRemove winrej(i,1):winrej(i,2)];
    end
    isBad(toRemove) = 1;

    % Which samples are associated with non-zero regressors
    isNonZero = any(x~=0,2);

    % Figure out which proportion of modelled samples are artifacts
    artifactProp = mean(isBad & isNonZero);
    allArtifacts(p) = artifactProp;
    disp(artifactProp);

    bEEG.data(:,toRemove) = [];
    x(toRemove,:) = [];
    xCollapsed(toRemove,:) = [];
    xPE(toRemove,:) = [];

    % CORRELATION
    for i = 1:30
        tempBeta = [];
        tempBetaPE = [];
        for j = 1:size(x,2)
            tempBeta = [tempBeta corr(bEEG.data(i,:)',x(:,j))];
        end
        tempBetaPE = [tempBetaPE corr(bEEG.data(i,:)',xPE)];
        allCoeff(p,:,i) = tempBeta;
        allCoeffPE(p,i) = tempBetaPE;
        allCoeffCollapsed(p,i) = corr(bEEG.data(i,:)',xCollapsed);
    end
    
end

save(fullfile(resultsFolder,'corr_data.mat'));

%% Artifacts
disp('artifact proportion');
meanArtifacts = mean(allArtifacts);
sdArtifacts = std(allArtifacts);
tval = tinv(1-0.025, length(ps)-1);
ciArtifacts = tval * sdArtifacts / sqrt(length(ps));
disp([meanArtifacts meanArtifacts-ciArtifacts meanArtifacts+ciArtifacts]);