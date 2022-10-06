% Analyze and plot behavioural data for the continuous feedback processing project

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Behavioural file headers
% [t thisGnomeType gnomeImageNumbers(thisGnomeType) thisTarget fixationTime rtStart rtEnd responseTime respX respY guessProp guessPoints totalPoints preFeedbackTime];
% 1 trial
% 2 gnome type: 1 = low, 2 = high, 3 = low/high, 4 = LOW/high, 5 = low/HIGH, 6 = uniform
% 3 image number
% 4 target level (proportion of bar, 0-1)
% 5 fixation time
% 6 rt start time
% 7 rt end time
% 8 total response time
% 9 response x
% 10 response y
% 11 partipant guess (proportion of bar, 0-1)
% 12 points this round
% 13 total points
% 14 pre-feedback time (time before bar starts to grow)

% For each participant and trial, compute the following time series (200 Hz):
% bar level
% reward level
% prediction error (using a normative model for now)

close all; clear all; clc;
if ispc
    projectFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_Gnomes_Hassall';
else
    projectFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_Gnomes_Hassall';
end
dataFolder = [projectFolder '/data'];
figuresFolder = [projectFolder '/figures'];
resultsFolder = [projectFolder '/analysis/results'];

if ~exist(resultsFolder)
    mkdir(resultsFolder);
end

ps = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'}; 
whichPs = [1:10 12:21]; % 11 noisy

% Swap the default order for plotting
whichTrials = 1:25;

allGuess = nan(length(ps),6,length(whichTrials));
meanGuess = nan(length(ps),6);
stdGuess = nan(length(ps),6);
allDist = nan(length(ps),6,length(whichTrials));
meanDist = nan(length(ps),6);
meanPoints = nan(length(ps),6);
allPoints = nan(length(ps),6,length(whichTrials));
targets = [1/3 2/3 1/2 1/3 2/3 1/2];
allTargets = nan(length(ps),6,length(whichTrials));
allBonus = [];
outcomes = [];
gnomeTypes = [];
allPointTotals = [];
for p = 1:length(ps)

    thisBehFolder = [dataFolder '/sub-' ps{p} '/beh'];
    thisBehFile = ['sub-' ps{p} '_task-gnomes_beh.tsv'];
    tempData = tdfread(fullfile(thisBehFolder,thisBehFile));
    participantData = [tempData.trial tempData.condition tempData.image tempData.target tempData.fixationTime tempData.rtStart tempData.rtEnd tempData.rt tempData.responseX tempData.responseY tempData.responseProp tempData.pointsThisRound tempData.pointsTotal tempData.preFeedbackTime];

    for g = 1:6
        thisData = participantData(participantData(:,2) == g,:);
        thisMeanOutcome = mean(thisData(whichTrials,4));
        meanGuess(p,g) = mean(thisData(whichTrials,11));
        stdGuess(p,g) = std(thisData(whichTrials,11));
        allGuess(p,g,:) = thisData(whichTrials,11);
        meanPoints(p,g) = 1 - (mean(abs(thisData(whichTrials,11) - mean(thisData(whichTrials,4)))));
        meanDist(p,g) = mean((thisData(whichTrials,11) - targets(g)));
        allTargets(p,g,:) = thisData(whichTrials,4);
        allPoints(p,g,:) = 1 - (abs(thisData(whichTrials,11) - thisData(whichTrials,4)));
        allDist(p,g,:) = abs(thisData(whichTrials,11) - thisMeanOutcome);
    end

    groupings = [1 2; 4 5; 3 6];
    boundaries = [];
    
    % One set of quantiles for all gnomes
    thesePoints = allPoints(p,:,:);
    thesePoints = thesePoints(:);
    boundaries(1,:) = median(thesePoints);
    boundaries(2,:) = boundaries(1,:);
    boundaries(3,:) = boundaries(1,:);

    % loop over trials
    for t = 1:150
        thisTrial = participantData(t,:);
        
        thisOutcome = 1 - abs(thisTrial(11) - thisTrial(4));

        % Collapse across gnomes
        switch thisTrial(2)
            case {1, 2}
                thisGnomeType = 1;
                thisBoundary = boundaries(1,:);
            case {4, 5}
                thisGnomeType = 2;
                thisBoundary = boundaries(2,:);
            case {3, 6}
                thisGnomeType = 3;
                thisBoundary = boundaries(3,:);
        end

        if thisOutcome < thisBoundary(1)
            outcomeCond = 1;
        else
            outcomeCond = 2;
        end

        outcomes(p,t) = outcomeCond;
        gnomeTypes(p,t) = thisGnomeType;
    end

    thesePoints = thisData(end,13);
    allPointTotals(p,1) = thesePoints;
    thisBonus = thesePoints * .0002;
    allBonus(p,1) = thisBonus;
end

% Save behavioural results
save(fullfile(resultsFolder,'behResults.mat'),'outcomes','gnomeTypes','allGuess','allDist');

disp(['Mean points ' num2str(mean(allPointTotals))]);
disp(['SD points ' num2str(std(allPointTotals))]);

disp(['Mean bonus ' num2str(mean(allBonus))]);
disp(['SD bonus ' num2str(std(allBonus))]);