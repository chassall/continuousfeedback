% Compute instantaneous reward signals, etc.
% 
% Other m-files required: 
% EEGLAB toolbox https://github.com/sccn/eeglab
% /private/getpds

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Behavioural file headers
% [t thisGnomeType gnomeImageNumbeinstReward(thisGnomeType) thisTarget fixationTime rtStart rtEnd rexpectancyponseTime rexpectancypX rexpectancypY guexpectancysProp guexpectancysPoints totalPoints preFeedbackTime];
% 1 trial
% 2 gnome type: 1 = low, 2 = high, 3 = low/high, 4 = LOW/high, 5 = low/HIGH, 6 = uniform
% 3 image number
% 4 target level (proportion of bar, 0-1)
% 5 fixation time
% 6 rt start time
% 7 rt end time
% 8 total rexpectancyponse time
% 9 rexpectancyponse x
% 10 rexpectancyponse y
% 11 partipant guexpectancys (proportion of bar, 0-1)
% 12 points this round
% 13 total points
% 14 pre-feedback time (time before bar starts to grow)

close all; clear all; clc;
if ispc
    projectFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_Gnomes_Hassall';
else
    projectFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_Gnomes_Hassall';
end
dataFolder = [projectFolder '/data'];

ps = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};

allGuexpectancys = nan(length(ps),6,25);
meanGuexpectancys = nan(length(ps),6);
pbarHeight = getpds(); % Probability density for each gnome 1-6

for p = 1:length(ps)
   
   
    % Load behavioural data
    thisBehFolder = [dataFolder '/sub-' ps{p} '/beh'];
    thisBehFile = ['sub-' ps{p} '_task-gnomes_beh.tsv'];
    tempData = tdfread(fullfile(thisBehFolder,thisBehFile));
    participantData = [tempData.trial tempData.condition tempData.image tempData.target tempData.fixationTime tempData.rtStart tempData.rtEnd tempData.rt tempData.responseX tempData.responseY tempData.responseProp tempData.pointsThisRound tempData.pointsTotal tempData.preFeedbackTime];


    % Load EEG
    thisEEGFolder = [dataFolder '/derivatives/eegprep/sub-' ps{p}];
    thisEEGFile = ['sub-' ps{p} '_task-gnomes_eegprep.mat'];
    load(fullfile(thisEEGFolder,thisEEGFile));

    timexpectancytep = 1000/EEG.srate;
    for g = 1:6
       thisData = participantData(participantData(:,2) == g,:);
       meanGuexpectancys(p,g) = mean(thisData(:,11));
       allGuexpectancys(p,g,:) = thisData(:,11);
    end
    
    t = 1; % trial counter
    
    % Initialize regrexpectancysoinstReward
    barHeight = zeros(1,EEG.pnts);
    instReward = zeros(1,EEG.pnts);
    instRewardFull = zeros(1,EEG.pnts); % trials for which the bar reached the guess
    expectancy = zeros(1,EEG.pnts);
    instPE = zeros(1,EEG.pnts);
    aveRew = zeros(1,EEG.pnts);
    
    barHeightSplit = zeros(3,EEG.pnts);
    instRewardSplit = zeros(3,EEG.pnts);
    
    barHeightRiseFall = zeros(2,EEG.pnts);
    instRewardRiseFall = zeros(2,EEG.pnts); 
    
    barHeightRiseFallSplit = zeros(3,2,EEG.pnts);
    instRewardRiseFallSplit = zeros(3,2,EEG.pnts);
    

    % Loop over events
    for e = 1:length(EEG.event)
        thisEvent = EEG.event(e);
        if ismember(thisEvent.type,{'S 41','S 42','S 43','S 44','S 45','S 46'})
            
            thisGnomeType = participantData(t,2);
            thisPD = pbarHeight{thisGnomeType};
            thisTrialTarget = participantData(t,4);
            thisTrialGuess = participantData(t,11);
            startTime = round(thisEvent.latency);
            endTime = round(EEG.event(e+1).latency);
            
            prevData= participantData(1:t,:);
            prevData = prevData(prevData(:,2)==thisGnomeType,:);
            numEncounters = size(prevData,1);
            if size(prevData,1) == 1 % This is the first encounter with this gnome
                thisAveRew = 0.5;
            else % Compute the average reward for this gnome
                thesePoints = 1 - abs(prevData(1:numEncounters-1,11) - prevData(1:numEncounters-1,4));
                thisAveRew = mean(thesePoints);
                % thisAveRew = mean(prevData(1:numEncounters-1,12))/100;
            end
            %thisPE = prevData(numEncounters,12) - aveR;
            %peByGnome(thisGnomeType,numEncounters) = thisPE;
            %plot(peByGnome'); drawnow(); pause();

            numTimesteps = endTime-startTime+1;
            barInc = thisTrialTarget/numTimesteps;
            
            thisTrialI = startTime:endTime;
            
            
            barHeight(thisTrialI) = cumsum(barInc*ones(1,numTimesteps));
            
            %if thisTrialGuess > thisTrialTarget
                instReward(thisTrialI) = 1 - abs(thisTrialGuess - barHeight(thisTrialI));
                instPE(thisTrialI) = instReward(thisTrialI) - thisAveRew;
                aveRew(thisTrialI) = thisAveRew;
            %end

            expectancy(thisTrialI) = pdf(thisPD,barHeight(thisTrialI)');
            
            isRising = barHeight(thisTrialI) <= thisTrialGuess;
            isFalling = barHeight(thisTrialI) > thisTrialGuess;
            
%             if any(isRising) && any(isFalling)
%                 return;
%             end

            barHeightRiseFall(1,thisTrialI(isRising)) = barHeight(thisTrialI(isRising));
            instRewardRiseFall(1,thisTrialI(isRising)) = instReward(thisTrialI(isRising));
            barHeightRiseFall(2,thisTrialI(isFalling)) = barHeight(thisTrialI(isFalling));
            instRewardRiseFall(2,thisTrialI(isFalling)) = instReward(thisTrialI(isFalling));
            
            if ismember(thisEvent.type,{'S 41','S 42'})
                barHeightSplit(1,thisTrialI) = cumsum(barInc*ones(1,numTimesteps));
                instRewardSplit(1,thisTrialI) = 1 - abs(thisTrialGuess - barHeight(thisTrialI));
                barHeightRiseFallSplit(1,1,thisTrialI(isRising)) = barHeight(thisTrialI(isRising));
                instRewardRiseFallSplit(1,1,thisTrialI(isRising)) = instReward(thisTrialI(isRising));
                barHeightRiseFallSplit(1,2,thisTrialI(isFalling)) = barHeight(thisTrialI(isFalling));
                instRewardRiseFallSplit(1,2,thisTrialI(isFalling)) = instReward(thisTrialI(isFalling));
            elseif ismember(thisEvent.type,{'S 44','S 45'})
                barHeightSplit(2,thisTrialI) = cumsum(barInc*ones(1,numTimesteps));
                instRewardSplit(2,thisTrialI) = 1 - abs(thisTrialGuess - barHeight(thisTrialI));
                barHeightRiseFallSplit(2,1,thisTrialI(isRising)) = barHeight(thisTrialI(isRising));
                instRewardRiseFallSplit(2,1,thisTrialI(isRising)) = instReward(thisTrialI(isRising));
                barHeightRiseFallSplit(2,2,thisTrialI(isFalling)) = barHeight(thisTrialI(isFalling));
                instRewardRiseFallSplit(2,2,thisTrialI(isFalling)) = instReward(thisTrialI(isFalling));
            elseif ismember(thisEvent.type,{'S 43','S 46'})
                barHeightSplit(3,thisTrialI) = cumsum(barInc*ones(1,numTimesteps));
                instRewardSplit(3,thisTrialI) = 1 - abs(thisTrialGuess - barHeight(thisTrialI));
                barHeightRiseFallSplit(3,1,thisTrialI(isRising)) = barHeight(thisTrialI(isRising));
                instRewardRiseFallSplit(3,1,thisTrialI(isRising)) = instReward(thisTrialI(isRising));
                barHeightRiseFallSplit(3,2,thisTrialI(isFalling)) = barHeight(thisTrialI(isFalling));
                instRewardRiseFallSplit(3,2,thisTrialI(isFalling)) = instReward(thisTrialI(isFalling));
            end           
            
            t = t+1;
        end
    end
    
    expectancy = expectancy ./ max(expectancy);
    
    thisRegFolder = [dataFolder '/derivatives/behmod/sub-' ps{p}];
    thisRegFile = ['sub-' ps{p} '_task-gnomes_reg.mat'];
    if ~isfolder(thisRegFolder)
        mkdir(thisRegFolder);
    end
        
    save(fullfile(thisRegFolder,thisRegFile),'barHeight','barHeightSplit','barHeightRiseFall','instReward','instPE','aveRew','instRewardRiseFall','instRewardSplit','expectancy','instRewardRiseFallSplit','barHeightRiseFallSplit');
end