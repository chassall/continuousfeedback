% GLM analysis for the Gnomes project
% Collapse across cue types to verify existence of signal
% AND include a regressor for the (changing) average reward for each cue
% 
% Other m-files required: 
% EEGLAB toolbox https://github.com/sccn/eeglab
% Unfold toolbox: https://github.com/unfoldtoolbox/unfold
% /private/num2bv.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

close all; clear all; init_unfold();

% Analysis settings
% 1,0,0 include bar height, no regularization, no CV
% 1,1,1 include bar height, regularization, CV
incBarHeight = 1;
useReg = 0;
runCV = 0;

if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_Gnomes_Hassall\data';
else
    dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_Gnomes_Hassall/data'; % iMac
end

% Participant numbers
ps = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};

% 1-2 predictable
% 3 low/high
% 4 LOW/high
% 5 low/HIGH
% 6 uniform

% Make some conditions: 1 = predictable, 2 = somewhat predictable, 3 =
% unpredictable
cond1 = [1 2];
cond2 = [4 5];
cond3 = [3 6];

response = num2bv(31:36);
animStart = num2bv(41:46);
animEnd = num2bv(51:56);
response1 = num2bv(30 + cond1);
animStart1 = num2bv(40+cond1);
animEnd1 = num2bv(50+cond1);
response2 = num2bv(30+cond2);
animStart2 = num2bv(40+cond2);
animEnd2 = num2bv(50+cond2);
response3 = num2bv(30+cond3);
animStart3 = num2bv(40+cond3);
animEnd3 = num2bv(50+cond3);


srate = 250;

respTimeLim = [-0.5,0];
respPntLim = srate * respTimeLim;
respPnts= respPntLim(1):respPntLim(2);
numRespPnts = length(respPnts);
respBL = [-0.8 -0.6]; % Baseline, in seconds

animStartTimeLim = [0,0.8];
animStartPntLim = srate * animStartTimeLim;
animStartPnts= animStartPntLim(1):animStartPntLim(2);
numAnimStartPnts = length(animStartPnts);
animStartBL = [-0.2 0]; % Baseline, in seconds

animEndTimeLim = [0,0.8];
animEndPntLim = srate * animEndTimeLim;
animEndPnts= animEndPntLim(1):animEndPntLim(2);
numAnimEndPnts = length(animEndPnts);
animEndBL = [-0.2 0]; % Baseline, in seconds


% if incBarHeight
%     % allBeta = nan(length(ps),13917,30);
%     allBeta = nan(length(ps),10182,30);
% else
%     allBeta = nan(length(ps),9729,30);
% end

allX = {};
allArtifactProp = [];

if useReg
    % Cross-validation
    if runCV
        lambdas = [0.001 0.01 0 1 10 100 1000 10000 100000 1000000 10000000];
        lambdas = [100 1000 10000 100000 1000000];
        k = 10;
        allCVErrors = [];
    else
        load('cv_results_collapsed.mat','lambdas','allCVErrors');
    end
else
    allCVErrors = [];
end

% Loop through participants
for p = 1:length(ps)
    disp(ps{p});
    
    % Load preprocessed EEG
    prepFile = ['sub-' ps{p} '_task-gnomes_eegprep.mat'];
    prepFolder = [dataFolder '/derivatives/eegprep/sub-' ps{p}];
    load(fullfile(prepFolder,prepFile));
    
    % Round latencies, as some may be non-integers due to resampling
    for i = 1:length(EEG.event)
        EEG.event(i).latency = round(EEG.event(i).latency);
    end
      %%

    % Load regressors
    load( fullfile([dataFolder '/derivatives/behmod/sub-' ps{p}],['/sub-' ps{p} '_task-gnomes_reg.mat']) ,'aveRew','barHeight','barHeightSplit','instReward','instRewardSplit','instRewardRiseFall','expectancy','instRewardSplit','expectancy','instRewardRiseFallSplit','barHeightRiseFallSplit')
     
    %actualRewards = aveRew(aveRew ~= 0);
    %aveRew = aveRew - mean(actualRewards);


    whichRewardRF = instRewardRiseFall; % All conditions
    whichRewardRF1 = squeeze(instRewardRiseFallSplit(1,:,:));
    whichRewardRF2 = squeeze(instRewardRiseFallSplit(2,:,:));
    whichRewardRF3 = squeeze(instRewardRiseFallSplit(3,:,:));

    whichBarRF = barHeightSplit;
    whichBarRF1 = squeeze(barHeightSplit(1,:,:));
    whichBarRF2 = squeeze(barHeightSplit(2,:,:));
    whichBarRF3 = squeeze(barHeightSplit(3,:,:));
    
    % Bar/reward step size
    barDeltas = nonzeros(diff(barHeight));
    barDelta = barDeltas(1);
    barDelta = 6.6861e-04;
    
    barSignal = 0:barDelta:1;
    
    fallingBar = 1:-barDelta:0;
    risingBar = flip(fallingBar);
    rewSignal = [risingBar(1:end-1) fallingBar];
    
    % Fixed-time components
    respX = sparse(EEG.pnts,numRespPnts);
    animStartX = sparse(EEG.pnts,numAnimStartPnts);
    animEndX = sparse(EEG.pnts,numAnimEndPnts);

    for iEvent = 1:length(EEG.event)
        thisLatency = EEG.event(iEvent).latency;
        switch EEG.event(iEvent).type
            case animStart
                for j = 1:numAnimStartPnts
                    animStartX(thisLatency+animStartPnts(j),j) = 1;
                end
            case animEnd
                for j = 1:numAnimEndPnts
                    animEndX(thisLatency+animEndPnts(j),j) = 1;
                end
        end
    end
    
    % Bar/Reward components
    barX = sparse(EEG.pnts,length(barSignal));
    rewRiseX = sparse(EEG.pnts,length(risingBar));
    rewFallX = sparse(EEG.pnts,length(fallingBar));
    aveRewRiseX = sparse(EEG.pnts,length(risingBar));
    aveRewFallX = sparse(EEG.pnts,length(fallingBar));

    % Bar height signal
    for i = 1:size(whichBarRF,2)
        if whichBarRF(1,i) ~= 0
            whichPoint = dsearchn(barSignal',whichBarRF(1,i));
            barX(i,whichPoint) = 1;
        end
    end
    
    % Reward signal
    for i = 1:size(whichRewardRF,2)
        if whichRewardRF(1,i) ~= 0
            whichPoint = dsearchn(risingBar',whichRewardRF(1,i));
            rewRiseX(i,whichPoint) = 1;
        end
    end
    for i = 1: size(whichRewardRF,2)
        if whichRewardRF(2,i) ~= 0
            whichPoint = dsearchn(fallingBar',whichRewardRF(2,i));
            rewFallX(i,whichPoint) = 1;
        end
    end

    % Ave reward signal
    for i = 1:size(whichRewardRF,2)
        if whichRewardRF(1,i) ~= 0
            whichPoint = dsearchn(risingBar',whichRewardRF(1,i));
            aveRewRiseX(i,whichPoint) = aveRew(i);
        end
    end
    for i = 1: size(whichRewardRF,2)
        if whichRewardRF(2,i) ~= 0
            whichPoint = dsearchn(fallingBar',whichRewardRF(2,i));
            aveRewFallX(i,whichPoint) = aveRew(i);
        end
    end

    if incBarHeight
        X = [rewRiseX rewFallX aveRewRiseX aveRewFallX animStartX animEndX];
    else
        X = [rewRiseX rewFallX aveRewRiseX aveRewFallX animEndX];
    end
    
    % Keep a record of all DMs, e.g. to calculate VIF
    allX{p} = X;

    % Indices into beta matrix
    rewSignalLength = length(risingBar) + length(fallingBar);
    betaI = {};
    
    if incBarHeight
        betaI{1}= 1:rewSignalLength;
        betaI{2} = (betaI{1}(end)+1):(betaI{1}(end)+rewSignalLength);
        betaI{3}= (betaI{2}(end)+1):(betaI{2}(end)+numAnimStartPnts);
        betaI{4}= (betaI{3}(end)+1):(betaI{3}(end)+numAnimEndPnts);
    else
        betaI{1}= 1:rewSignalLength;
        betaI{2}= (betaI{1}(end)+1):(betaI{1}(end)+rewSignalLength);
        betaI{3}= (betaI{2}(end)+1):(betaI{2}(end)+numAnimEndPnts);
    end
    if incBarHeight
        breakPoints = [betaI{1}(end) betaI{2}(end) betaI{3}(end)];
    else
        breakPoints = [betaI{1}(end) betaI{2}(end)];
    end

   
    % Remove zero rows
    nonZero = any(X,2);
    isZero = ~nonZero;
    
    % Check for artifacts
    isArtifact = zeros(size(isZero));
    winrej = uf_continuousArtifactDetect(EEG,'amplitudeThreshold',150,'windowsize',2000,'stepsize',100,'combineSegments',[]);
    
    % Remove bad samples from model
    toRemove = [];
    for i = 1:size(winrej,1)
        toRemove = [toRemove winrej(i,1):winrej(i,2)];
    end  
    isArtifact(toRemove) = 1;
    
    % Number of artifact as a proportion of samples of interest
    allArtifactProp(p) = mean(isArtifact & ~isZero)

    % Remove artifacts and non-zero rows
    X(isArtifact | isZero,:) = [];
    EEG.data(:,isArtifact | isZero) = [];
    EEG.pnts = size(EEG.data,2);

    % De-mean the average reward signal. This is a bit complicated because
    % the reward signals are not aligned (see design matrix)
     aveRewX =full(X(:,betaI{2})); % Ave reward part, convert to full
     isZero = aveRewX == 0; 
     aveRewX(isZero) = NaN; % Replace 0 with Nan
     meanAveRewX = mean(aveRewX,'omitnan'); % Mean of each column
     aveRewX = aveRewX - meanAveRewX; % Demean
     aveRewX(isnan(aveRewX)) = 0; % Replace NaN with 0
     aveRewX = sparse(aveRewX); % Back to sparse
    
     % Replace ave reward regressors with demeaned ave reward regressors
     X(:,betaI{2}) = aveRewX;

%     isNonZero =X(:,betaI{1})==1;
%     for i = 1:length(betaI{2})
%         thisColI = full(isNonZero(:,i));
%         thisCol = aveRewX(:,i);
%         thisMean = mean(thisCol(thisColI));
%         if ~isnan(thisMean)
%         thisCol = thisCol - thisMean;
%         aveRewX(:,i) = thisCol;
%         end
%     end

    if useReg
        % Solve with regularization
        % Need to split by condition to compute
        % Should be OK as conditions don't overlap
        disp('solving GLM with regularization');
        tic;

        if runCV
            regtype = 'onediff';
            condIs = {1:size(X,2)};
            whichBreakpoints = {breakPoints};
            [theseErrors,bestBeta]  = doRegCV(EEG.data,X,regtype,condIs,whichBreakpoints,lambdas,k);
            plot(theseErrors); drawnow();
            allBeta(p,:,:) = bestBeta;
            allCVErrors(p,:) = theseErrors;
        else
            theseErrors = allCVErrors(p,:);
            [~,bestI] = min(theseErrors);
            bestLambda = lambdas(bestI);
            whichBreakpoints = breakPoints;
            thisPDM = pinv_reg(X,bestLambda,'onediff',whichBreakpoints);
            allBeta(p,:) = thisPDM * EEG.data';
        end


        toc
    else
        lsmriterations = 400;
        [allBeta(p,:,:),~,~] = lsmr(X,double(EEG.data'),[],10^-8,10^-8,[],lsmriterations);
    end

    % Save this participant's data
    saveFile = ['sub-' ps{p} '_task-gnomes_glmcollapsedaverew_' num2str(incBarHeight) '_' num2str(useReg) '.mat'];
    saveFolder = [dataFolder '/derivatives/glmres/sub-' ps{p}];
    if ~exist(saveFolder)
        mkdir(saveFolder)
    end

    % To save
    chanlocs = EEG.chanlocs;
    srate = EEG.srate;
    beta = squeeze(allBeta(p,:,:));
    artifactProp = allArtifactProp(p);
    if useReg
        cvErrors = allCVErrors(p,:);
    else
        cvErrors = [];
    end
    X = allX{p};
    save(fullfile(saveFolder,saveFile),'chanlocs','srate','betaI','X','beta','artifactProp','cvErrors');

end