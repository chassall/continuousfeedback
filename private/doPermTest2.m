function [clusterInfo] = doPermTest2(data,alpha,neighbours,EEG)
%DOPERMTEST Do permuation testing
%   Do a t-test at each point. Look for contiguous sig.
%   values according to the given alpha value. Sum the t-values to get the
%   cluster mass. Do this for the actual data, then for permuted data.
%   Permute data by randomly inverting (in vertical axis) each
%   participant's waveform. 
% data: participants X channels X samples
% alpha: sign. level
%
% permClusterMasses: max sig. cluster value after permutation
% actualclusterMasses: cluster values for actual data
% window: start/end point of max cluster (greatest cluster mass)
% 

% Permutations (for 10 or fewer participants just do every combination)
if size(data,1) <= 10
    allPerms = dec2bin(0:2^10-1)-'0'; % Makes all binary combinations
    allPerms(allPerms == 0) = -1;
    allPerms = allPerms';
else
    numPerms = 1000;
    allPerms = -1 + 2*(randi(2,size(data,1),numPerms)-1);
end

% Compute actual t values for each channel and time
[actualhs, ps, cis, stats] = ttest(data,zeros(size(data)),'Alpha',alpha);
actualhs = squeeze(actualhs);
actualts = squeeze(stats.tstat);

% Now loop through each channel and compute a cluster statistic based on
% the temporal and spatial adjacency
clusterInfo = {};
for i = 1:length(neighbours)

    % This channel
    thisLabel = neighbours(i).label;
    clusterInfo{i,1} = thisLabel;
    thisI = eeg_chaninds(EEG,thisLabel);
    thisH = actualhs(thisI,:);
    thisT = actualts(thisI,:);

    % Neighbouring channels
    theseNeighbours = neighbours(i).neighblabel;
    neighbI = eeg_chaninds(EEG,theseNeighbours);
    neighbH = actualhs(neighbI,:);
    neighbT = actualts(neighbI,:);

    allH = [thisH; neighbH];
    commonH = all(allH,1);

    allT = [thisT; neighbT];

    % Find temporal clusters
    isDiff = [true; diff(commonH(:)) ~= 0];
    commonH = [commonH 0];
    changeInd = find([isDiff' true]);
    rejectNull = commonH(changeInd);
    actualClusterExtents = diff(changeInd);
    actualClusterExtents = actualClusterExtents(rejectNull(1:end-1) == 1); % Only count clusters where h = 1
    
    actualclusterStartI = changeInd(rejectNull == 1);
    actualclusterEndI = actualclusterStartI + actualClusterExtents - 1;
    actualclusterMasses = [];
    for j = 1:length(actualClusterExtents)
        actualclusterMasses(j) = mean(sum(allT(:,actualclusterStartI(j):actualclusterEndI(j)),2),1);
    end
    [actualMaxCluster,maxI] = max(abs(actualclusterMasses));
    scaledWindow = [actualclusterStartI(maxI) actualclusterEndI(maxI)];

    clusterInfo{i,2} = actualMaxCluster;
    clusterInfo{i,3} = scaledWindow;
    clusterInfo{i,4} = scaledWindow./size(data,3);
end

%% Permutation tests (new)

allPerms = -1 + 2*(randi(2,size(data,1),numPerms)-1);
clusterMasses = [];
for i = 1:size(allPerms,2)
    thisPerm = allPerms(:,i);
    flippedBetas = data .* thisPerm;
    [hs, ps, cis, stats] = ttest(flippedBetas,zeros(size(flippedBetas)),'Alpha',alpha);
    
    hs = squeeze(hs);
    tstats = squeeze(stats.tstat);

    for i = 1:length(neighbours)

        % This channel
        thisLabel = neighbours(i).label;
        thisI = eeg_chaninds(EEG,thisLabel);
        thisH = hs(thisI,:);
        thisT = actualts(thisI,:);

        % Neighbouring channels
        theseNeighbours = neighbours(i).neighblabel;
        neighbI = eeg_chaninds(EEG,theseNeighbours);
        neighbH = hs(neighbI,:);

        allH = [thisH; neighbH];
        commonH = all(allH,1);

        allT = [thisT; neighbT];

        % Find temporal clusters
        isDiff = [true; diff(commonH(:)) ~= 0];
        commonH = [commonH 0];
        changeInd = find([isDiff' true]);
        rejectNull = commonH(changeInd);
        actualClusterExtents = diff(changeInd);
        actualClusterExtents = actualClusterExtents(rejectNull(1:end-1) == 1); % Only count clusters where h = 1

        actualclusterStartI = changeInd(rejectNull == 1);
        actualclusterEndI = actualclusterStartI + actualClusterExtents - 1;
        theseClusterMasses = [];
        for j = 1:length(actualClusterExtents)
            theseClusterMasses(j) = mean(sum(allT(:,actualclusterStartI(j):actualclusterEndI(j)),2),1);
        end
        [actualMaxCluster,maxI] = max(abs(theseClusterMasses));

        if isempty(actualMaxCluster)
            clusterMasses = [clusterMasses 0];
        else
            clusterMasses = [clusterMasses actualMaxCluster];
        end
    end

end

% Look at the distribution of clustermasses (unused)
figure();
histogram(clusterMasses);
for i = 1:size(clusterInfo,1)
   if ~isempty(clusterInfo{i,2})
   xline( clusterInfo{i,2}); hold on;
   end
end

% Check our clusters against the permutation
for i = 1:size(clusterInfo,1)
    if ~isempty(clusterInfo{i,2})
        clusterInfo{i,5} = sum(clusterMasses < clusterInfo{i,2}) / length(clusterMasses);
        clusterInfo{i,6} = 1 - clusterInfo{i,5};
    end 
end

% Display p value of the target clusgter.
% Can't use the index from the EEG file because it might be
% different from what's in clusterInfo.
% thisTSIndex = find(strcmp(tsChannelString,clusterInfo(:,1))); 
% disp('Cluster Info');
% disp(clusterInfo(thisTSIndex,:));

% disp('Neighbours');
% disp(neighbours(thisTSIndex).label);
% disp(neighbours(thisTSIndex).neighblabel);
% clusterIs = [eeg_chaninds(EEG,neighbours(thisTSIndex).label) eeg_chaninds(EEG,neighbours(thisTSIndex).neighblabel)];


% % Actual t-stats
% [actualhs, ps, cis, stats] = ttest(data,zeros(size(data)),'Alpha',alpha);
% tstats = stats.tstat;
% isSig = find(actualhs);
% isDiff = [true; diff(actualhs(:)) ~= 0];
% actualhs = [actualhs 0];
% changeInd = find([isDiff' true]);
% rejectNull = actualhs(changeInd);
% actualClusterExtents = diff(changeInd);
% actualClusterExtents = actualClusterExtents(rejectNull(1:end-1) == 1); % Only count clusters where h = 1
% 
% actualclusterStartI = changeInd(rejectNull == 1);
% actualclusterEndI = actualclusterStartI + actualClusterExtents - 1;
% actualclusterMasses = [];
% actualClusterHeights = [];
% for i = 1:length(actualClusterExtents)
%     actualclusterMasses(i) = sum(tstats(actualclusterStartI(i):actualclusterEndI(i)));
%     actualClusterHeights(i) = max(abs(tstats(actualclusterStartI(i):actualclusterEndI(i))));
% end
% [actualMaxCluser,maxI] = max(abs(actualclusterMasses));
% %window = [actualclusterStartI(maxI) actualclusterEndI(maxI)];
% window = [actualclusterStartI; actualclusterEndI];
% 
% [actualMaxClusterHeight,maxI] = max(actualClusterHeights);
% 
% % subplot(1,2,1); plot(ts); hold on; plot(isSig,zeros(1,length(isSig)),'.');
% 
% % Permutations (for 10 or fewer participants just do every combination)
% if size(data,1) <= 10
%     allPerms = dec2bin(0:2^10-1)-'0'; % Makes all binary combinations
%     allPerms(allPerms == 0) = -1;
%     allPerms = allPerms';
% else
%     numPerms = 1000;
%     allPerms = -1 + 2*(randi(2,size(data,1),numPerms)-1);
% end
% 
% permClusterExtents = [];
% permClusterMasses = [];
% permClusterHeights = [];
% for i = 1:size(allPerms,2)
%     thisPerm = allPerms(:,i);
%     flippedBetas = data .* thisPerm;
%     [hs, ps, cis, stats] = ttest(flippedBetas,zeros(size(flippedBetas)),'Alpha',alpha);
%     tstats = stats.tstat;
%     
%     isDiff = [true; diff(hs(:)) ~= 0];
%     hs = [hs 0];
%     changeInd = find([isDiff' true]);
%     rejectNull = hs(changeInd);
%     runLengths = diff(changeInd);
%     runLengths = runLengths(rejectNull(1:end-1) == 1); % Only count clusters where h = 1
%     permClusterExtents = [permClusterExtents runLengths];
%     
%     clusterStartI = changeInd(rejectNull == 1);
%     clusterEndI = clusterStartI + runLengths - 1;
%     clusterMasses = [];
%     clusterHeights = [];
%     for i = 1:length(runLengths)
%         clusterMasses(i) = sum(tstats(clusterStartI(i):clusterEndI(i)));
%         clusterHeights(i) = max(abs(tstats(clusterStartI(i):clusterEndI(i))));
%     end
%     
%     [~,J] = max(abs(clusterMasses));
%     if isempty(clusterMasses)
%         permClusterMasses = [permClusterMasses 0];
%     else
%         permClusterMasses = [permClusterMasses clusterMasses(J)];
%     end
%     
%     [~,J] = max(clusterHeights);
%     permClusterHeights = [permClusterHeights clusterHeights(J)];
% end

% % Check our clusters against the permutation
% p = [];
% for i = 1:length(actualclusterMasses)
%     p(i) = 1 - sum(abs(permClusterMasses) < abs(actualclusterMasses(i))) / length(permClusterMasses);
% end
% disp(['Cluster ts: ' num2str(abs(actualclusterMasses))]);
% disp(['Cluster ps: ' num2str(p)]);
% 
% disp('(Points)');
% disp(num2str(window));
% 
% % Greatest t cluster?
% disp(['(Points) From ' num2str(actualclusterStartI(maxI)) ' to ' num2str(actualclusterEndI(maxI))]);

end

