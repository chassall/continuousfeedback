% Figures for the Gnomes project
% 
% Other m-files required: 
% figset.m (common stats and figure settings)
% EEGLAB toolbox https://github.com/sccn/eeglab
% brewermap: https://github.com/DrosteEffect/BrewerMap
% subtightplot: https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
% cspy: https://uk.mathworks.com/matlabcentral/fileexchange/46551-cspy-m
% /private/doPermTest.m
% boundedline.m https://github.com/kakearney/boundedline-pkg
% /private/makefigure.m
% /private/formatNBP.m
% https://uk.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m
% distributionPlot.m https://uk.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m
% https://github.com/raacampbell/notBoxPlot

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com


%% Figure 1c - choice distributions
close all; clear all; figset;
load(fullfile(resultsFolder,'behResults.mat'),'allGuess','allDist');

gnomeOrder = [1 2 4 5 3 6]; % Swap gnome order for plotting (group by expectancy)
gnomeLabels = {'low','high','low/high','LOW/high','low/HIGH','uniform'};

makefigure(9,6);
combinedGuesses = [];
for p = 1:length(whichPs)
    thisP = whichPs(p);
    combinedGuesses = [combinedGuesses squeeze(allGuess(thisP,:,:))];
end
combinedGuesses = combinedGuesses';
xs = 1:6;
distX = 0.55:1:5.55;
halfX = 0.5;
% myColormap = cbrewer('seq','GnBu',256);
myColormap = brewermap(256,'GnBu');

myColormap = flip(myColormap);
dp = distributionPlot(combinedGuesses(:,gnomeOrder),'xValues',xs,'colormap',myColormap,'showMM',0,'histOri','right','histOpt',1);
ylim([0,1]);

hold on;

probs = [1 0; 0 1; 0.5 0.5; 0.8 0.2; 0.2 0.8; nan nan];

yline(1/3,'LineStyle',':');
yline(2/3,'LineStyle',':');

ax = gca;
ax.XLabel.String = 'Gnome type';
ax.YLabel.String = 'Guess (proportion of bar)';
ax.XTickLabel = gnomeLabels(gnomeOrder);
ax.XTickLabel = [];
ax.XLabel.Visible = 0;
ax.FontName = fontName;
ax.FontSize = fontSize;

pause(0.1);
print(fullfile(figuresFolder,'fig_01c_choicedist.tiff'),'-dtiff','-r600');

%% Figure 1d - distance by trial and gnome
close all; clear all; figset;
load(fullfile(resultsFolder,'behResults.mat'),'allGuess','allDist');

plotLineColours = allColours([2 5 4 3],:);
lineWidth = 1.5;
fontSize = 8;
fontName = 'Arial';

groupings = [1 2; 4 5; 3 6]; % Group by expectancy
mmDist = movmean(allDist,5,3);

makefigure(8,6);
collapsedMean = [];
for i = 1:size(groupings,1)

    thisMean = squeeze(mean(mmDist(whichPs,groupings(i,:),:),2));
    collapsedMean(:,i) = mean(thisMean,2);
    thisGrandMean = mean(thisMean,1);
    thisGrandSD = std(thisMean,[],1);
    tval = tinv(0.025,length(whichPs)-1);
    thisGrandCI = tval * thisGrandSD ./sqrt(length(whichPs));
    plot(thisGrandMean,'LineWidth',lineWidth,'Color',plotLineColours(i,:)); hold on;
    %boundedline(1:25,thisGrandMean,thisGrandCI,'cmap',plotLineColours(i,:),'alpha'); hold on;

end
xlim([1 25]);
xlabel('Trial number');
ylabel('Distance to M_o_u_t_c_o_m_e (proportion of bar)');
legend('high','medium','low','Box','off','Location','SouthWest');
box off;

ax = gca;
ax.FontName = fontName;
ax.FontSize = fontSize;

% Stats
disp('Behavioural ANOVA')
between = array2table(collapsedMean,'VariableNames',{'V1','V2','V3'});
pCond = {'1'; '2'; '3'};
within = table(pCond,'VariableNames',{'P'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V3 ~ 1','WithinDesign',within); % Use ~1 since no betwee-subject variable
[ranovatbl,A,C,D] = ranova(rm,'WithinModel','P');
disp(ranovatbl);
etap = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4))
etag = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4) + ranovatbl.SumSq(2))

meanDist = mean(collapsedMean,1);
stdDist = std(collapsedMean,[],1);
tval = tinv(1-0.025,length(whichPs)-1);
ciDist = tval*stdDist./sqrt(length(whichPs));
disp('mean and ci');
disp(meanDist);
disp(meanDist - ciDist);
disp(meanDist + ciDist);

print(fullfile(figuresFolder,'fig_01d_distbytrial.tiff'),'-dtiff','-r600');

%% Fig 2a - sample reward signal for correlational analysis
close all; clear all; figset;
load(fullfile(dataFolder,'derivatives','behmod','sub-01','sub-01_task-gnomes_reg.mat'),'instReward');

pnts = 10000:20000;
times = (pnts / 250);
makefigure(13,3);
plot(times,100*instReward(pnts),'LineWidth',lineWidth,'Color','k');
ylabel('Reward (points)');
xlabel('Time (s)');
set(gca,'FontSize',fontSize);
set(gca,'FontName',fontName);
set(gca,'Box','off');
print(fullfile(figuresFolder,'fig_02a_instreward.tiff'),'-dtiff','-r600');

%% Fig 2b - correlational topography
close all; clear all; figset;
load(fullfile(resultsFolder,'corr_data.mat'),'allCoeffCollapsed','allArtifacts');
load(fullfile(dataFolder,'derivatives','eegprep','sub-01','sub-01_task-gnomes_eegprep.mat'),'EEG');

% Artifacts
disp(['P11 (excluded) artifacts: ' num2str(allArtifacts(11))]);
disp('Remaining artifacts mean and sd');
mean(allArtifacts(whichPs))
std(allArtifacts(whichPs))

meanAllCoeffCollapsed = squeeze(mean(allCoeffCollapsed(whichPs,:),1));

[~,minI] = min(meanAllCoeffCollapsed);
disp('min correlation');
disp(EEG.chanlocs(minI).labels);
maxCoeff = max(max(meanAllCoeffCollapsed));
minCoeff = min(min(meanAllCoeffCollapsed));
topoLim = max([abs(maxCoeff) abs(minCoeff)]);

makefigure(5,4);
tp = topoplot(meanAllCoeffCollapsed,EEG.chanlocs,'maplimits' ,[-topoLim,topoLim],'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
colormap(myColormap);
% colormap('parula');
c = colorbar();
%c.Location = 'EastOutside';
c.Label.String = 'Correlation coeff. (r)';
c.Label.FontName = fontName;
c.Label.FontSize = fontSize;
set(gca,'FontSize',fontSize);
set(gca,'FontName',fontName);
set(gca,'Box','off');

% 
disp('mean and CI of Pearson''s r')
meanR = mean(allCoeffCollapsed(whichPs,minI));
stdR = std(allCoeffCollapsed(whichPs,minI));
tval = tinv(1-0.025,length(whichPs)-1);
ciR = tval * stdR ./ sqrt(length(whichPs));
disp(meanR);
disp(meanR - ciR);
disp(meanR + ciR);

% t-test at every electrode
allT = [];
allP = [];
for i = 1:30
    [h,p, ci, stats] = ttest(allCoeffCollapsed(whichPs,i));

    allT(i) = stats.tstat;
    allP(i) = p;
end

disp('stats at peak')
disp(allT(minI));
disp(allP(minI));
disp('Cohen''s d');
disp(meanR/stdR);
bonfAlpha = 0.001/30

print(fullfile(figuresFolder,'fig_02b_corrtopo.tiff'),'-dtiff','-r600');

%% Fig 2c - sample design matrix
close all; clear all; figset;
subplot = @(m,n,i) subtightplot (m, n, i, [0.01 0.01], [0.01 0.01], [0.01 0.01]);

% Load and plot a sample design matrix and some other info
sub1data = load(fullfile(dataFolder,'derivatives/glmres/sub-01/sub-01_task-gnomes_glmcollapsedaverew_1_1.mat'));
load(fullfile(dataFolder,'derivatives/eegprep/sub-01/sub-01_task-gnomes_eegprep.mat'),'EEG');

halfRewTime = 0.5*(sub1data.betaI{1}(end)-1) / 250;
rewTimes = -halfRewTime:(1/250):halfRewTime;
xs = 1:(sub1data.betaI{1}(end)-1);
pnts = 9000:100000;

axs = [];

makefigure(6,8);
axs(1) = subplot(1,6,1);
plot(EEG.data(20,pnts),pnts,'Color','k');
ax = gca;
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'off';

axs(2) = subplot(1,6,2:6);
cspy(sub1data.X(pnts,xs),'markersize',5,'colormap','turbo','Levels',100);
ax = gca;
currentTicks = ax.XTick;
newTicksS = [-4 -2 0 2 4];
newTicksP = dsearchn(rewTimes',newTicksS');
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
set(gca,'FontSize',fontSize);
set(gca,'FontName',fontName);
set(gca,'Box','off');
linkaxes(axs,'y');

print(fullfile(figuresFolder,'fig_02c_sampledesignmatrix.tiff'),'-dtiff','-r600');


%% Fig 4 - sample design matrix number 2
close all; clear all; figset;
sub1data = load(fullfile(dataFolder,'derivatives/glmres/sub-01/sub-01_task-gnomes_glmcollapsedaverew_1_1.mat'));

% Settings
pnts = 9000:100000;
myColourMap = brewermap(10,'YlGnBu');

makefigure(12,12);
xs = [sub1data.betaI{1} sub1data.betaI{2}];
cspy(sub1data.X(pnts,xs),'markersize',5,'colormap','turbo','Levels',100);
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
set(gca,'FontSize',fontSize);
set(gca,'FontName',fontName);
set(gca,'Box','off');
print(fullfile(figuresFolder,'fig_04_sampledesignmatrix2.tiff'),'-dtiff','-r600');

%% Load CV results and display best lambda
close all; clear all; figset();
load(fullfile(resultsFolder,'cv_results.mat'));
theseCVErrors = allCVErrors(whichPs,:);
[~,i] = min(theseCVErrors,[],2)

% theseCVErrors = theseCVErrors - mean(theseCVErrors,2);
theseCVErrors = normalize(theseCVErrors,2);

makefigure(8,6);
nbp = notBoxPlot(theseCVErrors);
formatNBP(nbp);
xticklabels({'100', '1000', '10000', '100000', '1000000'});
xlabel('\lambda');
ylabel('Normalized CV Error');

ax = gca;
ax.FontSize = fontSize;
ax.FontName = fontName;
print(fullfile(figuresFolder,'sfig_04_cverrors.tiff'),'-dtiff','-r600');

%% Fig 2d-e - reward signal after collapsing across cues (plus topo)
close all; clear all; figset;

% Settings
incAnimStart = 1;
excludeEarly = 0; % Sanity check, if needed

% Analysis electrodes
sRew = 'P4';
sRews = {'P4', 'CP2','CP6'};

% Analysis Windows
plotClusterNoBL = [-4 0.2];
plotClusterBL = [-1.2320    0.2];
rewTimeWin = plotClusterBL; % Window for the topography and effect size
animStartTimeLim = [0,0.8]; % From analysis
animEndTimeLim = [0,0.8]; % From analysis

% Load betas;
load(fullfile(resultsFolder,['deconv_collapsed_' num2str(incAnimStart) '_' num2str(useReg)  '_' num2str(excludeEarly) '.mat'])); 
load(fullfile(dataFolder,'derivatives/eegprep/sub-01/sub-01_task-gnomes_eegprep.mat'),'EEG');
allBeta = double(allBeta(whichPs,:,:));

startTimes = animStartTimeLim(1) : (1/250) : animStartTimeLim(2);
halfRewTime = 0.5*length(betaI{1}) / 250;
rewTimes = -halfRewTime:(1/250):halfRewTime-1/250;
endTimes = animEndTimeLim(1) : (1/250) : animEndTimeLim(2);

iRew = eeg_chaninds(EEG,sRew);
iRews = eeg_chaninds(EEG,sRews);

% Display Windows
startDispWin = [startTimes(1) startTimes(end)];
rewDispWin = [-4 1.5];
endDispWin = [endTimes(1) endTimes(end)];

meanBeta = squeeze(mean(allBeta,1));
stdBeta = squeeze(std(allBeta,[],1));
tval = tinv(0.025,length(whichPs)-1);
ciBeta = tval * stdBeta / sqrt(length(whichPs));

startBL = [0 0.1];
startBLPnts = dsearchn(startTimes',startBL');
if incAnimStart
    rewERP = meanBeta(betaI{1},:);
    startERP = meanBeta(betaI{2},:);
    endERP = meanBeta(betaI{3},:);
else
    startERP= [];
    rewERP = meanBeta(betaI{1},:);
    endERP = meanBeta(betaI{2},:);
end
rewCI = ciBeta(betaI{1},:);
startCI = ciBeta(betaI{2},:);
endCI = ciBeta(betaI{3},:); 

% Permutation testing
analysisWindowS = [-4 0.2];
analysisWindowI = dsearchn(rewTimes',analysisWindowS');
analysisWindowIs = analysisWindowI(1):analysisWindowI(2);

% Permutation testing - no baseline
testBeta = squeeze(allBeta(:,betaI{1},:));
testBeta = testBeta(:,analysisWindowIs,:); % Restrict analysis window
testBeta = permute(testBeta,[1 3 2]);
[clusterInfo] = doPermTest2(testBeta,0.05,neighbours,EEG);
clusterMasses = clusterInfo(:,2);
idx = cellfun('isempty',clusterMasses);
clusterMasses(idx) = {NaN};
clusterMasses = cell2mat(clusterMasses);
[maxCluster, maxClusterI] = max(clusterMasses);
disp('Permutation test - no BL');
disp(['Max cluster centered on ' clusterInfo{maxClusterI,1}]);
disp(['Neighbours: ' [neighbours(maxClusterI).neighblabel{:}]]);
clusterPoints = [clusterInfo{maxClusterI,3}]
clusterTimes = rewTimes(analysisWindowIs(clusterPoints))
clusterP = clusterInfo{maxClusterI,6}
erpScores = mean(mean(testBeta(:,iRews,clusterPoints),2),3);
erpD = mean(erpScores) / std(erpScores);
disp(['Cohen''s D: ' num2str(erpD)]);

% Permutation test - with baseline
rewBL = [-4 -1]; 
rewBLPnts = dsearchn(rewTimes',rewBL');
testBeta = squeeze(allBeta(:,betaI{1},:));
testBeta = testBeta - mean(testBeta(:,rewBLPnts(1):rewBLPnts(2),:),2); % BL correction
testBeta = testBeta(:,analysisWindowIs,:); % Restrict analysis window
testBeta = permute(testBeta,[1 3 2]);
[clusterInfoBL] = doPermTest2(testBeta,0.05,neighbours,EEG);
clusterMasses = clusterInfoBL(:,2);
idx = cellfun('isempty',clusterMasses);
clusterMasses(idx) = {NaN};
clusterMasses = cell2mat(clusterMasses);
[maxCluster, maxClusterI] = max(clusterMasses);
disp('Permutation test - with BL');
disp(['Max cluster centered on ' clusterInfoBL{maxClusterI,1}]);
disp(['Neighbours: ' [neighbours(maxClusterI).neighblabel{:}]]);
clusterPoints = [clusterInfoBL{maxClusterI,3}]
clusterTimes = rewTimes(analysisWindowIs(clusterPoints))
clusterP = clusterInfoBL{maxClusterI,6}
erpScores = mean(mean(testBeta(:,iRews,clusterPoints),2),3);
erpD = mean(erpScores) / std(erpScores);
disp(['Cohen''s D: ' num2str(erpD)]);

makefigure(6,4);
boundedline(rewTimes,mean(rewERP(:,iRews),2),rewCI(:,iRew),'alpha','cmap',plotLineColours(1,:)); hold on;
xlim(rewDispWin);
ax = gca;
ylims = ax.YLim;

% Indicate significant cluster if no baseline correction were applied
areaAlpha = 0.1;
plot(plotClusterNoBL,[0 0],'Color',[0,0,0,areaAlpha],'LineWidth',3);

% Indicate significant cluster were a -4 to -1 baseline correction applied 
areaAlpha = 0.5;
plot(plotClusterBL,[0 0],'Color',[0,0,0,areaAlpha],'LineWidth',3);

ylabel('Voltage (\muV)');
xlabel('Time (s)');
ylim(ylims);
% t = text(ax.XLim(1),ax.YLim(2),'  P4');
t = text(ax.XLim(1),ax.YLim(2),'  P4, CP2, CP6');
t.FontName = fontName;
t.FontSize = fontSize;

set(gca,'FontSize',fontSize);
set(gca,'FontName',fontName);
set(gca,'Box','off');

print(fullfile(figuresFolder,'fig_02d_collapsedrerp.tiff'),'-dtiff','-r600');

% Topo
rewTimeWinPnts = dsearchn(rewTimes',rewTimeWin');
rewTopo = mean(rewERP(rewTimeWinPnts(1):rewTimeWinPnts(2),:),1);
[minV,~] = min(rewTopo);
[maxV,~] = max(rewTopo);
topoLim = max([abs(minV) abs(maxV)]);
makefigure(5,4);
tp2 = topoplot(rewTopo,EEG.chanlocs,'maplimits' ,[-topoLim,topoLim],'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp2.Parent.XLim = [-0.6 0.6];
tp2.Parent.YLim = [-0.6 0.6];
colormap(myColormap);
c = colorbar();
c.Label.String = 'Voltage (\muV)';
c.Label.FontName = fontName;
c.Label.FontSize = fontSize;
set(gca,'FontSize',fontSize);
set(gca,'FontName',fontName);
set(gca,'Box','off');
print(fullfile(figuresFolder,'fig_02e_collapsedtopo.tiff'),'-dtiff','-r600');


%% Supplementary Fig 4 - Animation Start/Stop and Reward Signals
close all; clear all; clc; figset;

incAnimStart = 1;
useReg = 1;
excludeEarly= 0;
load(fullfile(resultsFolder,['deconv_collapsed_' num2str(incAnimStart) '_' num2str(useReg)  '_' num2str(excludeEarly) '.mat'])); 
load(fullfile(dataFolder,'derivatives/eegprep/sub-01/sub-01_task-gnomes_eegprep.mat'),'EEG');

animStartTimeLim = [0,0.8]; % From analysis
animEndTimeLim = [0,0.8]; % From analysis

startTimes = animStartTimeLim(1) : (1/250) : animStartTimeLim(2);
halfRewTime = 0.5*length(betaI{1}) / 250;
rewTimes = -halfRewTime:(1/250):halfRewTime-1/250;
endTimes = animEndTimeLim(1) : (1/250) : animEndTimeLim(2);

meanBeta = squeeze(mean(allBeta,1));
stdBeta = squeeze(std(allBeta,[],1));
tval = tinv(0.025,length(whichPs)-1);
ciBeta = tval * stdBeta / sqrt(length(whichPs));

rewCI = ciBeta(betaI{1},:);
startCI = ciBeta(betaI{2},:);
endCI = ciBeta(betaI{3},:); 

rewERP = meanBeta(betaI{1},:);
startERP = meanBeta(betaI{2},:);
endERP = meanBeta(betaI{3},:);

% Baseline correction
startBL = [0 0.1];
startBLPnts = dsearchn(startTimes',startBL');
endBL = [0 0.1];
endBLPnts = dsearchn(endTimes',endBL');
startERP = startERP - mean(startERP(startBLPnts(1):startBLPnts(2),:),1);
endERP = endERP - mean(endERP(endBLPnts(1):endBLPnts(2),:),1);

% Display Windows
startDispWin = [startTimes(1) startTimes(end)];
rewDispWin = [-4 1.5];
endDispWin = [endTimes(1) endTimes(end)];

subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.2 0.1], [0.01 0.12], [0.08 0.04]);

sStart = 'O2';
sRew = 'P4';
sEnd = 'Cz';

iStart = eeg_chaninds(EEG,sStart);
iRew = eeg_chaninds(EEG,sRew);
iEnd = eeg_chaninds(EEG,sEnd);

makefigure(18,10);

eegAxs = {};
topoAxs = {};

% Find peaks based on min/max
[~,startMinI] = min(startERP(:,iStart))
[~,rewMinI] = min(rewERP(:,iRew))
[~,endMaxI] = max(endERP(:,iEnd))

startMinTime = startTimes(startMinI);
rewMinTime = 0;
endMaxTime = endTimes(endMaxI);

eegAxs{1} = subplot(2,3,1);
boundedline(startTimes,startERP(:,iStart),startCI(:,iStart),'alpha','cmap',plotLineColours(1,:)); xlim(startDispWin); hold on;
xline(startMinTime(1),'Color',[0.75 0.75 0.75]); title({'Start of Animation',''},'FontWeight','normal');
ax = gca; text(ax.XLim(1),ax.YLim(2),['  ' sStart],'FontSize',fontSize,'FontWeight','normal','FontName',fontName);
eegAxs{2} = subplot(2,3,2);
boundedline(rewTimes,rewERP(:,iRew),rewCI(:,iRew),'alpha','cmap',plotLineColours(1,:)); xlim(rewDispWin); hold on;
xline(rewMinTime(1),'Color',[0.75 0.75 0.75]); title({'During Animation',''},'FontWeight','normal');
ax = gca; text(ax.XLim(1),ax.YLim(2),['  ' sRew],'FontSize',fontSize,'FontWeight','normal','FontName',fontName);
eegAxs{3} = subplot(2,3,3);
boundedline(endTimes,endERP(:,iEnd),endCI(:,iEnd),'alpha','cmap',plotLineColours(1,:)); xlim(endDispWin); hold on;
xline(endMaxTime(1),'Color',[0.75 0.75 0.75]); title({'End of Animation',''},'FontWeight','normal');
ax = gca; text(ax.XLim(1),ax.YLim(2),['  ' sEnd],'FontSize',fontSize,'FontWeight','normal','FontName',fontName);

% Topos
startTimeWin = [0.24 0.24];
startTimeWinPnts = dsearchn(startTimes',startTimeWin');
startTopo = mean(startERP(startTimeWinPnts(1):startTimeWinPnts(2),:),1);
[minV,minI] = min(startTopo);
[maxV,maxI] = max(startTopo);
disp(['min: ' EEG.chanlocs(minI).labels]);
disp(['max: ' EEG.chanlocs(maxI).labels]);
topoLim = max([abs(minV) abs(maxV)]);
subplot(2,3,4);
topoAxs{1} = topoplot(startTopo,EEG.chanlocs,'maplimits' ,[-topoLim,topoLim],'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');

rewTimeWin = [0.014 0.014];
rewTimeWin = [ -0.6760    0.2000];
rewTimeWinPnts = dsearchn(rewTimes',rewTimeWin');
rewTopo = mean(rewERP(rewTimeWinPnts(1):rewTimeWinPnts(2),:),1);
[minV,minI] = min(rewTopo);
[maxV,maxI] = max(rewTopo);
disp(['min: ' EEG.chanlocs(minI).labels]);
disp(['max: ' EEG.chanlocs(maxI).labels]);
topoLim = max([abs(minV) abs(maxV)]);
subplot(2,3,5);
topoAxs{2} = topoplot(rewTopo,EEG.chanlocs,'maplimits' ,[-topoLim,topoLim],'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');

endTimeWin = [0.432 0.432];
endTimeWinPnts = dsearchn(endTimes',endTimeWin');
endTopo = mean(endERP(endTimeWinPnts(1):endTimeWinPnts(2),:),1);
[minV,minI] = min(endTopo);
[maxV,maxI] = max(endTopo);
disp(['min: ' EEG.chanlocs(minI).labels]);
disp(['max: ' EEG.chanlocs(maxI).labels]);
topoLim = max([abs(minV) abs(maxV)]);
subplot(2,3,6);
topoAxs{3} = topoplot(endTopo,EEG.chanlocs,'maplimits' ,[-topoLim,topoLim],'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');

for i = 1:3
    eegAxs{i}.XLabel.String = 'Time (s)';
    eegAxs{i}.YLabel.String = 'Voltage (\muV)';
    topoAxs{i}.Parent.XLim = [-0.6 0.6];
    topoAxs{i}.Parent.YLim = [-0.6 0.6];
    topoAxs{i}.Parent.Colormap = myColormap;
end

print(fullfile(figuresFolder,'sfig_05_allcollapsed.tiff'),'-dtiff','-r600');

%% Supplementary Figure 3 - Total Regressor Count
close all; clear all; figset;

% Load and plot a sample design matrix and some other info
sampleP = 'sub-01';
sub1data = load([projectFolder '/data/derivatives/glmres/' sampleP '/' sampleP '_task-gnomes_glmcollapsed_1_0_0.mat']);
halfRewTime = 0.5*(sub1data.betaI{1}(end)-1) / 250;
rewTimes = -halfRewTime:(1/250):halfRewTime;

% Combined design matrices and count each regressp (sum since they are 1's)
allX = {};
for p = 1:length(whichPs)
    thisP = ['sub-' num2str(whichPs(p),'%0.2i')];
    load([projectFolder '/data/derivatives/glmres/' thisP '/' thisP '_task-gnomes_glmcollapsed_1_0_0.mat'],'X');
    allX{p} = X;
end

makefigure(8,5);
whichIs = [1];
axs = [];
for i = 1:length(whichIs)
    allSum = [];
    for p = 1:length(allX)
        allSum(p,:) = full(sum(allX{p}(:,sub1data.betaI{whichIs(i)}),1));
    end

    meanSum = mean(allSum,1);
    sdSum = std(allSum,[],1);
    tval = abs(tinv(0.025,size(allSum,1)));
    ciSum = tval * sdSum ./ sqrt(size(allSum,1));

    % plot(rewTimes,mean(allSum,1),'LineWidth',lineWidth);
    [hl, hp] = boundedline(rewTimes,meanSum,ciSum,'cmap',plotLineColours(i,:),'alpha'); 
    axs(i) = hl;
    ylabel('Total samples');
    ax = gca;
    % ax.XAxis.TickLabels = [];
    %ax.XAxis.TickValues = [sub1data.betaI{1}(end)];
    set(gca,'FontSize',fontSize);
    set(gca,'Box','off');

    hold on;
end
xlabel('Time (s)');

print(fullfile(figuresFolder,'sfig_03_regcountcollapsed.tiff'),'-dtiff','-r600');

%% UNUSED - Regressor count by condition
close all; clear all; figset;

% Load and plot a sample design matrix and some other info
sampleP = 'sub-01';
sub1data = load([projectFolder '/data/derivatives/glmres/' sampleP '/' sampleP '_task-gnomes_glm_1_0.mat']);
halfRewTime = 0.5*(sub1data.betaI{1}(end)-1) / 250;
rewTimes = -halfRewTime:(1/250):halfRewTime;

% xs = 1:(sub1data.betaI{1}(end)-1);
% cspy(X(:,xs));

% Combined design matrices and count each regressp (sum since they are 1's)
ps = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
allX = {};
for p = 1:length(whichPs)
    thisP = ['sub-' ps{whichPs(p)}];
    load([projectFolder '/data/derivatives/glmres/' thisP '/' thisP '_task-gnomes_glm_1_0.mat'],'X');
    allX{p} = X;
end

makefigure(14,5);
whichIs = [1 4 7];
axs = [];
for i = 1:length(whichIs)
    allSum = [];
    for p = 1:length(allX)
        allSum(p,:) = full(sum(allX{p}(:,sub1data.betaI{whichIs(i)}),1));
    end

    meanSum = mean(allSum,1);
    sdSum = std(allSum,[],1);
    tval = abs(tinv(0.025,size(allSum,1)));
    ciSum = tval * sdSum ./ sqrt(size(allSum,1));

    % plot(rewTimes,mean(allSum,1),'LineWidth',lineWidth);
    [hl, hp] = boundedline(rewTimes,meanSum,ciSum,'cmap',plotLineColours(i,:),'alpha'); 
    axs(i) = hl;
    ylabel('Total samples');
    ax = gca;
    % ax.XAxis.TickLabels = [];
    %ax.XAxis.TickValues = [sub1data.betaI{1}(end)];
    set(gca,'FontSize',fontSize);
    set(gca,'Box','off');

    hold on;
end
xlabel('Time (s)');
legend(axs,{'high','medium','low'},'Box','off');


%% Fig 3 - Reward Signals by Condition
close all; clear all; clc; figset;

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

% Load results
incAnimStart = 1;
useReg = 1;
load(fullfile(resultsFolder,['deconv_' num2str(incAnimStart) '_' num2str(useReg)  '.mat'])); 

% Analysis electrodes
iRew = eeg_chaninds(EEG,'Pz'); % P4
iRews = eeg_chaninds(EEG,{'Pz', 'CP1', 'CPz', 'CP2'});

% Analysis Windows
rewTimeWin = [-0.5160    -0.164];
halfRewTime = 0.5*length(betaI{1}) / 250;
rewTimes = -halfRewTime:(1/250):(halfRewTime - 1/250);

% Reshape beta: condition (high, med, low) X participant X time
clear rBeta;

    rBeta(1,:,:,:)  = allBeta(whichPs,[betaI{1} betaI{2} betaI{3}],:);
    rBeta(2,:,:,:)  = allBeta(whichPs,[betaI{4} betaI{5} betaI{6}],:);
    rBeta(3,:,:,:)  = allBeta(whichPs,[betaI{7} betaI{8} betaI{9}],:);

rewBeta = rBeta(:,:,betaI{1},:); % just reward part
allDiff = squeeze(rBeta(3,:,:,:) - rBeta(1,:,:,:));

meanBeta = squeeze(mean(rBeta,2));

rewERP = meanBeta(:,betaI{1},:);

% Reward Plot
makefigure(8,7);
rewDispWin = [-2 1.5];
sps{1} = subplot(2,4,1:4);
% plot(rewTimes,rewERP(:,:,iRew), 'LineWidth',lineWidth); box off; hold on;
plot(rewTimes,mean(rewERP(:,:,iRews),3), 'LineWidth',lineWidth); box off; hold on;
ax = gca;
ax.ColorOrder = plotLineColours;
plot(rewTimeWin,[0 0],'Color',[0,0,0,0.5],'LineWidth',3);
xlim(rewDispWin);
ylims = sps{1}.YLim;
%area(sps{1},rewTimeWin, [sps{1}.YLim(1) sps{1}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
%area(sps{1},rewTimeWin, [sps{1}.YLim(2) sps{1}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');
l = legend(sps{1},'high', 'med','low','Box','off','Location','SouthEast');
l.Position = [0.6902    0.5949    0.2478    0.1641];
ylim(ylims);
ax = gca;
set(ax,'FontSize',fontSize);
set(ax,'FontName',fontName);
set(ax,'Box','off');
t = text(ax.XLim(1),ax.YLim(2),' Pz, CP1, CPz, CP2');
t.FontName = fontName;
t.FontSize = fontSize;

rewTimeWinPnts = dsearchn(rewTimes',rewTimeWin');
meanRewERP = squeeze(mean(rewERP(:,rewTimeWinPnts(1):rewTimeWinPnts(2),:),2));
maxTopo = max(max(abs(meanRewERP)));
topoLimits = [-maxTopo maxTopo];

titles = {'high','medium','low'};
for i = 1:3
    subplot(2,4,4+i);
    tp2 = topoplot(meanRewERP(i,:),EEG.chanlocs, 'maplimits',topoLimits,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
    tp2.Parent.XLim = [-0.6 0.6];
    tp2.Parent.YLim = [-0.6 0.6];
    colormap(myColormap);
    set(gca,'FontSize',fontSize);
    set(gca,'FontName',fontName);
    set(gca,'Box','off');

    t = title(titles{i});
    t.FontSize = fontSize;
    t.FontName = fontName;
    t.FontWeight = 'normal';
    t.Position(2) = -0.9;
end

% High-low difference topography (revision 1)
subplot(2,4,8);
highMinusLowDiff =  meanRewERP(1,:) - meanRewERP(3,:);
tp3 = topoplot(highMinusLowDiff,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp3.Parent.XLim = [-0.6 0.6];
tp3.Parent.YLim = [-0.6 0.6];
colormap(myColormap);
set(gca,'FontSize',fontSize);
set(gca,'FontName',fontName);
set(gca,'Box','off');
t = title('high-low');
t.FontSize = fontSize;
t.FontName = fontName;
t.FontWeight = 'normal';
t.Position(2) = -0.9;

print(fullfile(figuresFolder,'fig_03_rerpsbycondition.tiff'),'-dtiff','-r600');
data = permute(rewBeta,[2,1,4,3]);
testData = squeeze(data(:,1,:,:) - data(:,3,:,:));
[clusterInfo] = doPermTest2(testData,0.05,neighbours,EEG);

clusterMasses = clusterInfo(:,2);
idx = cellfun('isempty',clusterMasses);
clusterMasses(idx) = {NaN};
clusterMasses = cell2mat(clusterMasses);
[maxCluster, maxClusterI] = max(clusterMasses);
disp(['Max cluster centered on ' clusterInfo{maxClusterI,1}]);
disp(['Neighbours: ' [neighbours(maxClusterI).neighblabel{:}]]);

clusterPoints = [clusterInfo{maxClusterI,3}]
clusterTimes = rewTimes(clusterPoints)
clusterP = clusterInfo{maxClusterI,6}

% Compute effect size
erpScores = mean(mean(allDiff(:,rewTimeWinPnts(1):rewTimeWinPnts(2),iRews),2),3);
erpD = mean(erpScores) / std(erpScores);

disp(['Cohen''s d: ' num2str(erpD)]);

%% Fig 4 - effect of expectancy 

close all; clear all; clc; figset;
load(fullfile(resultsFolder,'deconv_collapsedaverew_1_1.mat'));

% Rescale "expectancy" beta so that points go 1-100 not 0.01-1
allBeta(:,betaI{2},:) = allBeta(:,betaI{2},:)/100; 
allBeta = allBeta(whichPs,:,:); % Restrict to included Ps only

% Summary stats
meanBeta = squeeze(mean(allBeta,1));
stdBeta = squeeze(std(allBeta,[],1));
tval = tinv(0.025,length(whichPs)-1);
ciBeta = tval * stdBeta / sqrt(length(whichPs));

% * * * Settings * * *

% Do baseline correction for the average reward signal?
% doRewBL = 1;

% Electrodes of intrerest (for plots and for computing effect size)
sRew = 'P4';
sRews = {'P4','CP2','CP6'};
sAveRew = 'F4';
sAveRews = {'F4','FC2','FC6'};
sEnd = 'Cz';
iRew = eeg_chaninds(EEG,sRew);
iRews = eeg_chaninds(EEG,sRews);
iExp = eeg_chaninds(EEG,sAveRew);
iExps = eeg_chaninds(EEG,sAveRews);
iEnd = eeg_chaninds(EEG,sEnd);

% Compute cluster ERPs (average across cluster)
rewBeta = allBeta(:,betaI{1},:);
meanRewBeta = squeeze(mean(mean(rewBeta(:,:,iRews),3),1));
stdRewBeta = squeeze(std(mean(rewBeta(:,:,iRews),3),[],1));
ciRewBeta = tval * stdRewBeta / sqrt(length(whichPs));
expBeta = allBeta(:,betaI{2},:);
meanExpBeta = squeeze(mean(mean(expBeta(:,:,iExps),3),1));
stdExpBeta = squeeze(std(mean(expBeta(:,:,iExps),3),[],1));
ciExpBeta = tval * stdExpBeta / sqrt(length(whichPs));


% Analysis windows (for plots and for computing effect size)
barTimeWin = [0.220 0.275];
rewTimeWin = [-0.468    0.192];
rewTimeWinNoBL = [-4 0.2];
aveRewTimeWin = [-1.6920   -0.5440];

% Plot display windows
rewDispWin = [-3 1.5];
aveRewDispWin = rewDispWin;
endTimeWin = [0.425 0.550];

% Generate time/point vectors
halfRewTime = 0.5*length(betaI{1}) / 250;
rewTimes = -halfRewTime:(1/250):halfRewTime - 1/250;
expTimes = rewTimes;
animEndTimeLim = [0,0.8]; % From the analysis
endTimes = animEndTimeLim(1) : (1/250) : animEndTimeLim(2);
rewTimeWinPnts = dsearchn(rewTimes',rewTimeWin');
rewTimeWinNoBLPnts = dsearchn(rewTimes',rewTimeWinNoBL');
aveRewTimeWinPnts = dsearchn(expTimes',aveRewTimeWin');
endTimeWinPnts = dsearchn(endTimes',endTimeWin');

% Compute ERPs 
rewERP = meanBeta(betaI{1},:);
aveRewERP = meanBeta(betaI{2},:);
endERP = meanBeta(betaI{3},:);

%  Compute topos based on windows defined above
rewTopo = mean(rewERP(rewTimeWinPnts(1):rewTimeWinPnts(2),:),1);
aveRewTopo = mean(aveRewERP(aveRewTimeWinPnts(1):aveRewTimeWinPnts(2),:),1);
animEndTopo = mean(endERP(endTimeWinPnts(1):endTimeWinPnts(2),:),1);

% Average Effect NO baseline correctionth 
aveRewBeta = squeeze(allBeta(:,betaI{1},:));
aveRewBeta = permute(aveRewBeta,    [1 3 2]);
analysisWindowS = [-4 0.2];
analysisWindowI = dsearchn(rewTimes',analysisWindowS');
analysisWindowIs = analysisWindowI(1):analysisWindowI(2);
[averRewClusterInfo] = doPermTest2(aveRewBeta(:,:,analysisWindowIs),0.05,neighbours,EEG);
clusterMasses = averRewClusterInfo(:,2);
idx = cellfun('isempty',clusterMasses);
clusterMasses(idx) = {NaN};
clusterMasses = cell2mat(clusterMasses);
[maxCluster, maxClusterI] = max(clusterMasses);
disp('Average Effect - No Baseline');
disp(['Max cluster centered on ' averRewClusterInfo{maxClusterI,1}]);
disp(['Neighbours: ' [neighbours(maxClusterI).neighblabel{:}]]);
clusterPoints = [averRewClusterInfo{maxClusterI,3}]
clusterTimes = rewTimes(analysisWindowIs(clusterPoints))
clusterP = averRewClusterInfo{maxClusterI,6}
rewERPScoresNoBL = mean(mean(aveRewBeta(:,iRews,rewTimeWinNoBLPnts(1):rewTimeWinNoBLPnts(2)),2),3);
rewNoBLD = mean(rewERPScoresNoBL) / std(rewERPScoresNoBL);
disp(['Cohen''s D: ' num2str(rewNoBLD)]);

% Average Effect WITH baseline correction
rewBL = [-4 -1]; % In line with earlier figure
rewBLPnts = dsearchn(rewTimes',rewBL');
aveRewBetaBL = squeeze(allBeta(:,betaI{1},:));
aveRewBetaBL = aveRewBetaBL - mean(aveRewBetaBL(:,rewBLPnts(1):rewBLPnts(2),:),2);
aveRewBetaBL = permute(aveRewBetaBL,    [1 3 2]);
analysisWindowS = [-4 0.2];
analysisWindowI = dsearchn(rewTimes',analysisWindowS');
analysisWindowIs = analysisWindowI(1):analysisWindowI(2);
[averRewClusterInfo] = doPermTest2(aveRewBetaBL(:,:,analysisWindowIs),0.05,neighbours,EEG);
clusterMasses = averRewClusterInfo(:,2);
idx = cellfun('isempty',clusterMasses);
clusterMasses(idx) = {NaN};
clusterMasses = cell2mat(clusterMasses);
[maxCluster, maxClusterI] = max(clusterMasses);
disp('Average Effect');
disp(['Max cluster centered on ' averRewClusterInfo{maxClusterI,1}]);
disp(['Neighbours: ' [neighbours(maxClusterI).neighblabel{:}]]);
clusterPoints = [averRewClusterInfo{maxClusterI,3}]
clusterTimes = rewTimes(analysisWindowIs(clusterPoints))
clusterP = averRewClusterInfo{maxClusterI,6}
rewERPScores = mean(mean(aveRewBetaBL(:,iRews,rewTimeWinPnts(1):rewTimeWinPnts(2)),2),3);
rewD = mean(rewERPScores) / std(rewERPScores);
disp(['Cohen''s D: ' num2str(rewD)]);

% Expectancy Effect
expBeta = squeeze(allBeta(:,betaI{2},:));
expBeta = permute(expBeta,[1 3 2]);
[clusterInfo] = doPermTest2(expBeta,0.05,neighbours,EEG);
clusterMasses = clusterInfo(:,2);
idx = cellfun('isempty',clusterMasses);
clusterMasses(idx) = {NaN};
clusterMasses = cell2mat(clusterMasses);
[maxCluster, maxClusterI] = max(clusterMasses);
disp('Expectancy Effect');
disp(['Max cluster centered on ' clusterInfo{maxClusterI,1}]);
disp(['Neighbours: ' [neighbours(maxClusterI).neighblabel{:}]]);
clusterPoints = [clusterInfo{maxClusterI,3}]
clusterTimes = rewTimes(clusterPoints)
clusterP = clusterInfo{maxClusterI,6}
expERPScores = mean(mean(expBeta(:,iExps,aveRewTimeWinPnts(1):aveRewTimeWinPnts(2)),2),3);
expD = mean(expERPScores) / std(expERPScores);
disp(['Cohen''s D: ' num2str(expD)]);


% Plots
makefigure(8,6);
% boundedline(rewTimes,meanBeta(betaI{1},iRew),ciBeta(betaI{1},iRew),'alpha','cmap',plotLineColours(1,:)); hold on;\
boundedline(rewTimes,meanRewBeta,ciRewBeta,'alpha','cmap',plotLineColours(1,:)); hold on;
plot(rewTimeWinNoBL,[0 0],'Color',[0,0,0,0.1],'LineWidth',3)
plot(rewTimeWin,[0 0],'Color',[0,0,0,0.5],'LineWidth',3);
xlim(rewDispWin);
ax = gca;
xlabel('Time (s)'); ylabel('Voltage (\muV)'); 
ax.FontName = fontName;
ax.FontSize = fontSize;
% t = text(ax.XLim(1),ax.YLim(2),['  ' sRew]);
t = text(ax.XLim(1),ax.YLim(2),['  ' sRews{1} ', ' sRews{2} ', ' sRews{3}]);
t.FontName = fontName;
t.FontSize = fontSize;
print(fullfile(figuresFolder,'fig_04NW_aveeffect.tiff'),'-dtiff','-r600');

makefigure(8,6);
% boundedline(expTimes,meanBeta(betaI{2},iExp),ciBeta(betaI{2},iExp),'alpha','cmap',plotLineColours(1,:)); hold on;
boundedline(expTimes,meanExpBeta,ciExpBeta,'alpha','cmap',plotLineColours(1,:)); hold on;
plot(clusterTimes,[0 0],'Color',[0,0,0,0.5],'LineWidth',3);
xlim(aveRewDispWin);
ax = gca;
xlabel('Time (s)'); ylabel('Beta (\muV/point)'); 
ax.FontName = fontName;
ax.FontSize = fontSize;
% t = text(ax.XLim(1),ax.YLim(2),['  ' sAveRew]);
t = text(ax.XLim(1),ax.YLim(2),['  ' sAveRews{1} ', ' sAveRews{2} ', ' sAveRews{3} ]);
t.FontName = fontName;
t.FontSize = fontSize;
print(fullfile(figuresFolder,'fig_04NE_expeffect.tiff'),'-dtiff','-r600');

makefigure(8,6);
tp2 = topoplot(rewTopo,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp2.Parent.XLim = [-0.6 0.6];
tp2.Parent.YLim = [-0.6 0.6];
colormap(myColormap);
c = colorbar();
c.Label.String = 'Voltage (\muV)';
c.FontName = fontName;
c.FontSize = fontSize;
print(fullfile(figuresFolder,'fig_04SW_aveeffecttopo.tiff'),'-dtiff','-r600');

makefigure(8,6);
tp2 = topoplot(aveRewTopo,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp2.Parent.XLim = [-0.6 0.6];
tp2.Parent.YLim = [-0.6 0.6];
colormap(myColormap);
c = colorbar();
c.Label.String = 'Beta (\muV/point)';
c.FontName = fontName;
c.FontSize = fontSize;
print(fullfile(figuresFolder,'fig_04SE_expeffecttopo.tiff'),'-dtiff','-r600');


%%
figure();
topoplot([],readlocs('ccnlabactichamp.locs'),'style','blank','electrodes','labelpoint','headrad','rim','whitebk','on');
