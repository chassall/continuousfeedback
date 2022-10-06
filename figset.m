% Common settings for stats and figures
%
% Other m-files required: 
% brewermap: https://github.com/DrosteEffect/BrewerMap

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Participants to include
whichPs = [1:10 12:21]; % P11 noisy

% Set project folder
projectFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_Gnomes_Hassall';
dataFolder = [projectFolder '/data'];
resultsFolder = [projectFolder '/analysis/results'];
figuresFolder = [projectFolder '/figures']; 

% Plot settings
lineWidth = 1;
fontSize = 8;
fontName = 'Arial';
allColours = brewermap(9,'Set1');
plotLineColours = allColours([2 5 4 3],:);
myColormap = brewermap(256,'RdBu');
myColormap = flip(myColormap);

% Load and modify neighbours file (defines electrode adjacency)
load(fullfile(dataFolder,'derivatives/eegprep/sub-01/sub-01_task-gnomes_eegprep.mat'),'EEG');
load('elec1010_neighb.mat');
ourChannels = {EEG.chanlocs.labels};
toDelete = [];
for i = 1:length(neighbours)
    invNeighb = ~ismember(neighbours(i).neighblabel,ourChannels);
    neighbours(i).neighblabel(invNeighb) = [];
    thisLabel = neighbours(i).label;
    if ~ismember(thisLabel,ourChannels)
        toDelete = [toDelete i];
    end
end
neighbours(toDelete) = []; % Remove invalid channels