function makeFigurePanels
% makeFigurePanels uses the simulated data to reproduce the panels used in
% the paper
% 
% Chen Chen

warning('off', 'MATLAB:MKDIR:DirectoryExists');
warning('off', 'MATLAB:rmpath:DirNotFound')
%% Specify Target Figure to Plot (Change this as needed)
% Target figure panel
%    Choose one of the following:
%       'fig1'   -  panels for figure 1
%       'fig2'   -  panels for figure 2
%       'fig3'   -  panels for figure 3, note that figure 2 simulation is
%                   required to use USE_PUBLISHED_DATASET = 0, and panel C 
%                   relys on the simulated data from sm-fig1, please run 
%                   sm-fig1 figure code to get fig3C reproduced. 
%                   (It will be copied over to fig3 folder automatically.)
%      'sm-fig1' -  panels for supplement figure 1.
%      'sm-fig2' -  panels for supplement figure 2, fig2 simulation data 
%                   are required to use USE_PUBLISHED_DATASET = 0
%      'sm-fig3' -  panels for supplement figure 3, note that the EIH 
%                   simulation for figure 2 is required to use
%                   USE_PUBLISHED_DATASET = 0
%      'sm-fig4' -  panels for supplement figure 4.  Simulation will be 
%                   submitted if no available locally simulated data were 
%                   found when using USE_PUBLISHED_DATASET = 0, or forced
%                   running new simulation using USE_PUBLISHED_DATASET = 2
%                   The simulation consists a total of 80 trials and will
%                   normally take around 12 hours to finish.
%      'sm-fig5' -  panels for supplement figure 5
%      'sm-fig6' -  panels for supplement figure 6. Note that this figure 
%                   does not require any simulation data and therefore
%                   USE_PUBLISHED_DATASET will be ignored
% 
targetFig = 'fig6';

% Control whether or not to use previously simulated dataset
% Use flag (USE_PUBLISHED_DATASET = flag)
%   2 | (only available for sm-fig4) force sm-fig4 to submit new simulation
%   1 | use previouly published dataset (default)
%   0 | use locally simulated data if possible, otherwise proceed with new
%       simulation (sm-fig4)
USE_PUBLISHED_DATASET = 0;

%% Internal parameters (do not change)
if USE_PUBLISHED_DATASET == 1
    FIG_DATA_PATH = './PublishedData/';
else
    FIG_DATA_PATH = '../SimulationCode/SimData/';
end
FIG_OUTPUT_PATH = sprintf('./FigureOutput/%s/', targetFig);

% add figure code and common path
addpath('./FigureCode');
addpath('./FigureCode/common')
addpath('./FigureCode/common/boundedline/boundedline/');
addpath('./FigureCode/common/boundedline/Inpaint_nans/');

%% Make production figure panels (do not change)
mkdir(FIG_OUTPUT_PATH);
switch targetFig
    case 'fig1'
        makeFig1Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig2'
        makeFig2Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig3'
        makeFig3Plots(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig4'
        makeFig4Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig5'
        makeFig5Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig6'
        if ~USE_PUBLISHED_DATASET
            Fig6ProcessData(FIG_DATA_PATH, FIG_OUTPUT_PATH);
        end
        makeFig6Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig7'
        makeFig7Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'sm-fig1'
        makeFigS1Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'sm-fig2'
        makeFigS2Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    otherwise
        error('Target figure not found!');
end

function cpySimDataFiles(srcPattern, dstDir)
fileList = dir(srcPattern);
if isempty(fileList)
    error('No matching data found, have you chosen the right parameters?');
end
for i = 1:length(fileList)
    copyfile(fullfile(fileList(i).folder,fileList(i).name), ...
        fullfile(dstDir, fileList(i).name));
end