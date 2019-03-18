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
%       'all'    -  reproduce all figures
%       'fig2'   -  panels for figure 2
%       'fig3'   -  panels for figure 3
%       'fig4'   -  panels for figure 4
%       'fig5'   -  panels for figure 5
%       'fig6'   -  panels for figure 6
%       'fig7'   -  panels for figure 7
%       'sm-fig1'-  panels for figure S1
%       'sm-fig2'-  panels for figure S2
%       'sm-fig3'-  panels for figure S3
%       'sm-fig4'-  panels for figure S4
%       'sm-fig5'-  panels for figure S5
%       'sm-fig6'-  panels for figure S6
%
targetFig = 'sm-fig6';

% Control whether or not to use previously simulated dataset
% Use flag (USE_PUBLISHED_DATASET = flag)
%   1 | use previouly published dataset (default)
%   0 | use locally simulated data
USE_PUBLISHED_DATASET = 1;

%% Internal parameters (do not change)
FIG_DATA_PATH = './PublishedData/';
% locally simulated data, only available when USE_PUBLISHED_DATASET == 1
FIG_DATA_PATH_LOCAL = '../SimulationCode/SimData/';
FIG_OUTPUT_PATH = sprintf('./FigureOutput/%s/', targetFig);

flagNeedSimData = @(s) isempty(find(ismember({'sm-fig4', 'sm-fig6'}, s), 1));
if ~USE_PUBLISHED_DATASET && flagNeedSimData(targetFig)
    FIG_DATA_PATH = FIG_DATA_PATH_LOCAL;
end

% add figure code and common path
addpath('./FigureCode');
addpath('./FigureCode/common')
addpath('./FigureCode/common/boundedline/boundedline/');
addpath('./FigureCode/common/boundedline/Inpaint_nans/');

%% Make production figure panels (do not change)
if strcmp(targetFig, 'all')
    allFigs = {'fig2', 'fig3', 'fig4', 'fig5', 'fig6', 'fig7', ...
        'sm-fig1', 'sm-fig2', 'sm-fig3', 'sm-fig4', 'sm-fig5', 'sm-fig6'};
    for i = 1:length(allFigs)
        targetFig = allFigs{i};
        FIG_OUTPUT_PATH = sprintf('./FigureOutput/%s/', targetFig);
        if ~USE_PUBLISHED_DATASET && flagNeedSimData(targetFig)
            FIG_DATA_PATH = FIG_DATA_PATH_LOCAL;
        else
            FIG_DATA_PATH = './PublishedData/';
        end
        callFigureCode(targetFig, FIG_DATA_PATH, FIG_OUTPUT_PATH, USE_PUBLISHED_DATASET);
    end
else
    callFigureCode(targetFig, FIG_DATA_PATH, FIG_OUTPUT_PATH, USE_PUBLISHED_DATASET)
end

function callFigureCode(targetFig, FIG_DATA_PATH, FIG_OUTPUT_PATH, USE_PUBLISHED_DATASET)
mkdir(FIG_OUTPUT_PATH);
switch targetFig
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
    case 'sm-fig3'
        makeFigS3Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH, USE_PUBLISHED_DATASET);
    case 'sm-fig4'
        makeFigS4Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'sm-fig5'
        makeFigS5Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'sm-fig6'
        makeFigS6Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    otherwise
        error('Target figure not found!');
end
