function makeFigurePanels
%% Main function to reproduce published figures and their supplements
% for "Tuning movement for sensing in an uncertain world"

close all;
warning('off', 'MATLAB:MKDIR:DirectoryExists');
warning('off', 'MATLAB:rmpath:DirNotFound')

%% Specify Target Figure to Plot (Change this as needed)
% Target figure panel
%    Choose one of the following:
%       'all'    -  reproduce all figures
%       'fig1'   -  panels for figure 1
%       'fig2'   -  panels for figure 2
%       'fig2s1' -  panels for figure 2---figure supplement 1
%       'fig2s2' -  panels for figure 2---figure supplement 2
%       'fig2s3' -  panels for figure 2---figure supplement 3
%       'fig3'   -  panels for figure 3
%       'fig4'   -  panels for figure 4
%       'fig4s1' -  panels for figure 4---figure supplement 1
%       'fig5'   -  panels for figure 5
%       'fig6'   -  panels for figure 6
%       'fig6s1' -  panels for figure 6---figure supplement 1
%       'fig6s2' -  panels for figure 6---figure supplement 2
%       'fig6s3' -  panels for figure 6---figure supplement 3
%       'fig6s4' -  panels for figure 6---figure supplement 4
%       'fig6s5' -  panels for figure 6---figure supplement 5
%       'fig7'   -  panels for figure 7
targetFig = 'all';

% Control whether or not to use previously simulated dataset
%   1 | use previouly published dataset (default)
%   0 | use locally simulated data
USE_PUBLISHED_DATASET = 1;

%% Internal parameters (do not change)
FIG_DATA_PATH = './PublishedData/';
% locally simulated data, only available when USE_PUBLISHED_DATASET == 1
FIG_DATA_PATH_LOCAL = '../SimulationCode/SimData/';
FIG_OUTPUT_PATH = sprintf('./FigureOutput/%s/', targetFig);
NeedSimData = @(s) isempty(find(ismember({'fig2s3', 'fig4s1'}, s), 1));
if ~USE_PUBLISHED_DATASET && NeedSimData(targetFig)
    FIG_DATA_PATH = FIG_DATA_PATH_LOCAL;
end

% add figure code and common path
addpath('./FigureCode');
addpath('./FigureCode/common')
addpath('./FigureCode/common/boundedline/boundedline/');
addpath('./FigureCode/common/boundedline/Inpaint_nans/');

%% Make production figure panels (do not change)
allFigs = {...
    'fig1',...
    'fig2', 'fig2s1', 'fig2s2', 'fig2s3', ...
    'fig3', ...
    'fig4', 'fig4s1', ...
    'fig5', ...
    'fig6', 'fig6s1', 'fig6s2', 'fig6s3', 'fig6s4', 'fig6s5', ...
    'fig7'};
if strcmp(targetFig, 'all')
    for i = 1:length(allFigs)
        targetFig = allFigs{i};
        FIG_OUTPUT_PATH = sprintf('./FigureOutput/%s/', targetFig);
        if ~USE_PUBLISHED_DATASET && NeedSimData(targetFig)
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
rng(0);
switch targetFig
    case 'fig1'
        makeFigure1Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig2'
        makeFigure2Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig2s1'
        makeFigure2S1Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig2s2'
        makeFigure2S2Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig2s3'
        makeFigure2S3Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig3'
        makeFigure3Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig4'
        if ~USE_PUBLISHED_DATASET
            Figure45ProcessData(FIG_DATA_PATH, FIG_OUTPUT_PATH);
        end
        makeFigure4Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig4s1'
        makeFigure4S1Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig5'
        if ~USE_PUBLISHED_DATASET
            Figure45ProcessData(FIG_DATA_PATH, FIG_OUTPUT_PATH);
        end
        makeFigure5Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig6'
        makeFigure6Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig6s1'
        makeFigure6S1Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH, USE_PUBLISHED_DATASET);
    case 'fig6s2'
        makeFigure6S2Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig6s3'
        makeFigure6S3Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig6s4'
        makeFigure6S4Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig6s5'
        makeFigure6S5Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig7'
        makeFigure7Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    otherwise
        error('Target figure not found!');
end
