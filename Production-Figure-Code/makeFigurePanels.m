function makeFigurePanels
% makeFigurePanels uses the simulated data to reproduce the panels used in
% the paper
% 
% Chen Chen

%% Specify Target Figure to Plot
% Target figure panel
%    Choose one of the following:
%       'fig1'   -  panels for figure 1
%       'fig2'   -  panels for figure 2
% 
%      'sm-fig4' - this is slightly more complicated than the rest.
%                  In order to reproduce supplement figure 4, 
targetFig = 'fig1';

% Choose whether or not to use previously simulated dataset
% Use
%   USE_PREV_DATASET = 1; % if local simulation step is skipped
%   USE_PREV_DATASET = 0; % if local simulation is done
USE_PREV_DATASET = 1;

%% Internal parameters (do not change)
if USE_PREV_DATASET
    DATA_PATH = './FigureCode/';
    FIG_DATA_PATH = sprintf([DATA_PATH,'%s/Data/'], targetFig);
else
    DATA_PATH = '../SimulationCode/SimData/';
    FIG_DATA_PATH = sprintf([DATA_PATH,'%s/'], targetFig);
end
FIG_OUTPUT_PATH = sprintf('./FigureOutput/%s/', targetFig);
FIG_CODE_PATH = sprintf('./FigureCode/%s/', targetFig);

%% Make production figure panels
addpath(FIG_CODE_PATH);
mkdir(FIG_OUTPUT_PATH);
switch targetFig
    case 'fig1'
        makeFig1Plot(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'fig2'
        makeFig2Plots(FIG_DATA_PATH, FIG_OUTPUT_PATH);
    case 'sm-fig4'
        
    otherwise
        error('Target figure not found!');
end