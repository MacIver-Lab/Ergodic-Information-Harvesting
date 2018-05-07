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
targetFig = 'fig2';

%% Internal parameters (do not change)
DATA_PATH = '../SimulationCode/SimData/';
FIG_OUTPUT_PATH = sprintf('./FigureOutput/%s/', targetFig);
FIG_CODE_PATH = sprintf('./FigureCode/%s/', targetFig);
FIG_DATA_PATH = sprintf([DATA_PATH,'%s/'], targetFig);

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