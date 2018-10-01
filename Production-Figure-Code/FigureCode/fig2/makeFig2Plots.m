function makeFig2Plots(dataPath, savePath)
%% Plot individual panels for figure 2 insets
% Note that due to the complexity of this figure, each animal's inset will 
% be plotted separately into individual PDF files
% 
% Chen Chen

close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, fname);
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile(pwd, 'FigureCode', 'fig2', 'BehaviorData', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
barColor = [72, 110, 181;...
    50, 180, 74; ...
    236, 29, 36] / 255;
% Lambda function handle for computing cumulative 1D distance travelled
cumDist = @(x) sum(abs(diff(x)));
% Whether nor not to plot EER band
% set to 1 to enable EER plot overlay
% note that due to the complexity of EER bands, it's can be fairly slow 
% to plot the EER bands
global PLOT_EER_BAND
PLOT_EER_BAND = 0;
% Use split plot method to generate vector graphic plots
% This is a workaround for the buffer issue in MATLAB due
% to the EER patch is too complex for the interal save
% function to save as a vector graphic PDF
SPLIT_PLOT = 0;

%% Electric Fish Simulation
% Load data
EH_lSNR = load(GEN_DATA_PATH('fig2-ErgodicHarvest-ElectricFish-SNR-25-RandSeed-1.mat'), ...
    'oTrajList', 'sTrajList', 'dt', 'phi');
EH_hSNR = load(GEN_DATA_PATH('fig2-ErgodicHarvest-ElectricFish-SNR-60-RandSeed-1.mat'), ...
    'oTrajList', 'sTrajList', 'dt', 'phi');
EH_lSNR.eidList = flattenResultList(EH_lSNR.phi(:,:,1:end-1))';
EH_hSNR.eidList = flattenResultList(EH_hSNR.phi(:,:,1:end-1))';

% Load fish behavioral data
fish.hSNR = load(GEN_BEHAVIOR_DATA_PATH('ElectricFish-StrongSignal-Sine.mat'));
fish.lSNR = load(GEN_BEHAVIOR_DATA_PATH('ElectricFish-WeakSignal-Sine.mat'));

% Compute relative exploration
fish.hSNR.sDist = cumDist(fish.hSNR.fishTraj);
fish.hSNR.oDist = cumDist(fish.hSNR.refugeTraj);
fish.lSNR.sDist = cumDist(fish.lSNR.fishTraj);
fish.lSNR.oDist = cumDist(fish.lSNR.refugeTraj);

% Filter trajectory
simTrajHighCutFreq = 2.10;
EH_lSNR.sTrajList = LPF(EH_lSNR.sTrajList, 1/EH_lSNR.dt, simTrajHighCutFreq);
EH_hSNR.sTrajList = LPF(EH_hSNR.sTrajList, 1/EH_hSNR.dt, simTrajHighCutFreq);

% Cumulative 1D distance traveled
EH_hSNR.sDist = cumDist(EH_hSNR.sTrajList(1:end));
EH_hSNR.oDist = cumDist(EH_hSNR.oTrajList(1:end));
EH_lSNR.sDist = cumDist(EH_lSNR.sTrajList(1:end));
EH_lSNR.oDist = cumDist(EH_lSNR.oTrajList(1:end));

%--------- Fish Sinusoidal Tracking ---------%
% Strong Signal Trajectory plot
figure(1);clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
plot(fish.hSNR.refugeTraj, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(fish.hSNR.fishTraj, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(fish.hSNR.refugeTraj)];
opt.YLim = [mean(fish.hSNR.refugeTraj)-300, mean(fish.hSNR.refugeTraj)+300];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor(1:2,:);
setPlotProp(opt);
legend(gca, 'off');
set(gca,  'Position', [1    4    2.8320    1.7700]);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.6];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
plot(fish.lSNR.refugeTraj, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(fish.lSNR.fishTraj, ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(fish.lSNR.refugeTraj)];
opt.YLim = [mean(fish.lSNR.refugeTraj)-200, mean(fish.lSNR.refugeTraj)+200];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3],:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.6];
set(gca, 'Position', axesPosition);
% Relative Exploration bar plot
barAxes = axes; hold on;
barData = [1, fish.hSNR.sDist/fish.hSNR.oDist, fish.lSNR.sDist/fish.lSNR.oDist];
for i = 1:3
    bar(i, barData(i), 0.4, 'BaseValue', 0, ...
        'FaceColor', barColor(i,:));
end
opt = [];
opt.BoxDim = [8,5] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [1, 2, 3];
opt.YTick = [1, 3];
opt.YLim = [0, 3];
opt.FontSize = 8;
opt.FontName = 'Helvetica';
setAxesProp(opt, barAxes);
set(barAxes,'YTickLabel', {'1x', '3x'});
set(barAxes,'XTickLabel', {'Target', 'Strong Signal', 'Weak Signal'});
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.6];
set(gca, 'Position', axesPosition);

%--------- Ergodic Harvesting Simulation Trajectory ---------%
% Strong Signal Trajectory plot
trajAxes = axes; hold on;
plot(EH_hSNR.oTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(EH_hSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(EH_hSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor(1:2,:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(EH_hSNR);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.2, 0.3];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
plot(EH_lSNR.oTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(EH_lSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(EH_lSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3],:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(EH_lSNR);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.3];
set(gca, 'Position', axesPosition);
% Bar plot
barAxes = axes; hold on;
barData = [1, EH_hSNR.sDist/EH_hSNR.oDist, EH_lSNR.sDist/EH_lSNR.oDist];
for i = 1:3
    bar(i, barData(i), 0.4, 'BaseValue', 0, ...
        'FaceColor', barColor(i,:));
end
opt = [];
opt.BoxDim = [8,5] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [1, 2, 3];
opt.YTick = [1, 3];
opt.YLim = [0, 3];
opt.FontSize = 8;
opt.FontName = 'Helvetica';
setAxesProp(opt, barAxes);
set(barAxes,'YTickLabel', {'1x', '3x'});
set(barAxes,'XTickLabel', {'Target', 'Strong Signal', 'Weak Signal'});
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.3];
set(barAxes, 'Position', axesPosition);

% All set, now print the first section into PDF
if SPLIT_PLOT
    splitprint(gcf,... %separate the current figure
        GEN_SAVE_PATH('fig2-ElectricFish'),... %filenames will begin with 'disp2'
        {{'line';'text'},{'surface';'patch';'image'}}, ...% types of objects
        {'-dpdf','-dtiff'},... %file formats
        0,... %alignment mark will not be added
        [1 0],... %axes in first figure will be visible
        {'','-r400'});
else
    print(GEN_SAVE_PATH('fig2-ElectricFish.pdf'),'-dpdf');
end

%% Rat Odor Tracking
% Load Data
lSNR = load(GEN_DATA_PATH('fig2-ErgodicHarvest-Rat-WeakSignal.mat'));
hSNR = load(GEN_DATA_PATH('fig2-ErgodicHarvest-Rat-StrongSignal.mat'));
lSNR.eidList = flattenResultList(lSNR.phi(:,:,1:end-1))';
hSNR.eidList = flattenResultList(hSNR.phi(:,:,1:end-1))';

khan = load(GEN_BEHAVIOR_DATA_PATH('/Khan12a_fig2.mat'));
khan.lSNR.sTraj = khan.fig2b_nose;
khan.lSNR.oTraj = khan.fig2b_trail;
khan.hSNR.sTraj = khan.fig2a_nose;
khan.hSNR.oTraj = khan.fig2a_trail;

% Adjust time horizon to fit into the actual data length (provided in Khan12a)
lSNR.dt = 7.8947 / length(lSNR.sTrajList);
hSNR.dt = 6.8421 / length(hSNR.sTrajList);
khan.lSNR.dt = 7.8947 / length(khan.lSNR.sTraj);
khan.hSNR.dt = 6.8421 / length(khan.hSNR.sTraj);

% Compute cumulative 1D distance travelled
% Exclude initial global search before converge to ensure data consistency
% The criteria is to crop the initial searching trajectory until the sensor
% crosses the target for the first time
dist_hSNR_Sensor = cumDist(hSNR.sTrajList(1:end));
dist_hSNR_Trail = cumDist(hSNR.oTrajList(1:end));
dist_lSNR_Sensor = cumDist(lSNR.sTrajList(1:end));
dist_lSNR_Trail = cumDist(lSNR.oTrajList(1:end));
dist_rat_hSNR_Sensor = cumDist(khan.hSNR.sTraj);
dist_rat_hSNR_Trail = cumDist(khan.hSNR.oTraj);
dist_rat_lSNR_Sensor = cumDist(khan.lSNR.sTraj);
dist_rat_lSNR_Trail = cumDist(khan.lSNR.oTraj);

%--------- Rat Odor Tracking data from Khan12a ---------%
% Strong Signal Trajectory plot
figure(2);clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
plot(khan.hSNR.oTraj, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(khan.hSNR.sTraj, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(khan.hSNR.oTraj)];
opt.YLim = [mean(khan.hSNR.oTraj)-50, mean(khan.hSNR.oTraj)+50];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor(1:2,:);
setPlotProp(opt);
legend(gca, 'off');
set(gca,  'Position', [1    4    2.8320    1.7700]);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.6];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
plot(khan.lSNR.oTraj, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(khan.lSNR.sTraj, ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(khan.lSNR.oTraj)];
opt.YLim = [mean(khan.lSNR.oTraj)-50, mean(khan.lSNR.oTraj)+50];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3],:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.6];
set(gca, 'Position', axesPosition);
% Relative Exploration bar plot
barAxes = axes; hold on;
barData = [1, dist_rat_hSNR_Sensor/dist_rat_hSNR_Trail, dist_rat_lSNR_Sensor/dist_rat_lSNR_Trail];
for i = 1:3
    bar(i, barData(i), 0.4, 'BaseValue', -1, ...
        'FaceColor', barColor(i,:));
end
opt = [];
opt.BoxDim = [8,5] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [1, 2, 3];
opt.YTick = [1, 5];
opt.YLim = [-1, 5];
opt.FontSize = 8;
opt.FontName = 'Helvetica';
setAxesProp(opt, barAxes);
set(gca,'YTickLabel', {'1x', '5x'});
set(gca,'XTickLabel', {'Target', 'Strong Signal', 'Weak Signal'});
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.6];
set(gca, 'Position', axesPosition);

%--------- Ergodic Harvesting Simulation Trajectory ---------%
% Strong Signal Trajectory plot
trajAxes = axes; hold on;
plot(hSNR.oTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(hSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(hSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor(1:2,:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(hSNR);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.2, 0.3];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
plot(lSNR.oTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(lSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(lSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3],:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(lSNR);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.3];
set(gca, 'Position', axesPosition);
% Bar plot
barAxes = axes; hold on;
barData = [1, dist_hSNR_Sensor/dist_hSNR_Trail, dist_lSNR_Sensor/dist_lSNR_Trail];
for i = 1:3
    bar(i, barData(i), 0.4, 'BaseValue', -1, ...
        'FaceColor', barColor(i,:));
end
opt = [];
opt.BoxDim = [8,5] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [1, 2, 3];
opt.YTick = [1, 5];
opt.YLim = [-1, 5];
opt.FontSize = 8;
opt.FontName = 'Helvetica';
setAxesProp(opt, barAxes);
set(gca,'YTickLabel', {'1x', '5x'});
set(gca,'XTickLabel', {'Target', 'Strong Signal', 'Weak Signal'});
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.3];
set(gca, 'Position', axesPosition);

% All set, now print the first section into PDF
if SPLIT_PLOT
    splitprint(gcf,... %separate the current figure
        GEN_SAVE_PATH('fig2-Rat'),... %filenames will begin with 'disp2'
        {{'line';'text'},{'surface';'patch';'image'}}, ...% types of objects
        {'-dpdf','-dtiff'},... %file formats
        0,... %alignment mark will not be added
        [1 0],... %axes in first figure will be visible
        {'','-r400'});
else
    print(GEN_SAVE_PATH('fig2-Rat.pdf'),'-dpdf');
end

%% Mole Odor Localization
% Load Data
mole.hSNR = load(GEN_BEHAVIOR_DATA_PATH('Mole-StrongSignal.mat'));
mole.lSNR = load(GEN_BEHAVIOR_DATA_PATH('Mole-WeakSignal.mat'));
hSNR = load(GEN_DATA_PATH('fig2-ErgodicHarvest-Mole-StrongSignal-RandSeed-1.mat'));
lSNR = load(GEN_DATA_PATH('fig2-ErgodicHarvest-Mole-WeakSignal-RandSeed-1.mat'));

lSNR.eidList = flattenResultList(lSNR.phi(:,:,1:end))';
hSNR.eidList = flattenResultList(hSNR.phi(:,:,1:end))';

% Relative exploration
cumAngularDist = @(x) sum(abs(diff(x)));
hSNR.moleDist = cumAngularDist(hSNR.sTrajList);
lSNR.moleDist = cumAngularDist(lSNR.sTrajList);
mole.hSNR.moleDist = cumAngularDist(mole.hSNR.angleData);
mole.hSNR.refDist = cumAngularDist(mole.hSNR.angleData);
mole.lSNR.moleDist = cumAngularDist(mole.lSNR.angleData);
mole.lSNR.refDist = cumAngularDist(mole.lSNR.angleData);

hSNR.sTrajList = hSNR.sTrajList(1:1200);
lSNR.sTrajList = lSNR.sTrajList(1:1200);
%--------- Mole behavioral tracking data from Cata13a ---------%
% Strong Signal Trajectory plot
figure(3);clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
hLine = line([0, length(mole.hSNR.angleData)], [0, 0], ...
    'LineStyle', '--', ...
    'LineWidth', 2, ...
    'Color', 'k');
plot(mole.hSNR.angleData, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time');
ylabel('Angular Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(mole.hSNR.angleData)];
% opt.YLim = [mean(mole.hSNR.angleData)-50, mean(mole.hSNR.angleData)+50];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(2,:)];
setPlotProp(opt);
legend(gca, 'off');
hLine.LineStyle = '--';
set(gca,  'Position', [1    4    2.8320    1.7700]);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.6];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
hLine = line([0, length(mole.lSNR.angleData)], [0, 0], ...
    'LineStyle', '--', ...
    'LineWidth', 2, ...
    'Color', 'k');
plot(mole.lSNR.angleData, ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Time');
ylabel('Angular Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(mole.lSNR.angleData)];
% opt.YLim = [mean(mole.hSNR.angleData)-50, mean(mole.hSNR.angleData)+50];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(3,:)];
setPlotProp(opt);
legend(gca, 'off');
hLine.LineStyle = '--';
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.6];
set(gca, 'Position', axesPosition);
% Relative Exploration bar plot
barAxes = axes; hold on;
barData = [0, 1, mole.lSNR.moleDist/mole.hSNR.moleDist];
for i = 1:3
    bar(i, barData(i), 0.4, 'BaseValue', 0, ...
        'FaceColor', barColor(i,:));
end
opt = [];
opt.BoxDim = [8,5] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [1, 2, 3];
opt.YTick = [1, 2];
opt.YLim = [0.5, 2];
opt.FontSize = 8;
opt.FontName = 'Helvetica';
setAxesProp(opt, barAxes);
set(gca,'YTickLabel', {'1x', '2x'});
set(gca,'XTickLabel', {'Target', 'Strong Signal', 'Weak Signal'});
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.6];
set(gca, 'Position', axesPosition);

%--------- Ergodic Harvesting Simulation Trajectory ---------%
% Strong Signal Trajectory plot
trajAxes = axes; hold on;
hLine = line([0, length(hSNR.sTrajList)], [0.8, 0.8], ...
    'LineStyle', '--', ...
    'LineWidth', 2, ...
    'Color', 'k');
plot(hSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time');
ylabel('Angular Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(hSNR.sTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(2,:)];
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(hSNR);
hLine.LineStyle = '--';
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.2, 0.3];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
hLine = line([0, length(lSNR.sTrajList)], [0.5, 0.5], ...
    'LineStyle', '--', ...
    'LineWidth', 2, ...
    'Color', 'k');
plot(lSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Time');
ylabel('Angular Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(lSNR.sTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(3,:)];
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(lSNR);
hLine.LineStyle = '--';
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.3];
set(gca, 'Position', axesPosition);
% Bar plot
barAxes = axes; hold on;
barData = [0, 1, lSNR.moleDist/hSNR.moleDist];
for i = 1:3
    bar(i, barData(i), 0.4, 'BaseValue', 0, ...
        'FaceColor', barColor(i,:));
end
opt = [];
opt.BoxDim = [8,5] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [1, 2, 3];
opt.YTick = [1, 2];
opt.YLim = [0.5, 2];
opt.FontSize = 8;
opt.FontName = 'Helvetica';
setAxesProp(opt, barAxes);
set(gca,'YTickLabel', {'1x', '2x'});
set(gca,'XTickLabel', {'Target', 'Strong Signal', 'Weak Signal'});
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.3];
set(gca, 'Position', axesPosition);

% All set, now print the first section into PDF
if SPLIT_PLOT
    splitprint(gcf,... %separate the current figure
        GEN_SAVE_PATH('fig2-Mole'),... %filenames will begin with 'disp2'
        {{'line';'text'},{'surface';'patch';'image'}}, ...% types of objects
        {'-dpdf','-dtiff'},... %file formats
        0,... %alignment mark will not be added
        [1 0],... %axes in first figure will be visible
        {'','-r400'});
else
    print(GEN_SAVE_PATH('fig2-Mole.pdf'),'-dpdf');
end

fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));

%% Cockroach odor source localization
cockroachData = load(GEN_BEHAVIOR_DATA_PATH('cockroach_data.mat'));
% 4mm
[head, ~, src] = parseDataStream(cockroachData.cockroach_c4);
distRef = max([head(:, 1); src']) - max(head(1, 1), 0.0);
dist4 = cockroachData.cockroach_c4_table.walk_dist / distRef;
% Plot
figure(4); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
% Reference path
line([head(1, 1), src(1)], [0, 0], ...
    'LineStyle', '--', 'LineWidth', 4, 'Color', barColor(1, :));
% Strong signal
plot(head(:,1), head(:,2), ...
    'LineWidth', 2, ...
    'Color', barColor(2, :));
% Configure figure
xlabel('Position (X)');
ylabel('Position (Y)');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XLim = [head(1, 1), head(end, 1)];
opt.YLim = [-10, 10];
opt.XTick = [];
opt.YTick = [];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,2], :);
setPlotProp(opt);
legend(gca, 'off');
set(gca,  'Position', [1    4    2.8320    1.7700]);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.6];
set(gca, 'Position', axesPosition);
% 1mm
[head, ~, src] = parseDataStream(cockroachData.cockroach_c1);
distRef = max([head(:, 1); src']) - max(head(107, 1), 0.0);
dist1 = cockroachData.cockroach_c1_table.walk_dist / distRef;
trajAxes = axes; hold on;
% Reference path
line([head(107, 1), src(1)], [0, 0], ...
    'LineStyle', '--', 'LineWidth', 4, 'Color', barColor(1, :));
% Weak signal
plot(head(107:end-2,1), head(107:end-2,2), ...
    'LineWidth', 2, ...
    'Color', barColor(3, :));
% Configure figure
xlabel('Position (X)');
ylabel('Position (Y)');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XLim = [head(107, 1), head(end, 1)];
opt.YLim = [-10, 10];
opt.XTick = [];
opt.YTick = [];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3], :);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.6];
set(gca, 'Position', axesPosition);
% Barplot
barAxes = axes; hold on;
dist = [1, 1, dist1/dist4];
for i = 1:3
    bar(i, dist(i), 0.4, 'BaseValue', 0, ...
        'FaceColor', barColor(i,:));
end
opt = [];
opt.BoxDim = [8,5] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [1, 2, 3];
opt.YTick = [1, 3];
opt.YLim = [0, 3];
opt.FontSize = 8;
opt.FontName = 'Helvetica';
setAxesProp(opt, barAxes);
set(barAxes,'YTickLabel', {'1x', '3x'});
set(barAxes,'XTickLabel', {'Odor', 'Strong Signal', 'Weak Signal'});
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.6];
set(barAxes, 'Position', axesPosition);

%--------- Ergodic Harvesting Simulation Trajectory ---------%
lsnr = load(GEN_DATA_PATH('fig2-ErgodicHarvest-Cockroach-WeakSignal-SNR-25-RandSeed-1.mat'), ...
    'dt', 'sTrajList', 'oTrajList', 'phi');
hsnr = load(GEN_DATA_PATH('fig2-ErgodicHarvest-Cockroach-StrongSignal-SNR-80-RandSeed-1.mat'), ...
    'dt', 'sTrajList', 'oTrajList', 'phi');
trajLen = 1400;
hsnr.oTrajList = hsnr.oTrajList(1:trajLen);
hsnr.sTrajList = LPF(hsnr.sTrajList(1:trajLen), 1/hsnr.dt, 2);
lsnr.oTrajList = lsnr.oTrajList(1:trajLen);
lsnr.sTrajList = LPF(lsnr.sTrajList(1:trajLen), 1/lsnr.dt, 2);
hsnr.eidList = flattenResultList(hsnr.phi(:,:,1:end-1))';
lsnr.eidList = flattenResultList(lsnr.phi(:,:,1:end-1))';
hsnr.eidList = hsnr.eidList(:, 1:trajLen);
lsnr.eidList = lsnr.eidList(:, 1:trajLen);
% Plot
trajAxes = axes; hold on;
% Reference path
line([0, length(hsnr.oTrajList)], [0.5, 0.5], ...
    'LineStyle', '--', 'LineWidth', 4, 'Color', barColor(1, :));
% Strong signal
plot(hsnr.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(2, :));
% Configure figure
xlabel('Time');
ylabel('Lateral Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XLim = [0, length(hsnr.oTrajList)];
opt.YLim = [0.2, 0.8];
opt.XTick = [];
opt.YTick = [];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,2], :);
setPlotProp(opt);
legend(gca, 'off');
mPlotContinuousEID(hsnr);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.2, 0.3];
set(gca, 'Position', axesPosition);
% Weak Signal
trajAxes = axes; hold on;
% Reference path
line([0, length(hsnr.oTrajList)], [0.5, 0.5], ...
    'LineStyle', '--', 'LineWidth', 4, 'Color', barColor(1, :));
% Weak signal
plot(lsnr.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(3, :));
% Configure figure
xlabel('Time');
ylabel('Lateral Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XLim = [0, length(hsnr.oTrajList)];
opt.YLim = [0.2, 0.8];
opt.XTick = [];
opt.YTick = [];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3], :);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(lsnr);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.3];
set(gca, 'Position', axesPosition);
% Barplot
calcDist = @(x) sum(abs(diff(x)));
barAxes = axes; hold on;
dist = [1, 1, calcDist(lsnr.sTrajList) / calcDist(hsnr.sTrajList)];
for i = 1:3
    bar(i, dist(i), 0.4, 'BaseValue', 0, ...
        'FaceColor', barColor(i,:));
end
opt = [];
opt.BoxDim = [8,5] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [1, 2, 3];
opt.YTick = [1, 3];
opt.YLim = [0, 3];
opt.FontSize = 8;
opt.FontName = 'Helvetica';
setAxesProp(opt, barAxes);
set(barAxes,'YTickLabel', {'1x', '3x'});
set(barAxes,'XTickLabel', {'Odor', 'Strong Signal', 'Weak Signal'});
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.3];
set(barAxes, 'Position', axesPosition);

% All set, now print the first section into PDF
if SPLIT_PLOT
    splitprint(gcf,... %separate the current figure
        GEN_SAVE_PATH('fig2-Cockroach'),... %filenames will begin with 'disp2'
        {{'line';'text'},{'surface';'patch';'image'}}, ...% types of objects
        {'-dpdf','-dtiff'},... %file formats
        0,... %alignment mark will not be added
        [1 0],... %axes in first figure will be visible
        {'','-r400'});
else
    print(GEN_SAVE_PATH('fig2-Cockroach.pdf'),'-dpdf');
end



function mPlotContinuousEID(dat)
global PLOT_EER_BAND
if ~PLOT_EER_BAND
    return;
end
%% Plot Parameters
tScale = 5;   % Interval of EID plot update, set to 1 will plot all of the EID map
nBins = 80;   % Color resolution in the y axis
alpha = 0.5;  % Transparency of the EID color
% cmap = lines(10);
% cmap = cmap(7, :);
cmap = [0.7 0 0.4];

eidList = dat.eidList;
tRes = length(dat.oTrajList) / (size(eidList,2)-1);
sRes = size(eidList,1);
s = 1 / sRes;
faces = 1:4;

idxList = 1:tScale:floor(length(dat.oTrajList) / tRes);
for idx = 1:length(idxList)
    i = idxList(idx);
    [~,~,bin] = histcounts(eidList(:,i), nBins);
    for k = 1:sRes
        if bin(k) <= 2
            continue;
        end
        verts = [(i-tScale)*tRes, (k-1)*s;...
            (i-0)*tRes, (k-1)*s;...
            (i-0)*tRes, (k-0)*s;...
            (i-tScale)*tRes, (k-0)*s];
        patch('Faces',faces,'Vertices',verts,... 'FaceColor', [0.7 0 0.4],... 'FaceAlpha', alpha*bin(k)/nBins,...
            'FaceColor', cmap,...
            'FaceAlpha', alpha*bin(k)/nBins,...
            'EdgeColor', 'none');
%             'FaceColor', cmap(nBins-bin(k)+1,:),...
%             'FaceAlpha', alpha*bin(k)/nBins,...
%             'EdgeColor', 'none');
    end
    drawnow;
end


function outList = flattenResultList(list)
outList = zeros(size(list,2)*size(list,3), size(list,1));
for i = 1:size(list,3)
    for j = 1:size(list,2)
        outList((i-1)*size(list,2) + j,:) = list(:,j,i)';
    end
end

function [head, tail, src] = parseDataStream(dat)
head = dat(:, 1:2);
tail = dat(:, 3:4);
src = dat(1, 5:6);

function [m, sem] = calc95CI(x)
m = mean(x);
sem = 1.96 * std(x) / sqrt(length(x));

function re = calcRE(dat)
% Computes relative exploration
len = size(dat, 1);
re = zeros(len, 1, 'double');
for i = 1:len
    re(i) = dat.walk_dist(i) / dat.distRef(i);
end
