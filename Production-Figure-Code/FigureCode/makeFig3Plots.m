function makeFig3Plots(dataPath, savePath)
%% Plot individual panels for figure 3 insets
% Note that due to the complexity of this figure, each animal's inset will 
% be plotted separately into individual PDF files
% 
% Chen Chen

close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, 'behavior_sim', fname);
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile('./PublishedData/', 'animal_behavior_data', fname);
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
PLOT_EER_BAND = 1;
% Use split plot method to generate vector graphic plots
% This is a workaround for the buffer issue in MATLAB due
% to the EER patch is too complex for the interal save
% function to save as a vector graphic PDF
SPLIT_PLOT = 0;

%% Electric Fish Simulation
% Load data
EH_lSNR = load(GEN_DATA_PATH('EIH-ElectricFish-WeakSignal-RandSeed-1.mat'), ...
    'oTrajList', 'sTrajList', 'dt', 'phi');
EH_hSNR = load(GEN_DATA_PATH('EIH-ElectricFish-StrongSignal-RandSeed-1.mat'), ...
    'oTrajList', 'sTrajList', 'dt', 'phi');
EH_lSNR.eidList = flattenResultList(EH_lSNR.phi(:,:,1:end-1))';
EH_hSNR.eidList = flattenResultList(EH_hSNR.phi(:,:,1:end-1))';


% Load fish behavioral data
ss = load(GEN_BEHAVIOR_DATA_PATH('/ElectricFish/Eigenmannia-sp-StrongSignal.mat'), 'strongSigData');
fish.hSNR.fishTraj = ss.strongSigData{8, 1};
fish.hSNR.refugeTraj = ss.strongSigData{8, 2};
dt = 1 / 60;
trajSegStr = 31.5;
trajSegLen = 30;
trajSegStrIdx = trajSegStr * 60;
trajSegIdx = (trajSegStr+trajSegLen) * 60;
timeIdx = 0:dt:trajSegLen*60*dt;
fish.hSNR.fishTraj = fish.hSNR.fishTraj(trajSegStrIdx:trajSegIdx);
fish.hSNR.fishTraj = fish.hSNR.fishTraj - mean(fish.hSNR.fishTraj);
fish.hSNR.refugeTraj = fish.hSNR.refugeTraj(trajSegStrIdx:trajSegIdx);
fish.hSNR.refugeTraj = fish.hSNR.refugeTraj - mean(fish.hSNR.refugeTraj);

ws = load(GEN_BEHAVIOR_DATA_PATH('/ElectricFish/Eigenmannia-sp-WeakSignal.mat'), 'weakSigData');
fish.lSNR.fishTraj = ws.weakSigData{1, 1};
fish.lSNR.fishTraj = fish.lSNR.fishTraj - mean(fish.lSNR.fishTraj);
fish.lSNR.refugeTraj = ws.weakSigData{1, 2};
fish.lSNR.refugeTraj = fish.lSNR.refugeTraj - mean(fish.lSNR.refugeTraj);
% Filter trajectory
simTrajHighCutFreq = 2.10;
EH_lSNR.sTrajList = LPF(EH_lSNR.sTrajList, 1/EH_lSNR.dt, simTrajHighCutFreq);
EH_hSNR.sTrajList = LPF(EH_hSNR.sTrajList, 1/EH_hSNR.dt, simTrajHighCutFreq);

% calibrate distance unit
pix2cm = 1.7 / (max(fish.hSNR.refugeTraj) - min(fish.hSNR.refugeTraj));

%--------- Fish Sinusoidal Tracking ---------%
% Strong Signal Trajectory plot
figure(1);clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
plot(timeIdx, fish.hSNR.refugeTraj*pix2cm, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(timeIdx, fish.hSNR.fishTraj*pix2cm, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time (s)');
ylabel('Position (cm)');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [0, 30];
opt.YTick = -2:2:2;
opt.YLim = [-2, 2];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor(1:2,:);
setAxesProp(opt, gca);
legend(gca, 'off');
set(gca,  'Position', [1    4    2.8320    1.7700]);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.6];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
plot(timeIdx, fish.lSNR.refugeTraj*pix2cm, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(timeIdx, fish.lSNR.fishTraj*pix2cm, ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Time (s)');
ylabel('Position (cm)');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XTick = [0, 30];
opt.YTick = -2:2:2;
opt.YLim = [-2, 2];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3],:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.6];
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

% Statistics plot
FPS = 60;
cumDist = @(x) sum(abs(diff(x)));
disp('*****************************************************');
disp('Electric fish behavioral tracking data');
% Strong Signal, n = 10 trials, with volitational motion cropped
load(GEN_BEHAVIOR_DATA_PATH('/ElectricFish/Eigenmannia-sp-StrongSignal.mat'), 'averageTrialLenInSec', 'trialLength', 'strongSigData');
nTrials = length(trialLength);
strongSig = zeros(nTrials, 2, 'double');
for i = 1:nTrials
    strongSig(i,1) = cumDist(strongSigData{i,1});
    strongSig(i,2) = cumDist(strongSigData{i,2});
end
fprintf('Strong signal trials, n = %d, average trial length = %.2f seconds\n', ...
    length(trialLength), averageTrialLenInSec);

% Weak Signal, n = 11
load(GEN_BEHAVIOR_DATA_PATH('/ElectricFish/Eigenmannia-sp-WeakSignal.mat'), 'averageTrialLenInSec', 'trialLength', 'weakSigData');
nTrials = length(trialLength);
weakSig = zeros(nTrials, 2, 'double');
for i = 1:nTrials
    weakSig(i,1) = cumDist(weakSigData{i,1});
    weakSig(i,2) = cumDist(weakSigData{i,2});
end
fprintf('Weak signal trials, n = %d, average trial length = %.2f seconds\n', ...
    length(trialLength), averageTrialLenInSec);

% Statistical Analysis
% Compute relative exploration (re)
reStrongSignalFish = strongSig(:,1) ./ strongSig(:, 2);
reWeakSignalFish = weakSig(:,1) ./ weakSig(:, 2);

% Ergodic Harvesting Data
% Load Ergodic data
EH_lSNR_files = dir(GEN_DATA_PATH('EIH-ElectricFish-WeakSignal-RandSeed-*.mat'));
EH_hSNR_files = dir(GEN_DATA_PATH('EIH-ElectricFish-StrongSignal-RandSeed-*.mat'));
IT_lSNR_files = dir(GEN_DATA_PATH('Infotaxis-ElectricFish-WeakSignal-RandSeed-*.mat'));
IT_hSNR_files = dir(GEN_DATA_PATH('Infotaxis-ElectricFish-StrongSignal-RandSeed-*.mat'));

reErgSS = zeros(1, length(EH_lSNR_files), 'double');
reErgWS = zeros(1, length(EH_hSNR_files), 'double');
reInfSS = zeros(1, length(IT_lSNR_files), 'double');
reInfWS = zeros(1, length(IT_hSNR_files), 'double');

for i = 1:length(EH_lSNR_files)
    EH_lSNR = load(GEN_DATA_PATH(EH_lSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    EH_hSNR = load(GEN_DATA_PATH(EH_hSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    IT_lSNR = load(GEN_DATA_PATH(IT_lSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    IT_hSNR = load(GEN_DATA_PATH(IT_hSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    
    % Filter trajectory
    trajHighCutFreq = 2.10;
    EH_lSNR.sTrajList = LPF(EH_lSNR.sTrajList, 1/EH_lSNR.dt, trajHighCutFreq);
    EH_hSNR.sTrajList = LPF(EH_hSNR.sTrajList, 1/EH_hSNR.dt, trajHighCutFreq);
    IT_lSNR.sTrajList = LPF(IT_lSNR.sTrajList, 1/IT_lSNR.dt, trajHighCutFreq);
    IT_hSNR.sTrajList = LPF(IT_hSNR.sTrajList, 1/IT_hSNR.dt, trajHighCutFreq);
    
    % Cumulative 1D distance traveled
    EH_hSNR.sPower = cumDist(EH_hSNR.sTrajList(100:end));
    EH_hSNR.oPower = cumDist(EH_hSNR.oTrajList(100:end));
    EH_lSNR.sPower = cumDist(EH_lSNR.sTrajList(100:end));
    EH_lSNR.oPower = cumDist(EH_lSNR.oTrajList(100:end));
    IT_hSNR.sPower = cumDist(IT_hSNR.sTrajList(100:end));
    IT_hSNR.oPower = cumDist(IT_hSNR.oTrajList(100:end));
    IT_lSNR.sPower = cumDist(IT_lSNR.sTrajList(100:end));
    IT_lSNR.oPower = cumDist(IT_lSNR.oTrajList(100:end));
    
    reErgSS(i) = EH_hSNR.sPower/EH_hSNR.oPower;
    reErgWS(i) = EH_lSNR.sPower/EH_lSNR.oPower;
    reInfSS(i) = IT_hSNR.sPower/IT_hSNR.oPower;
    reInfWS(i) = IT_lSNR.sPower/IT_lSNR.oPower;
end
reErgSS_Mean = mean(reErgSS);
reErgWS_Mean = mean(reErgWS);
reInfSS_Mean = mean(reInfSS);
reInfWS_Mean = mean(reInfWS);
reErgSS_SEM = 1.96 * std(reErgSS) / sqrt(length(reErgSS));
reErgWS_SEM = 1.96 * std(reErgWS) / sqrt(length(reErgWS));
reInfSS_SEM = 1.96 * std(reInfSS) / sqrt(length(reInfSS));
reInfWS_SEM = 1.96 * std(reInfWS) / sqrt(length(reInfWS));

fprintf(['---------------------------------------------------------\n',...
    '|         Relative exploration (mean +/- 95%% CI)        |\n', ...
    '|--------------------------------------------------------|',...
    '\n|> Ergodic Harvesting (strong signal) = %.3f +/- %.3f <|', ...
    '\n|> Ergodic Harvesting (weak signal)   = %.3f +/- %.3f <|',...
    '\n|>     Infotaxis (strong signal)      = %.3f +/- %.3f <|', ...
    '\n|>     Infotaxis (weak signal)        = %.3f +/- %.3f <|\n',...
    '----------------------------------------------------------\n'], ...
    reErgSS_Mean, reErgSS_SEM, reErgWS_Mean, reErgWS_SEM, ...
    reInfSS_Mean, reInfSS_SEM, reInfWS_Mean, reInfWS_SEM);


% Compute statistics
% Ranksum test, right tail = median of A is greater than median of B
[P, ~, Stats] = ranksum(reWeakSignalFish, reStrongSignalFish,...
    'tail', 'right');
fprintf('Behavior statistics: Wilcoxon rank sum test (one-sided) - p = %.4f (n = %d)\n', ...
    P, length([reWeakSignalFish; reStrongSignalFish]));
[P, ~, Stats] = ranksum(reErgWS, reErgSS,...
    'tail', 'right');
fprintf('EIH statistics: Wilcoxon rank sum test (one-sided) - p = %.4f (n = %d)\n', ...
    P, length([reErgWS, reErgSS]));


% Plot group data
axes, hold on;
hL = line([0,3], [1,1], 'LineStyle', '--', ...
    'Color', [140,140,140]/255.0, 'LineWidth',4);
hBoxPlot = notBoxPlot([reStrongSignalFish',reWeakSignalFish'], ...
    [1.1*ones(1,length(reStrongSignalFish)), 1.9*ones(1,length(reWeakSignalFish))], ...
    'jitter', 0.04, 'alpha', 0.5, 'dotSize', 24);
opt = [];
opt.BoxDim = [8,5] * 0.354;
opt.YLabel = 'Relative Exploration'; % ylabel
opt.XLim = [0.5, 2.5];
opt.YLim = [0.75, 3.5];
opt.YTick = [1, 2, 3];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontSize = 13;
opt.FontName = 'Helvetica';
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hL.Color = [140,140,140]/255.0;
hL.LineWidth = 4;
hL.LineStyle = '--';
hLine(1).LineWidth = 3;
hLine(2).LineWidth = 3;
hLine(1).Color = [162,0,0]/255.0;
hLine(2).Color = [50,180,74]/255.0;
legend(gca, 'off');
set(gca,'XTickLabel',{...
    '\color[rgb]{0.1961,0.7059,0.2902}Strong Signal',...
    '\color[rgb]{0.6353,0,0}Weak Signal'})
set(gca,'YTickLabel',{'1x', '2x', '3x'})
set(gca, 'units', 'normalized');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.6];
set(gca, 'Position', axesPosition);

% EIH and Infotaxis
axes; hold on;
hL = line([-1,5], [1,1], 'LineStyle', '--', ...
    'Color', [140,140,140]/255.0, 'LineWidth',4);
hBoxPlot = notBoxPlot([reErgSS, reInfSS, reErgWS, reInfWS], ...
    [0.5*ones(1,length(reErgSS)), ...
    0.75*ones(1,length(reInfSS)), ...
    1.25*ones(1,length(reErgWS)),...
    1.5*ones(1,length(reInfWS))], ...
    'jitter', 0.04, 'alpha', 0.5, 'dotSize', 24);
cmap = [...
    85, 1, 159; ...
    255, 196, 3;...
    85, 1, 159; ...
    255, 196, 3;...
    ] ./ 255;
opt = [];
opt.BoxDim = [8,5]*0.354;
opt.YLabel = 'Relative Exploration'; % ylabel
opt.YLim = [0.75, 3.5];
opt.YTick = [1, 2, 3];
opt.XLim = [0,2];
opt.XTick = [0.625, 1.375];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontSize = 13;
opt.FontName = 'Helvetica';
opt.Colors = cmap;
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hL.Color = [140,140,140]/255.0;
hL.LineWidth = 4;
hL.LineStyle = '--';
hLine(1).LineWidth = 3;
hLine(2).LineWidth = 3;
hLine(3).LineWidth = 3;
hLine(4).LineWidth = 3;
hLine(1).Color = cmap(1, :);
hLine(2).Color = cmap(2, :);
hLine(3).Color = cmap(1, :);
hLine(4).Color = cmap(2, :);
legend(gca, 'off');
set(gca,'XTickLabel',{...
    '\color[rgb]{0.1961,0.7059,0.2902}Strong Signal',...
    '\color[rgb]{0.6353,0,0}Weak Signal'})
set(gca,'YTickLabel',{'1x', '2x', '3x'})
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.6, 0.3];
set(gca, 'Position', axesPosition);


% All set, now print the first section into PDF
if SPLIT_PLOT
    splitprint(gcf,... %separate the current figure
        GEN_SAVE_PATH('fig3-ElectricFish'),... %filenames will begin with 'disp2'
        {{'line';'text'},{'surface';'patch';'image'}}, ...% types of objects
        {'-dpdf','-dtiff'},... %file formats
        0,... %alignment mark will not be added
        [1 0],... %axes in first figure will be visible
        {'','-r400'});
else
    print(GEN_SAVE_PATH('fig3-ElectricFish.pdf'),'-dpdf');
end

%% Mole Odor Localization
% Load Data
mole.hSNR = load(GEN_BEHAVIOR_DATA_PATH('/Mole/fig2-normal-right-green.mat'));
mole.lSNR = load(GEN_BEHAVIOR_DATA_PATH('/Mole/fig2-rBlock-blue.mat'));
mole.hSNR = processMoleTraj(mole.hSNR);
mole.lSNR = processMoleTraj(mole.lSNR);

% EIH
hSNR = load(GEN_DATA_PATH('EIH-Mole-StrongSignal-RandSeed-1.mat'));
lSNR = load(GEN_DATA_PATH('EIH-Mole-WeakSignal-RandSeed-1.mat'));
lSNR.eidList = flattenResultList(lSNR.phi(:,:,1:end))';
hSNR.eidList = flattenResultList(hSNR.phi(:,:,1:end))';

hSNR.sTrajList = hSNR.sTrajList(1:1200);
lSNR.sTrajList = lSNR.sTrajList(1:1200);
%--------- Mole behavioral tracking data from Cata13a ---------%
% y axis ratio [cm/pts]
yAxisUnitRatio =  6.89 / (max(mole.lSNR.molePath(end:-1:1, 1)) - min(mole.lSNR.molePath(end:-1:1, 1)));
% Strong Signal Trajectory plot
figure(3); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
hLine = line([0, length(mole.hSNR.rotPath(:, 1))], ...
    [mole.hSNR.rotRef(1, 1), mole.hSNR.rotRef(2, 1)]*yAxisUnitRatio, ...
    'LineStyle', '--', ...
    'LineWidth', 2, ...
    'Color', 'k');
timeStamps = linspace(0, 6, length(mole.hSNR.molePath(end:-1:1, 1)));
plot(timeStamps, mole.hSNR.molePath(end:-1:1, 1)*yAxisUnitRatio, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time (s)');
ylabel('Lateral Position (cm)');
ylim(([295, 365]-19)*yAxisUnitRatio);
yRange = ylim;
yTickInc = 5;
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [0, max(timeStamps)];
opt.YTick = [mean(yRange)-yTickInc, mean(yRange)+yTickInc];
opt.XLim = [0, max(timeStamps)];
opt.YLim = [mean(yRange)-yTickInc, mean(yRange)+yTickInc];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(2,:)];
setAxesProp(opt, gca);
legend(gca, 'off');
set(gca, 'YTickLabel', {'-5', '5'});
hLine.LineStyle = '--';
set(gca,  'Position', [1    4    2.8320    1.7700]);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.6];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
hLine = line([0, length(mole.lSNR.rotPath(:, 1))], ...
    ([mole.lSNR.rotRef(1, 1), mole.lSNR.rotRef(2, 1)]+19.5)*yAxisUnitRatio, ...
    'LineStyle', '--', ...
    'LineWidth', 2, ...
    'Color', 'k');
timeStamps = linspace(0, 9, length(mole.lSNR.molePath(end:-1:1, 1)));
plot(timeStamps, mole.lSNR.rotPath(end:-1:1, 1)*yAxisUnitRatio, ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Time (s)');
ylabel('Lateral Position (cm)');
ylim(([295, 365]+4)*yAxisUnitRatio);
yRange = ylim;
yTickInc = 5;
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [0, max(timeStamps)];
opt.YTick = [mean(yRange)-yTickInc, mean(yRange)+yTickInc];
opt.XLim = [0, max(timeStamps)];
opt.YLim = [mean(yRange)-yTickInc, mean(yRange)+yTickInc];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(3,:)];
setAxesProp(opt, gca);
legend(gca, 'off');
set(gca, 'YTickLabel', {'-5', '5'});
hLine.LineStyle = '--';
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.6];
set(gca, 'Position', axesPosition);

%--------- Ergodic Harvesting Simulation Trajectory ---------%
% Strong Signal Trajectory plot
trajAxes = axes; hold on;
hLine = line([0, length(hSNR.sTrajList)], [0.5, 0.5], ...
    'LineStyle', '--', ...
    'LineWidth', 2, ...
    'Color', 'k');
plot(hSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time');
ylabel('Lateral Position');
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
ylabel('Lateral Position');
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

% Statictics plot
cumAngularDist = @(x) sum(abs(diff(x)));
disp('*****************************************************');
disp('Mole behavioral tracking data from Cata13a');
fNames = dir(GEN_BEHAVIOR_DATA_PATH('/Mole/fig*.mat'));
nFiles = length(fNames);
reNrm = [];
reBlk = [];
for i = 1:nFiles
    load(GEN_BEHAVIOR_DATA_PATH(['/Mole/', fNames(i).name]));
    
    % Process data
    if ~contains(fNames(i).name,'normal')
        %%%%%%- One-side nostril block - Low SNR data -%%%%%%
        % Relative exploration
        reBlk = [reBlk, calcSumLength2D(molePath)/calcSumLength2D(refPath)];
    else
        %%%%%%- Normal - High SNR data -%%%%%%
        % Relative exploration
        reNrm = [reNrm, calcSumLength2D(molePath)/calcSumLength2D(refPath)];
    end
end
fprintf('Strong signal trials, n = %d\n', length(reNrm));
fprintf('Weak signal trials, n = %d\n', length(reBlk));

% Mole Odor Localization
% Ergodic Harvesting
EH_lSNR_files = dir(GEN_DATA_PATH('EIH-Mole-WeakSignal-RandSeed-*.mat'));
EH_hSNR_files = dir(GEN_DATA_PATH('EIH-Mole-StrongSignal-RandSeed-*.mat'));
IT_lSNR_files = dir(GEN_DATA_PATH('Infotaxis-Mole-WeakSignal-RandSeed-*.mat'));
IT_hSNR_files = dir(GEN_DATA_PATH('Infotaxis-Mole-StrongSignal-RandSeed-*.mat'));

reErgSS = zeros(1, length(EH_lSNR_files), 'double');
reErgWS = zeros(1, length(EH_hSNR_files), 'double');
reInfSS = zeros(1, length(IT_lSNR_files), 'double');
reInfWS = zeros(1, length(IT_hSNR_files), 'double');


for i = 1:length(EH_lSNR_files)
    EH_lSNR = load(GEN_DATA_PATH(EH_lSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    EH_hSNR = load(GEN_DATA_PATH(EH_hSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    IT_lSNR = load(GEN_DATA_PATH(IT_lSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    IT_hSNR = load(GEN_DATA_PATH(IT_hSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    
    % Filter trajectory
    trajHighCutFreq = 3.0;
    EH_lSNR.sTrajList = LPF(EH_lSNR.sTrajList, 1/EH_lSNR.dt, trajHighCutFreq);
    EH_hSNR.sTrajList = LPF(EH_hSNR.sTrajList, 1/EH_hSNR.dt, trajHighCutFreq);
    IT_lSNR.sTrajList = LPF(IT_lSNR.sTrajList, 1/IT_lSNR.dt, trajHighCutFreq);
    IT_hSNR.sTrajList = LPF(IT_hSNR.sTrajList, 1/IT_hSNR.dt, trajHighCutFreq);
    
    % Reference distance (straight to the target)
    EH_lSNR.refDist = abs(EH_lSNR.sTrajList(1) - EH_lSNR.oTrajList(1));
    EH_hSNR.refDist = abs(EH_hSNR.sTrajList(1) - EH_hSNR.oTrajList(1));
    IT_lSNR.refDist = abs(IT_lSNR.sTrajList(1) - IT_lSNR.oTrajList(1));
    IT_hSNR.refDist = abs(IT_hSNR.sTrajList(1) - IT_hSNR.oTrajList(1));
    
    % Relative exploration effort - Angular distance traveled
    % crop to ignore the distance due to different initial 
    % position
    EH_hSNR.moleDist = cumAngularDist(EH_hSNR.sTrajList(100:end));
    EH_lSNR.moleDist = cumAngularDist(EH_lSNR.sTrajList(100:end));
    IT_hSNR.moleDist = cumAngularDist(IT_hSNR.sTrajList(100:end));
    IT_lSNR.moleDist = cumAngularDist(IT_lSNR.sTrajList(100:end));
    
    reErgSS(i) = EH_hSNR.moleDist;
    reErgWS(i) = EH_lSNR.moleDist;
    reInfSS(i) = IT_hSNR.moleDist;
    reInfWS(i) = IT_lSNR.moleDist;
end
reErgSS_Mean = mean(reErgSS);
reErgWS_Mean = mean(reErgWS);
reInfSS_Mean = mean(reInfSS);
reInfWS_Mean = mean(reInfWS);
reErgSS_SEM = 1.96 * std(reErgSS) / sqrt(length(reErgSS));
reErgWS_SEM = 1.96 * std(reErgWS) / sqrt(length(reErgWS));
reInfSS_SEM = 1.96 * std(reInfSS) / sqrt(length(reInfSS));
reInfWS_SEM = 1.96 * std(reInfWS) / sqrt(length(reInfWS));

fprintf(['---------------------------------------------------------\n',...
    '|         Relative exploration (mean +/- 95%% CI)        |\n', ...
    '|--------------------------------------------------------|',...
    '\n|> Ergodic Harvesting (strong signal) = %.3f +/- %.3f <|', ...
    '\n|> Ergodic Harvesting (weak signal)   = %.3f +/- %.3f <|',...
    '\n|>     Infotaxis (strong signal)      = %.3f +/- %.3f <|', ...
    '\n|>     Infotaxis (weak signal)        = %.3f +/- %.3f <|\n',...
    '----------------------------------------------------------\n'], ...
    reErgSS_Mean, reErgSS_SEM, reErgWS_Mean, reErgWS_SEM, ...
    reInfSS_Mean, reInfSS_SEM, reInfWS_Mean, reInfWS_SEM);

% Compute statistics
% Ranksum test, right tail = median of A is greater than median of B
[P, ~, Stats] = ranksum(reBlk, reNrm,...
    'tail', 'right');
fprintf('Behavior statistics: Wilcoxon rank sum test (one-sided) - p = %.4f (n = %d)\n', ...
    P, length([reBlk, reNrm]));
[P, ~, Stats] = ranksum(reErgWS, reErgSS,...
    'tail', 'right');
fprintf('EIH statistics: Wilcoxon rank sum test (one-sided) - p = %.4f (n = %d)\n', ...
    P, length([reErgWS, reErgSS]));


% Plot group data
axes; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
hBoxPlot = notBoxPlot([reBlk,reNrm], ...
    [1.9*ones(1,length(reBlk)), 1.1*ones(1,length(reNrm))], ...
    'jitter', 0.04, 'alpha', 0.5, 'DotSize', 24);
opt = [];
opt.BoxDim = [8,5]*0.354;
opt.YLabel = 'Relative Exploration'; % ylabel
opt.XLim = [0.5, 2.5];
opt.YLim = [1, 6];
opt.YTick = [1, 2, 4, 6];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 13;
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hLine(1).LineWidth = 3;
hLine(2).LineWidth = 3;
hLine(1).Color = [162,0,0]/255.0;
hLine(2).Color = [50,180,74]/255.0;
legend(gca, 'off');
set(gca,'XTickLabel',{...
    '\color[rgb]{0.1961,0.7059,0.2902}Strong Signal',...
    '\color[rgb]{0.6353,0,0}Weak Signal'})
set(gca,'YTickLabel',{'1x', '2x', '4x', '6x'})
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.6, 0.6];
set(gca, 'Position', axesPosition);

% EIH and Infotaxis data
axes; hold on;
hBoxPlot = notBoxPlot([reErgSS, reInfSS, reErgWS, reInfWS], ...
    [0.5*ones(1,length(reErgSS)), ...
    0.75*ones(1,length(reInfSS)), ...
    1.25*ones(1,length(reErgWS)),...
    1.5*ones(1,length(reInfWS))], ...
    'jitter', 0.04, 'alpha', 0.5, 'DotSize', 24);
cmap = [...
    85, 1, 159; ...
    255, 196, 3;...
    85, 1, 159; ...
    255, 196, 3;...
    ] ./ 255;
opt = [];
opt.BoxDim = [8,5]*0.354;
opt.YLabel = 'Exploration'; % ylabel
opt.YLim = [0, 6];
opt.YTick = [0, 2, 4, 6];
opt.XLim = [0.2, 1.8];
opt.XTick = [0.625, 1.375];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 13;
opt.Colors = cmap;
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hLine(1).LineWidth = 3;
hLine(2).LineWidth = 3;
hLine(3).LineWidth = 3;
hLine(4).LineWidth = 3;
hLine(1).Color = cmap(1, :);
hLine(2).Color = cmap(2, :);
hLine(3).Color = cmap(1, :);
hLine(4).Color = cmap(2, :);
legend(gca, 'off');
set(gca,'XTickLabel',{...
    '\color[rgb]{0.1961,0.7059,0.2902}Strong Signal',...
    '\color[rgb]{0.6353,0,0}Weak Signal'})
set(gca,'YTickLabel',{'0', '2', '4', '6'})
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.6, 0.3];
set(gca, 'Position', axesPosition);

% All set, now print the first section into PDF
if SPLIT_PLOT
    splitprint(gcf,... %separate the current figure
        GEN_SAVE_PATH('fig3-Mole'),... %filenames will begin with 'disp2'
        {{'line';'text'},{'surface';'patch';'image'}}, ...% types of objects
        {'-dpdf','-dtiff'},... %file formats
        0,... %alignment mark will not be added
        [1 0],... %axes in first figure will be visible
        {'','-r400'});
else
    print(GEN_SAVE_PATH('fig3-Mole.pdf'),'-dpdf');
end

fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));

%% Cockroach odor source localization
cockroachData = load(GEN_BEHAVIOR_DATA_PATH('Cockroach/cockroach_data.mat'));
% 4mm
[head, ~, src] = parseDataStream(cockroachData.cockroach_c4);
dt = cockroachData.cockroach_c1_table.time_to_source / size(head,1);
head = head(5:end, :);
timeStamps = 0:dt:(size(head,1)-1)*dt;
% Plot
figure(4); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
% Reference path
line([0, src(1)], [0, 0], ...
    'LineStyle', '--', 'LineWidth', 4, 'Color', barColor(1, :));
% Strong signal
plot(timeStamps, head(:,2)-mean(head(:,2)), ...
    'LineWidth', 2, ...
    'Color', barColor(2, :));
% Configure figure
xlabel('Time (s)');
ylabel('Lateral Position (cm)');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XLim = [0, 10];
opt.YLim = [-8, 8];
opt.XTick = [0, 10];
opt.YTick = [-8, 8];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,2], :);
setAxesProp(opt, gca);
legend(gca, 'off');
set(gca,  'Position', [1    4    2.8320    1.7700]);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.6];
set(gca, 'Position', axesPosition);
% 1mm
[head, ~, src] = parseDataStream(cockroachData.cockroach_c1);
dt = cockroachData.cockroach_c1_table.time_to_source / size(head,1);
head = head(107:end-2, :);
timeStamps = 0:dt:(size(head,1)-1)*dt;
trajAxes = axes; hold on;
% Reference path
line([0, src(1)], [0, 0], ...
    'LineStyle', '--', 'LineWidth', 4, 'Color', barColor(1, :));
% Weak signal
plot(timeStamps, head(:,2)-mean(head(:,2)), ...
    'LineWidth', 2, ...
    'Color', barColor(3, :));
% Configure figure
xlabel('Time (s)');
ylabel('Lateral Position (cm)');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XLim = [0, 3.5];
opt.YLim = [-8, 8];
opt.XTick = [0, 3.5];
opt.YTick = [-8, 8];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3], :);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.6];
set(gca, 'Position', axesPosition);

%--------- Ergodic Harvesting Simulation Trajectory ---------%
lsnr = load(GEN_DATA_PATH('EIH-Cockroach-WeakSignal-RandSeed-1.mat'), ...
    'dt', 'sTrajList', 'oTrajList', 'phi');
hsnr = load(GEN_DATA_PATH('EIH-Cockroach-StrongSignal-RandSeed-1.mat'), ...
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
setAxesProp(opt, trajAxes);
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

% Statictics plot
disp('*****************************************************');
disp('Cockroach behavioral tracking data from Lock15a');
cockroachData = load(GEN_BEHAVIOR_DATA_PATH('/Cockroach/cockroach_data.mat'));
% 4mm trials
re4 = calcRE(cockroachData.trial_c4);
% 1, 2 mm trials
re12 = calcRE(cockroachData.trial_c12);
fprintf('Strong signal trials, n = %d\n', length(re4));
fprintf('Weak signal trials, n = %d\n', length(re12));

% Mole Odor Localization
% Ergodic Harvesting
EH_lSNR_files = dir(GEN_DATA_PATH('EIH-Cockroach-WeakSignal*.mat'));
EH_hSNR_files = dir(GEN_DATA_PATH('EIH-Cockroach-StrongSignal*.mat'));
IT_lSNR_files = dir(GEN_DATA_PATH('Infotaxis-Cockroach-WeakSignal*.mat'));
IT_hSNR_files = dir(GEN_DATA_PATH('Infotaxis-Cockroach-StrongSignal*.mat'));

reErgSS = zeros(1, length(EH_lSNR_files), 'double');
reErgWS = zeros(1, length(EH_hSNR_files), 'double');
reInfSS = zeros(1, length(IT_lSNR_files), 'double');
reInfWS = zeros(1, length(IT_hSNR_files), 'double');


for i = 1:length(EH_lSNR_files)
    EH_lSNR = load(GEN_DATA_PATH(EH_lSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    EH_hSNR = load(GEN_DATA_PATH(EH_hSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    IT_lSNR = load(GEN_DATA_PATH(IT_lSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    IT_hSNR = load(GEN_DATA_PATH(IT_hSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    
    % Filter trajectory
    trajHighCutFreq = 3.0;
    EH_lSNR.sTrajList = LPF(EH_lSNR.sTrajList, 1/EH_lSNR.dt, trajHighCutFreq);
    EH_hSNR.sTrajList = LPF(EH_hSNR.sTrajList, 1/EH_hSNR.dt, trajHighCutFreq);
    IT_lSNR.sTrajList = LPF(IT_lSNR.sTrajList, 1/IT_lSNR.dt, trajHighCutFreq);
    IT_hSNR.sTrajList = LPF(IT_hSNR.sTrajList, 1/IT_hSNR.dt, trajHighCutFreq);
    
    % Relative exploration effort - Angular distance traveled
    % crop to ignore the distance due to different initial 
    % position
    EH_hSNR.moleDist = cumDist(EH_hSNR.sTrajList(100:end));
    EH_lSNR.moleDist = cumDist(EH_lSNR.sTrajList(100:end));
    IT_hSNR.moleDist = cumDist(IT_hSNR.sTrajList(100:end));
    IT_lSNR.moleDist = cumDist(IT_lSNR.sTrajList(100:end));
    
    reErgSS(i) = EH_hSNR.moleDist;
    reErgWS(i) = EH_lSNR.moleDist;
    reInfSS(i) = IT_hSNR.moleDist;
    reInfWS(i) = IT_lSNR.moleDist;
end
reErgSS_Mean = mean(reErgSS);
reErgWS_Mean = mean(reErgWS);
reInfSS_Mean = mean(reInfSS);
reInfWS_Mean = mean(reInfWS);
reErgSS_SEM = 1.96 * std(reErgSS) / sqrt(length(reErgSS));
reErgWS_SEM = 1.96 * std(reErgWS) / sqrt(length(reErgWS));
reInfSS_SEM = 1.96 * std(reInfSS) / sqrt(length(reInfSS));
reInfWS_SEM = 1.96 * std(reInfWS) / sqrt(length(reInfWS));

fprintf(['---------------------------------------------------------\n',...
    '|         Relative exploration (mean +/- 95%% CI)        |\n', ...
    '|--------------------------------------------------------|',...
    '\n|> Ergodic Harvesting (strong signal) = %.3f +/- %.3f <|', ...
    '\n|> Ergodic Harvesting (weak signal)   = %.3f +/- %.3f <|',...
    '\n|>     Infotaxis (strong signal)      = %.3f +/- %.3f <|', ...
    '\n|>     Infotaxis (weak signal)        = %.3f +/- %.3f <|\n',...
    '----------------------------------------------------------\n'], ...
    reErgSS_Mean, reErgSS_SEM, reErgWS_Mean, reErgWS_SEM, ...
    reInfSS_Mean, reInfSS_SEM, reInfWS_Mean, reInfWS_SEM);

% Compute statistics
% Ranksum test, right tail = median of A is greater than median of B
[P, ~, Stats] = ranksum(re12, re4,...
    'tail', 'right');
fprintf('Behavior statistics: Wilcoxon rank sum test (one-sided) - p = %.4f (n = %d)\n', ...
    P, length([re12; re4]));
[P, ~, Stats] = ranksum(reErgWS, reErgSS,...
    'tail', 'right');
fprintf('EIH statistics: Wilcoxon rank sum test (one-sided) - p = %.4f (n = %d)\n', ...
    P, length([reErgWS, reErgSS]));

% Plot group data
axes; hold on;
hL = line([0,3], [1,1], 'LineStyle', '--', ...
    'Color', [140,140,140]/255.0, 'LineWidth',4);
hBoxPlot = notBoxPlot([re4', re12'], ...
    [1.1*ones(1,length(re4)), 1.9*ones(1,length(re12))], ...
    'jitter', 0.04, 'alpha', 0.5, 'DotSize', 24);
opt = [];
opt.BoxDim = [8,5]*0.354;
opt.YLabel = 'Relative Exploration'; % ylabel
opt.XLim = [0.5, 2.5];
opt.YLim = [0, 12];
opt.YTick = [1, 3:3:12];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 13;
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hL.Color = [140,140,140]/255.0;
hL.LineWidth = 4;
hL.LineStyle = '--';
hLine(1).LineWidth = 3;
hLine(2).LineWidth = 3;
hLine(1).Color = [162,0,0]/255.0;
hLine(2).Color = [50,180,74]/255.0;
legend(gca, 'off');
set(gca,'XTickLabel',{...
    '\color[rgb]{0.1961,0.7059,0.2902}Strong Signal',...
    '\color[rgb]{0.6353,0,0}Weak Signal'})
set(gca,'YTickLabel',{'1x', '3x', '6x', '9x', '12x'})
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.6, 0.6];
set(gca, 'Position', axesPosition);

% EIH and Infotaxis data
axes; hold on;
hBoxPlot = notBoxPlot([reErgSS, reInfSS, reErgWS, reInfWS], ...
    [0.5*ones(1,length(reErgSS)), ...
    0.75*ones(1,length(reInfSS)), ...
    1.25*ones(1,length(reErgWS)),...
    1.5*ones(1,length(reInfWS))], ...
    'jitter', 0.04, 'alpha', 0.5, 'DotSize', 24);
cmap = [...
    85, 1, 159; ...
    255, 196, 3;...
    85, 1, 159; ...
    255, 196, 3;...
    ] ./ 255;
opt = [];
opt.BoxDim = [8,5]*0.354;
opt.YLabel = 'Exploration'; % ylabel
opt.YLim = [0, 6];
opt.YTick = [0, 2, 4, 6];
opt.XLim = [0.2, 1.8];
opt.XTick = [0.625, 1.375];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 13;
opt.Colors = cmap;
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hLine(1).LineWidth = 3;
hLine(2).LineWidth = 3;
hLine(3).LineWidth = 3;
hLine(4).LineWidth = 3;
hLine(1).Color = cmap(1, :);
hLine(2).Color = cmap(2, :);
hLine(3).Color = cmap(1, :);
hLine(4).Color = cmap(2, :);
legend(gca, 'off');
set(gca,'XTickLabel',{...
    '\color[rgb]{0.1961,0.7059,0.2902}Strong Signal',...
    '\color[rgb]{0.6353,0,0}Weak Signal'})
set(gca,'YTickLabel',{'0', '2', '4', '6'})
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.6, 0.3];
set(gca, 'Position', axesPosition);


% All set, now print the first section into PDF
if SPLIT_PLOT
    splitprint(gcf,... %separate the current figure
        GEN_SAVE_PATH('fig3-Cockroach'),... %filenames will begin with 'disp2'
        {{'line';'text'},{'surface';'patch';'image'}}, ...% types of objects
        {'-dpdf','-dtiff'},... %file formats
        0,... %alignment mark will not be added
        [1 0],... %axes in first figure will be visible
        {'','-r400'});
else
    print(GEN_SAVE_PATH('fig3-Cockroach.pdf'),'-dpdf');
end

function mPlotContinuousEID(dat)
global PLOT_EER_BAND
if ~PLOT_EER_BAND
    return;
end
%% Plot Parameters
tScale = 10;   % Interval of EID plot update, set to 1 will plot all of the EID map
nBins = 256;   % Color resolution in the y axis
alpha = 0.5;  % Transparency of the EID color
cmap = [0.7 0 0.4];
eidList = dat.eidList;
tRes = length(dat.oTrajList) / (size(eidList,2)-1);
sRes = size(eidList,1);
s = 1 / sRes;
faces = 1:4;

idxList = tScale:tScale:floor(length(dat.oTrajList) / tRes);
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
        patch('Faces',faces,'Vertices',verts,...
            'FaceColor', cmap,...
            'FaceAlpha', alpha*bin(k)/nBins,...
            'EdgeColor', 'none');
    end
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

function re = calcRE(dat)
    % Computes relative exploration
    len = size(dat, 1);
    re = zeros(len, 1, 'double');
    for i = 1:len
        re(i) = dat.walk_dist(i) / dat.distRef(i);
    end

function dist = calcSumLength2D(path)
    %% Calculates the cumulative 2D euclidean distance
    distMat = pdist(path, 'euclidean');
    distMat = squareform(distMat);

    dist = 0;
    for i = 2:size(path,1)
        dist = dist + distMat(i, i-1);
    end


function moleData = processMoleTraj(moleData)
    findAngle = @(a, b) acosd(min(1,max(-1, a(:).' * b(:) / norm(a) / norm(b) )));
    refVec = moleData.refPath(2, :) - moleData.refPath(1, :);
    refAngle = findAngle(refVec/norm(refVec), [0, 1]);
    if refVec(1) > 0
        refAngle = -refAngle;
    end
    moleData.rotPath = rotate2D(moleData.molePath, refAngle);
    moleData.rotRef = rotate2D(moleData.refPath, refAngle);

function x = rotate2D(x, theta)
    offset = x(1, :);
    x = x - offset;
    rot2x = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    for i = 1:size(x, 1)
        x(i, :) = x(i, :) * rot2x;
    end
    x = x + offset;