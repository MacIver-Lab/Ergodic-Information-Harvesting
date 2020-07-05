function makeFigure3Plot(dataPath, savePath)
%% Plot figure 3

warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, 'behavior_sim', fname);
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile('./PublishedData/', 'animal_behavior_data', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
barColor = [72, 110, 181;...
    50, 180, 74; ...
    236, 29, 36] / 255;

% Whether nor not to plot EER band
% set to 1 to enable EER plot overlay
% note that due to the complexity of EER bands, it's can be fairly slow 
% to plot the EER bands
global PLOT_EER_BAND
PLOT_EER_BAND = 1;

% Lock random number seed.
rng(0);

%% Electric Fish Simulation
% Load data
EH_lSNR = load(GEN_DATA_PATH('EIH-ElectricFish-WeakSignal-RandSeed-1.mat'), ...
    'oTrajList', 'sTrajList', 'dt', 'phi');
EH_hSNR = load(GEN_DATA_PATH('EIH-ElectricFish-StrongSignal-RandSeed-1.mat'), ...
    'oTrajList', 'sTrajList', 'dt', 'phi');
EH_lSNR.eidList = flattenResultList(EH_lSNR.phi(:,:,1:end-1))';
EH_hSNR.eidList = flattenResultList(EH_hSNR.phi(:,:,1:end-1))';
% Filter trajectory
simTrajHighCutFreq = 2.10;
EH_lSNR.sTrajList = LPF(EH_lSNR.sTrajList, 1/EH_lSNR.dt, simTrajHighCutFreq);
EH_hSNR.sTrajList = LPF(EH_hSNR.sTrajList, 1/EH_hSNR.dt, simTrajHighCutFreq);

% Load fish behavioral data
medFreqRange = [0.25, 1];
cmap = [56, 180, 74; ...
    161, 30, 34] / 255;
barColors = [[75, 111, 182]/255; cmap];
rng(0);
% Load behavior data
fish_HSNR = load(GEN_BEHAVIOR_DATA_PATH('ElectricFish/Eigenmannia-sp-StrongSignal.mat'), 'strongSigData');
fish_HSNR = fish_HSNR.strongSigData;
fish_LSNR = load(GEN_BEHAVIOR_DATA_PATH('ElectricFish/Eigenmannia-sp-WeakSignal.mat'), 'weakSigData');
fish_LSNR = fish_LSNR.weakSigData;
fish.hSNR.fishTraj = fish_HSNR{8, 1};
fish.hSNR.refugeTraj = fish_LSNR{8, 2};
FPS = 60;
dt = 1 / FPS;
trajSegStr = 31.5;
trajSegLen = 30;
trajSegStrIdx = trajSegStr * FPS;
trajSegIdx = (trajSegStr+trajSegLen) * FPS;
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

% calibrate distance unit
pix2cm = 1.7 / (max(fish.hSNR.refugeTraj) - min(fish.hSNR.refugeTraj));

for i = 1:length(fish_HSNR)
    % loop through strong signal trials
    fishTraj = fish_HSNR{i, 1};
    % FFT
    [freqTicks, fishMagHSNR(i, :)] = decomposeFourierMag(fishTraj, [], 1/60);
end
for i = 1:length(fish_LSNR)
    % loop through weak signal trials
    fishTraj = fish_LSNR{i, 1};
    % FFT
    [freqTicks, fishMagLSNR(i, :)] = decomposeFourierMag(fishTraj, [], 1/60);
end
% Compute averaged medium frequency range magnitude
freqIdxStart = find(freqTicks <= medFreqRange(1), 1,'last'); % starting index
freqIdxEnd = find(freqTicks >= medFreqRange(2), 1,'first');  % start index
freqIdxRange = freqIdxStart:freqIdxEnd;
medMagHSNR = mean(fishMagHSNR(:, freqIdxRange), 2);
medMagLSNR = mean(fishMagLSNR(:, freqIdxRange), 2);

magList = [medMagHSNR', medMagLSNR'];
magList = (magList - mean(medMagHSNR)) / max(magList);
magLabel = [ones(1,length(medMagHSNR)), 2*ones(1,length(medMagLSNR))];
% Behavior spectrum statistics
figure(1); clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [16 8]);
notBoxPlot(magList, magLabel, 'jitter', 0.05, 'alpha', 0.5);
ylim([-0.3, 1])
xlim([0.5, 2.5]);
set(gca, 'YTick', []);
set(gca, 'XTick', [1, 2]);
set(gca, 'XTickLabel', ...
    {'\color[rgb]{0.1961,0.7059,0.2902}Strong Signal', ...
    '\color[rgb]{0.6353,0,0}Weak Signal'});
set(gca, 'FontName', 'Helvetica');
ylabel('Fourier Magnitude');
hLine = findobj(gca,'Type','line');
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.YTick = [0, 1];
opt.FontName = 'Helvetica';
opt.FontSize = 12;
setAxesProp(opt);
legend('off');
hLine(1).LineWidth = 4;
hLine(1).Color = cmap(2, :);
hLine(2).LineWidth = 4;
hLine(2).Color = cmap(1, :);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.6, 0.6];
set(gca, 'Position', axesPosition);
title('Electric Fish (Behavior)', 'FontSize', 12)
% Test statistics (non-parametric ANOVA)
p = kruskalwallis(magList, magLabel, 'off');

% Behavior FFT
axes; hold on;
[~, refugeMag] = decomposeFourierMag(fish_HSNR{1, 2}, [], 1/60);
nrmFactor = max([fishMagHSNR(2, :),fishMagLSNR(1, :),refugeMag']);
plot(freqTicks, refugeMag/nrmFactor, ...
        'LineWidth', 2, 'Color', barColors(1, :));
plot(freqTicks, fishMagHSNR(2, :)/nrmFactor, ...
        'LineWidth', 2, 'Color', cmap(1, :));
plot(freqTicks, fishMagLSNR(1, :)/nrmFactor, ...
    'LineWidth', 2, 'Color', cmap(2, :));
ylim([0, 1]);
patch(...
    'Vertices', [0.25, 0; 1, 0; 1, 1.7e5; 0.25, 1.7e5], ...
    'Faces', [1, 2, 3, 4], ...
    'FaceColor', 'k', 'FaceAlpha', 0.1, ...
    'EdgeAlpha', 0)
xlabel('Frequency');
ylabel('Fourier Magnitude');
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [0, 0.25, 1];
opt.YTick = [0, 1];
opt.XLim = [0.0, 1];
opt.FontSize = 12;
opt.FontName = 'Helvetica';
opt.Colors = barColors;
setAxesProp(opt);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.4, 0.6];
set(gca, 'Position', axesPosition);
title('Electric Fish (Behavior)', 'FontSize', 12)

% Load simulation data
eihFilesHSNR = dir(GEN_DATA_PATH('EIH-ElectricFish-StrongSignal-RandSeed-*.mat'));
eihFilesLSNR = dir(GEN_DATA_PATH('EIH-ElectricFish-WeakSignal-RandSeed-*.mat'));
for i = 1:length(eihFilesHSNR)
    dat = load(GEN_DATA_PATH(eihFilesHSNR(i).name), 'oTrajList', 'sTrajList', 'dt');
    % loop through strong signal trials
    sensorTraj = dat.sTrajList;
    % FFT
    [freqTicks, sensorMagHSNR(i, :)] = decomposeFourierMag(sensorTraj, [], dat.dt/2);
end
for i = 1:length(eihFilesLSNR)
    dat = load(GEN_DATA_PATH(eihFilesLSNR(i).name), 'oTrajList', 'sTrajList', 'dt');
    % loop through strong signal trials
    sensorTraj = dat.sTrajList;
    % FFT
    [freqTicks, sensorMagLSNR(i, :)] = decomposeFourierMag(sensorTraj, [], dat.dt/2);
end
% Compute averaged medium frequency range magnitude
medMagHSNR_EIH = mean(sensorMagHSNR(:, freqIdxRange), 2);
medMagLSNR_EIH = mean(sensorMagLSNR(:, freqIdxRange), 2);

magList = [medMagHSNR_EIH', medMagLSNR_EIH'];
magList = (magList - mean(medMagHSNR_EIH)) / max(magList);
magLabel = [ones(1,length(medMagHSNR_EIH)), 2*ones(1,length(medMagLSNR_EIH))];

axes; hold on;
notBoxPlot(magList, magLabel, 'jitter', 0.05, 'alpha', 0.5);
ylim([-0.3, 1])
xlim([0.5, 2.5]);
set(gca, 'YTick', []);
set(gca, 'XTick', [1, 2]);
set(gca, 'XTickLabel', ...
    {'\color[rgb]{0.1961,0.7059,0.2902}Strong Signal', ...
    '\color[rgb]{0.6353,0,0}Weak Signal'});
set(gca, 'FontName', 'Helvetica');
ylabel('Fourier Magnitude');
hLine = findobj(gca,'Type','line');
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.YTick = [0, 1];
opt.FontName = 'Helvetica';
opt.FontSize = 12;
setAxesProp(opt);
legend('off');
hLine(1).LineWidth = 4;
hLine(1).Color = cmap(2, :);
hLine(2).LineWidth = 4;
hLine(2).Color = cmap(1, :);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.6, 0.3];
set(gca, 'Position', axesPosition);
title('Electric Fish (Simulation)', 'FontSize', 12)
% Test statistics (non-parametric ANOVA)
p = kruskalwallis(magList, magLabel, 'off');

% Plot FFT
axes; hold on;
dat = load(GEN_DATA_PATH(eihFilesHSNR(1).name), 'oTrajList', 'dt');
[~, trajMag] = decomposeFourierMag(dat.oTrajList, [], dat.dt/2);
nrmFactor = max([sensorMagHSNR(3, :),sensorMagLSNR(3, :),trajMag]);
plot(freqTicks, trajMag/nrmFactor, ...
        'LineWidth', 2, 'Color', [75, 111, 182]/255);
plot(freqTicks, sensorMagHSNR(3, :)/nrmFactor, ...
        'LineWidth', 2, 'Color', cmap(1, :));
plot(freqTicks, sensorMagLSNR(3, :)/nrmFactor, ...
    'LineWidth', 2, 'Color', cmap(2, :));
xlabel('Frequency');
ylabel('Fourier Magnitude');
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [0, 0.25, 1];
opt.YTick = [0, 1];
opt.XLim = [0, 1];
opt.YLim = [0, 1];
opt.FontSize = 12;
opt.FontName = 'Helvetica';
opt.Colors = barColors;
setAxesProp(opt, gca);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.4, 0.3];
set(gca, 'Position', axesPosition);
title('Electric Fish (Simulation)', 'FontSize', 12)
patch(...
    'Vertices', [0.25, 0; 1, 0; 1, 250; 0.25, 250], ...
    'Faces', [1, 2, 3, 4], ...
    'FaceColor', 'k', 'FaceAlpha', 0.1, ...
    'EdgeAlpha', 0)

% Statistics plot
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
[P, ~, ~] = ranksum(reWeakSignalFish, reStrongSignalFish,...
    'tail', 'right');
fprintf('Behavior statistics: Wilcoxon rank sum test (one-sided) - p = %.4f (n = %d)\n', ...
    P, length([reWeakSignalFish; reStrongSignalFish]));
[P, ~, ~] = ranksum(reErgWS, reErgSS,...
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
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.6];
set(gca, 'Position', axesPosition);
title('Electric Fish (Behavior)', 'FontSize', 12)

% EIH
axes; hold on;
hL = line([-1,5], [1,1], 'LineStyle', '--', ...
    'Color', [140,140,140]/255.0, 'LineWidth',4);
hBoxPlot = notBoxPlot([reErgSS, reErgWS], ...
    [1.1*ones(1,length(reErgSS)), 1.9*ones(1,length(reErgWS))], ...
    'jitter', 0.04, 'alpha', 0.5, 'dotSize', 24);
cmap = [85, 1, 159; ...
    85, 1, 159;] ./ 255;
opt = [];
opt.BoxDim = [8,5]*0.354;
opt.YLabel = 'Relative Exploration'; % ylabel
opt.YLim = [0.75, 3.5];
opt.YTick = [1, 2, 3];
opt.XLim = [0.5, 2.5];
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
hLine(1).Color = [162,0,0]/255.0;
hLine(2).Color = [50,180,74]/255.0;
legend(gca, 'off');
set(gca,'XTickLabel',{...
    '\color[rgb]{0.1961,0.7059,0.2902}Strong Signal',...
    '\color[rgb]{0.6353,0,0}Weak Signal'})
set(gca,'YTickLabel',{'1x', '2x', '3x'})
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.3];
set(gca, 'Position', axesPosition);
title('Electric Fish (Simulation)', 'FontSize', 12)

% All set, now print the first section into PDF
print(GEN_SAVE_PATH('fig3.pdf'),'-dpdf');

function mPlotContinuousEID(dat)
global PLOT_EER_BAND
if ~PLOT_EER_BAND
    return;
end
%% Plot Parameters
tScale = 5;   % Interval of EID plot update, set to 1 will plot all of the EID map
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

function [freqTicks, gainX] = decomposeFourierMag(x, highCut, varargin)
dFreq = 0.005;
maxFreq = 1;
if length(varargin) == 1
    dt = varargin{1};
elseif length(varargin) == 2
    dt = varargin{1};
    maxFreq = varargin{2};
elseif length(varargin) == 3
    dt = varargin{1};
    minFreq = varargin{2};
    maxFreq = varargin{3};
end
fftEndIdx = floor(maxFreq/dFreq) + 1;
freqTicks = dFreq:dFreq:maxFreq;
if isempty(highCut)
    highCutFreq = 2.10;
else
    highCutFreq = highCut;
end
nFreqSamps = @(t) ceil((dFreq*t)^-1);
% apply zero-phase low-pass filter to traj
fs = 1 / dt;
x = LPF(x, fs, highCutFreq);
% FFT
fftX = fft(x-mean(x), nFreqSamps(dt));
gainX = abs(fftX(2:fftEndIdx));

function locs = findLocalPeaks(x, refLocs)
winWidth = 3;
locs = zeros(size(refLocs));
for i = 1:length(refLocs)
    searchWin = (refLocs(i)-winWidth):(refLocs(i)+winWidth);
    [~, idx] = max(x(searchWin));
    locs(i) = searchWin(idx);
end