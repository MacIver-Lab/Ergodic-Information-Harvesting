function makeFig4Plot(dataPath, savePath)

close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, 'behavior_sim', fname);
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile(dataPath, 'animal_behavior_data', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
% Lock random number generator seed
rng(0);

%% Electric Fish Behavioral Trials
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
figure(1); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
hL = line([0,3], [1,1], 'LineStyle', '--', ...
    'Color', [140,140,140]/255.0, 'LineWidth',4);
hBoxPlot = notBoxPlot([reStrongSignalFish',reWeakSignalFish'], ...
    [1.1*ones(1,length(reStrongSignalFish)), 1.9*ones(1,length(reWeakSignalFish))], ...
    'jitter', 0.04, 'alpha', 0.5);
opt = [];
opt.BoxDim = [8,5]*0.6;
opt.YLabel = 'Relative Exploration'; % ylabel
opt.XLim = [0.5, 2.5];
opt.YLim = [0.75, 3.5];
opt.YTick = [1, 2, 3];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hL.Color = [140,140,140]/255.0;
hL.LineWidth = 4;
hL.LineStyle = '--';
hLine(1).LineWidth = 5;
hLine(2).LineWidth = 5;
hLine(1).Color = [162,0,0]/255.0;
hLine(2).Color = [50,180,74]/255.0;
legend(gca, 'off');
set(gca,'XTickLabel',{...
    '\color[rgb]{0.1961,0.7059,0.2902}Strong Signal',...
    '\color[rgb]{0.6353,0,0}Weak Signal'})
set(gca,'YTickLabel',{'1x', '2x', '3x'})
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.4];
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
    'jitter', 0.04, 'alpha', 0.5);
cmap = [...
    85, 1, 159; ...
    255, 196, 3;...
    85, 1, 159; ...
    255, 196, 3;...
    ] ./ 255;
opt = [];
opt.BoxDim = [8,5]*0.6;
opt.YLabel = 'Relative Exploration'; % ylabel
opt.YLim = [0.75, 3.5];
opt.YTick = [1, 2, 3];
opt.XLim = [0,2];
opt.XTick = [0.625, 1.375];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.Colors = cmap;
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hL.Color = [140,140,140]/255.0;
hL.LineWidth = 4;
hL.LineStyle = '--';
hLine(1).LineWidth = 5;
hLine(2).LineWidth = 5;
hLine(3).LineWidth = 5;
hLine(4).LineWidth = 5;
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
axesPosition(1:2) = [0.5, 0.4];
set(gca, 'Position', axesPosition);
print(gcf,'-dpdf',GEN_SAVE_PATH('fig4AB-electric-fish.pdf'));

%% Mole behavioral data
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
figure(3); clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
hBoxPlot = notBoxPlot([reBlk,reNrm], ...
    [1.9*ones(1,length(reBlk)), 1.1*ones(1,length(reNrm))], ...
    'jitter', 0.04, 'alpha', 0.5);
opt = [];
opt.BoxDim = [8,5]*0.6;
opt.YLabel = 'Relative Exploration'; % ylabel
opt.XLim = [0.5, 2.5];
opt.YLim = [1, 6];
opt.YTick = [1, 2, 4, 6];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hLine(1).LineWidth = 5;
hLine(2).LineWidth = 5;
hLine(1).Color = [162,0,0]/255.0;
hLine(2).Color = [50,180,74]/255.0;
legend(gca, 'off');
set(gca,'XTickLabel',{...
    '\color[rgb]{0.1961,0.7059,0.2902}Strong Signal',...
    '\color[rgb]{0.6353,0,0}Weak Signal'})
set(gca,'YTickLabel',{'1x', '2x', '4x', '6x'})
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.4];
set(gca, 'Position', axesPosition);

% EIH and Infotaxis data
axes; hold on;
% hL = line([-1,5], [1,1], 'LineStyle', '--', ...
%     'Color', [140,140,140]/255.0, 'LineWidth',4);
hBoxPlot = notBoxPlot([reErgSS, reInfSS, reErgWS, reInfWS], ...
    [0.5*ones(1,length(reErgSS)), ...
    0.75*ones(1,length(reInfSS)), ...
    1.25*ones(1,length(reErgWS)),...
    1.5*ones(1,length(reInfWS))], ...
    'jitter', 0.04, 'alpha', 0.5);
cmap = [...
    85, 1, 159; ...
    255, 196, 3;...
    85, 1, 159; ...
    255, 196, 3;...
    ] ./ 255;
opt = [];
opt.BoxDim = [8,5]*0.6;
opt.YLabel = 'Exploration'; % ylabel
opt.YLim = [0, 6];
opt.YTick = [0, 2, 4, 6];
opt.XLim = [0.2, 1.8];
opt.XTick = [0.625, 1.375];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.Colors = cmap;
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
% hL.Color = [140,140,140]/255.0;
% hL.LineWidth = 4;
% hL.LineStyle = '--';
hLine(1).LineWidth = 5;
hLine(2).LineWidth = 5;
hLine(3).LineWidth = 5;
hLine(4).LineWidth = 5;
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
axesPosition(1:2) = [0.5, 0.4];
set(gca, 'Position', axesPosition);
print(gcf,'-dpdf',GEN_SAVE_PATH('fig4CD-mole.pdf'));

%% Cockroach
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
figure(5); clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
hL = line([0,3], [1,1], 'LineStyle', '--', ...
    'Color', [140,140,140]/255.0, 'LineWidth',4);
hBoxPlot = notBoxPlot([re4', re12'], ...
    [1.1*ones(1,length(re4)), 1.9*ones(1,length(re12))]);
opt = [];
opt.BoxDim = [8,5]*0.6;
opt.YLabel = 'Relative Exploration'; % ylabel
opt.XLim = [0.5, 2.5];
opt.YLim = [0, 12];
opt.YTick = [1, 3:3:12];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hL.Color = [140,140,140]/255.0;
hL.LineWidth = 4;
hL.LineStyle = '--';
hLine(1).LineWidth = 5;
hLine(2).LineWidth = 5;
hLine(1).Color = [162,0,0]/255.0;
hLine(2).Color = [50,180,74]/255.0;
legend(gca, 'off');
set(gca,'XTickLabel',{...
    '\color[rgb]{0.1961,0.7059,0.2902}Strong Signal',...
    '\color[rgb]{0.6353,0,0}Weak Signal'})
set(gca,'YTickLabel',{'1x', '3x', '6x', '9x', '12x'})
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.4];
set(gca, 'Position', axesPosition);

% EIH and Infotaxis data
axes; hold on;
% hL = line([-1,5], [1,1], 'LineStyle', '--', ...
%     'Color', [140,140,140]/255.0, 'LineWidth',4);
hBoxPlot = notBoxPlot([reErgSS, reInfSS, reErgWS, reInfWS], ...
    [0.5*ones(1,length(reErgSS)), ...
    0.75*ones(1,length(reInfSS)), ...
    1.25*ones(1,length(reErgWS)),...
    1.5*ones(1,length(reInfWS))]);
cmap = [...
    85, 1, 159; ...
    255, 196, 3;...
    85, 1, 159; ...
    255, 196, 3;...
    ] ./ 255;
opt = [];
opt.BoxDim = [8,5]*0.6;
opt.YLabel = 'Exploration'; % ylabel
opt.YLim = [0, 6];
opt.YTick = [0, 2, 4, 6];
opt.XLim = [0.2, 1.8];
opt.XTick = [0.625, 1.375];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.Colors = cmap;
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
% hL.Color = [140,140,140]/255.0;
% hL.LineWidth = 4;
% hL.LineStyle = '--';
hLine(1).LineWidth = 5;
hLine(2).LineWidth = 5;
hLine(3).LineWidth = 5;
hLine(4).LineWidth = 5;
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
axesPosition(1:2) = [0.5, 0.4];
set(gca, 'Position', axesPosition);
print(gcf,'-dpdf',GEN_SAVE_PATH('fig4EF-cockroach.pdf'));


function dist = calcSumLength2D(path)
%% Calculates the cumulative 2D euclidean distance of a given 2D curve
distMat = pdist(path, 'euclidean');
distMat = squareform(distMat);

dist = 0;
for i = 2:size(path,1)
    dist = dist + distMat(i, i-1);
end

function re = calcRE(dat)
% Computes relative exploration
len = size(dat, 1);
re = zeros(len, 1, 'double');
for i = 1:len
    re(i) = dat.walk_dist(i) / dat.distRef(i);
end