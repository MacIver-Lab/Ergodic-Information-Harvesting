function makeFig3Plot(dataPath, savePath)

warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, fname);
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile(pwd, 'FigureCode', 'fig3', 'BehaviorData', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);

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

% Weak Signal, n = 12
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

% Compute statistics
% Ranksum test, right tail = median of rteWeakSignal is greater than median of rteStrongSignal
[P, ~] = ranksum(reWeakSignalFish, reStrongSignalFish,...
    'tail', 'right');
fprintf('Wilcoxon rank sum test (one-sided) - p = %.4f\n', P);

% Ergodic Harvesting Data
% Load Ergodic data
EH_lSNR = load(GEN_DATA_PATH('ErgodicHarvesting-WeakSignal-Sine.mat'), ...
    'oTrajList', 'sTrajList', 'dt');
EH_hSNR = load(['./Data/ElectricFish/', ...
    'ErgodicHarvesting-StrongSignal-Sine.mat'], ...
    'oTrajList', 'sTrajList', 'dt');

% Load Infotaxis data
IT_lSNR = load(GEN_DATA_PATH('fig3-Infotaxis-ElectricFish-SNR-30.mat'), ...
    'oTrajList', 'sTrajList', 'dt');
IT_hSNR = load(GEN_DATA_PATH('fig3-Infotaxis-ElectricFish-SNR-60.mat'), ...
    'oTrajList', 'sTrajList', 'dt');

% Filter trajectory
trajHighCutFreq = 2.10; 
EH_lSNR.sTrajList = LPF(EH_lSNR.sTrajList, 1/EH_lSNR.dt, trajHighCutFreq);
EH_hSNR.sTrajList = LPF(EH_hSNR.sTrajList, 1/EH_hSNR.dt, trajHighCutFreq);
IT_lSNR.sTrajList = LPF(IT_lSNR.sTrajList, 1/IT_lSNR.dt, trajHighCutFreq);
IT_hSNR.sTrajList = LPF(IT_hSNR.sTrajList, 1/IT_hSNR.dt, trajHighCutFreq);

% Cumulative 1D distance traveled with the initial condition cropped to 
% ensure consistency
EH_hSNR.sPower = cumDist(EH_hSNR.sTrajList(200:end));
EH_hSNR.oPower = cumDist(EH_hSNR.oTrajList(200:end));
EH_lSNR.sPower = cumDist(EH_lSNR.sTrajList(200:end));
EH_lSNR.oPower = cumDist(EH_lSNR.oTrajList(200:end));
IT_hSNR.sPower = cumDist(IT_hSNR.sTrajList(200:end));
IT_hSNR.oPower = cumDist(IT_hSNR.oTrajList(200:end));
IT_lSNR.sPower = cumDist(IT_lSNR.sTrajList(200:end));
IT_lSNR.oPower = cumDist(IT_lSNR.oTrajList(200:end));

reStrongSignalErg = EH_hSNR.sPower/EH_hSNR.oPower;
reWeakSignalErg = EH_lSNR.sPower/EH_lSNR.oPower;
reStrongSignalInfo = IT_hSNR.sPower/IT_hSNR.oPower;
reWeakSignalInfo = IT_lSNR.sPower/IT_lSNR.oPower;
fprintf(['\t\t\tRelative exploration\n', ...
    '---------------------------------------------',...
    '\n\t> Ergodic Harvesting (strong signal) = %.3f', ...
    '\n\t> Ergodic Harvesting (weak signal) = %.3f\n',...
    '\n\t> Infotaxis (strong signal) = %.3f', ...
    '\n\t> Infotaxis (weak signal) = %.3f\n'], ...
    reStrongSignalErg, reWeakSignalErg, reStrongSignalInfo, reWeakSignalInfo);

% Plot group data
figure(1); clf; hold on;
hBoxPlot = notBoxPlot([reStrongSignalFish',reWeakSignalFish'], ...
    [ones(1,length(reStrongSignalFish)), 2*ones(1,length(reWeakSignalFish))]);
opt = [];
opt.BoxDim = [8,5];
opt.YLabel = 'Relative Exploration'; % ylabel
opt.YLim = [0.75, 3.5];
opt.YTick = [1, 2, 3];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hLine(1).LineStyle = 'none';
hLine(3).LineStyle = 'none';
hLine(1).Marker = '.';
hLine(3).Marker = '.';
hLine(1).MarkerEdgeColor = 'b';
hLine(1).MarkerSize = 26;
hLine(3).MarkerSize = 26;
hLine = line([0,3],[1,1]);
hLine.LineStyle = '--';
hLine.LineWidth = 2;
legend(gca, 'off');
set(gca,'XTickLabel',{'Strong Signal','Weak Signal'})
set(gca,'YTickLabel',{'1x', '2x', '3x'})
print(gcf,'-dpdf','fish-RE.pdf');


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

% Compute statistics
% Ranksum test, right tail = median of reBlk is greater than 
% median of reNrm
[P, ~, ~] = ranksum(reBlk, reNrm,...
    'tail', 'right');
fprintf('Wilcoxon rank sum test (one-sided) - p = %.4f\n', P);


% Mole Odor Localization
% Ergodic Harvesting
EH_hSNR = load(GEN_DATA_PATH('ErgodicHarvesting-StrongSignal.mat'), ...
    'sTrajList', 'sTrajList', 'dt');
EH_lSNR = load(GEN_DATA_PATH('ErgodicHarvesting-WeakSignal.mat'), ...
    'sTrajList', 'sTrajList', 'dt');

IT_hSNR = load(GEN_DATA_PATH('fig3-Infotaxis-Mole-SNR-60.mat'), ...
    'sTrajList', 'sTrajList', 'dt');
IT_lSNR = load(GEN_DATA_PATH('fig3-Infotaxis-Mole-SNR-10.mat'), ...
    'sTrajList', 'sTrajList', 'dt');

% Filter trajectory - with moving average filter
EH_lSNR.sTrajList = movmean(EH_lSNR.sTrajList, 2);
EH_hSNR.sTrajList = movmean(EH_hSNR.sTrajList, 2);
IT_lSNR.sTrajList = movmean(IT_lSNR.sTrajList, 2);
IT_hSNR.sTrajList = movmean(IT_hSNR.sTrajList, 2);

% Relative exploration effort - Angular distance traveled
EH_hSNR.moleDist = cumAngularDist(EH_hSNR.sTrajList);
EH_lSNR.moleDist = cumAngularDist(EH_lSNR.sTrajList);
IT_hSNR.moleDist = cumAngularDist(IT_hSNR.sTrajList);
IT_lSNR.moleDist = cumAngularDist(IT_lSNR.sTrajList);

fprintf(['\t\t\tRelative exploration\n', ...
    '---------------------------------------------',...
    '\n\t> Ergodic Harvesting (strong signal) = %.3f', ...
    '\n\t> Ergodic Harvesting (weak signal) = %.3f\n',...
    '\n\t> Infotaxis (strong signal) = %.3f', ...
    '\n\t> Infotaxis (weak signal) = %.3f\n'], ...
    EH_hSNR.moleDist, EH_lSNR.moleDist, IT_hSNR.moleDist, IT_lSNR.moleDist);

% Plot group data
figure(2); clf;
hBoxPlot = notBoxPlot([reBlk,reNrm], ...
    [2*ones(1,length(reBlk)), 1*ones(1,length(reNrm))]);
opt = [];
opt.BoxDim = [8,5];
opt.YLabel = 'Relative Exploration'; % ylabel
opt.YLim = [1, 7];
opt.YTick = [1, 2, 4, 6];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
setAxesProp(opt);
hLine = findobj(gca,'Type','line');
hLine(1).LineStyle = 'none';
hLine(3).LineStyle = 'none';
hLine(1).Marker = '.';
hLine(3).Marker = '.';
hLine(1).MarkerEdgeColor = 'b';
hLine(1).MarkerSize = 26;
hLine(3).MarkerSize = 26;
% hold on;
% hLine = line([0,3],[1,1]);
% hLine.LineStyle = '--';
% hLine.LineWidth = 2;
legend(gca, 'off');
set(gca,'XTickLabel',{'Strong Signal','Weak Signal'})
set(gca,'YTickLabel',{'1x', '2x', '4x', '6x'})
print(gcf,'-dpdf','mole-RE.pdf');


function dist = calcSumLength2D(path)
%% Calculates the cumulative 2D euclidean distance of a given 2D curve
distMat = pdist(path, 'euclidean');
distMat = squareform(distMat);

dist = 0;
for i = 2:size(path,1)
    dist = dist + distMat(i, i-1);
end