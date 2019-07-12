function makeFig5Plot(dataPath, savePath)
%% Distance from Ergodicity vs Mean Tracking Performance under different wiggle attenuation (fig 5)
% Chen Chen

close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, 'wiggle_attenuation_sim', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
% Lock random seed
rng(0);

%% Load Dataset
load(GEN_DATA_PATH('PerfData.mat'));
targetTrajAmp = 0.2; % Amplitude of the target trajectory
gAtten = PerfData.gAtten;
% meanEstErrorMeanBelief 
% Unit: percentage tracking error with respect to target trajectory
%       amplitude
meanEstErrorMeanBelief = 100.0 * PerfData.meanEstErrorMeanBelief / targetTrajAmp;
meanErgodicity = PerfData.meanErgodicity;
rmsErgodicity = PerfData.rmsErgodicity;
clear PerfData

%% Categorize data fields according to wiggle attenuation gain
AttenConditions = unique(gAtten);
nAttenConditions = length(AttenConditions);

% Initialize data field
attenDataSet(nAttenConditions).gAtten = [];
attenDataSet(nAttenConditions).meanEstError = [];
attenDataSet(nAttenConditions).meanErgodicity = [];

for idx = 1:nAttenConditions
    attenDataSet(idx).gAtten = AttenConditions(idx);
    attenDataSet(idx).meanEstError = meanEstErrorMeanBelief(gAtten == AttenConditions(idx));
    attenDataSet(idx).meanErgodicity = meanErgodicity(gAtten == AttenConditions(idx));
    attenDataSet(idx).rmsErgodicity = rmsErgodicity(gAtten == AttenConditions(idx));
end

%% Statistics
% Threshold for low, medium, and high attenuation trials
%     Low    = [ 5 -  15] dB
%     Medium = (15 -  80] dB
%     High   = (80 - 150] dB
AttenThreshLow = 15;
AttenThreshHigh = 80;
% Exclude the reference
AttenConditions(1) = NaN;
refTrial.meanEstErrRaw = [attenDataSet(1).meanEstError];
refTrial.meanErgodicityRaw = [attenDataSet(1).meanErgodicity];
refTrial.meanEstErrMean = mean([attenDataSet(1).meanEstError]);
refTrial.meanEstErrVar = var([attenDataSet(1).meanEstError]);
refTrial.meanEstErrStd = std([attenDataSet(1).meanEstError]);
refTrial.meanErgodicityMean = mean([attenDataSet(1).meanErgodicity]);
refTrial.meanErgodicityVar = var([attenDataSet(1).meanErgodicity]);
refTrial.meanErgodicityStd = std([attenDataSet(1).meanErgodicity]);

% Low Attenuation Trials
% Statistics
lowAttenTrials.meanEstErrRaw = [attenDataSet(AttenConditions<=AttenThreshLow).meanEstError];
lowAttenTrials.meanErgodicityRaw = [attenDataSet(AttenConditions<=AttenThreshLow).meanErgodicity];
lowAttenTrials.meanEstErrMean = mean([attenDataSet(AttenConditions<=AttenThreshLow).meanEstError]);
lowAttenTrials.meanEstErrVar = var([attenDataSet(AttenConditions<=AttenThreshLow).meanEstError]);
lowAttenTrials.meanEstErrStd = std([attenDataSet(AttenConditions<=AttenThreshLow).meanEstError]);
lowAttenTrials.meanErgodicityMean = mean([attenDataSet(AttenConditions<=AttenThreshLow).meanErgodicity]);
lowAttenTrials.meanErgodicityVar = var([attenDataSet(AttenConditions<=AttenThreshLow).meanErgodicity]);
lowAttenTrials.meanErgodicityStd = std([attenDataSet(AttenConditions<=AttenThreshLow).meanErgodicity]);

% Medium Attenuation Trials
% Statistics
medAttenTrials.meanEstErrRaw = [attenDataSet(AttenConditions>AttenThreshLow & AttenConditions<=AttenThreshHigh).meanEstError];
medAttenTrials.meanErgodicityRaw = [attenDataSet(AttenConditions>AttenThreshLow & AttenConditions<=AttenThreshHigh).meanErgodicity];
medAttenTrials.meanEstErrMean = mean([attenDataSet(AttenConditions>AttenThreshLow & AttenConditions<=AttenThreshHigh).meanEstError]);
medAttenTrials.meanEstErrVar = var([attenDataSet(AttenConditions>AttenThreshLow & AttenConditions<=AttenThreshHigh).meanEstError]);
medAttenTrials.meanEstErrStd = std([attenDataSet(AttenConditions>AttenThreshLow & AttenConditions<=AttenThreshHigh).meanEstError]);
medAttenTrials.meanErgodicityMean = mean([attenDataSet(AttenConditions>AttenThreshLow & AttenConditions<=AttenThreshHigh).meanErgodicity]);
medAttenTrials.meanErgodicityVar = var([attenDataSet(AttenConditions>AttenThreshLow & AttenConditions<=AttenThreshHigh).meanErgodicity]);
medAttenTrials.meanErgodicityStd = std([attenDataSet(AttenConditions>AttenThreshLow & AttenConditions<=AttenThreshHigh).meanErgodicity]);


% High Attenuation Trials
% Statistics
highAttenTrials.meanEstErrRaw = [attenDataSet(AttenConditions>AttenThreshHigh).meanEstError];
highAttenTrials.meanErgodicityRaw = [attenDataSet(AttenConditions>AttenThreshHigh).meanErgodicity];
highAttenTrials.meanEstErrMean = mean([attenDataSet(AttenConditions>AttenThreshHigh).meanEstError]);
highAttenTrials.meanEstErrVar = var([attenDataSet(AttenConditions>AttenThreshHigh).meanEstError]);
highAttenTrials.meanEstErrStd = std([attenDataSet(AttenConditions>AttenThreshHigh).meanEstError]);
highAttenTrials.meanErgodicityMean = mean([attenDataSet(AttenConditions>AttenThreshHigh).meanErgodicity]);
highAttenTrials.meanErgodicityVar = var([attenDataSet(AttenConditions>AttenThreshHigh).meanErgodicity]);
highAttenTrials.meanErgodicityStd = std([attenDataSet(AttenConditions>AttenThreshHigh).meanErgodicity]);

% Show Plot
colorMap = lines(4);
colorMap(1, :) = [0, 0, 0];
colorMap(3, :) = colorMap(2, :);
colorMap(4, :) = colorMap(2, :);
figure(1); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [12 8]);
% Reference
legHdl(1) = line([refTrial.meanErgodicityMean,refTrial.meanErgodicityMean], [-100, 100],...
    'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');
line([0, 1], [refTrial.meanEstErrMean, refTrial.meanEstErrMean],...
    'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');


% Plot raw data
plot(refTrial.meanErgodicityRaw, refTrial.meanEstErrRaw, 'o', 'color', colorMap(1,:),...
    'markerfacecolor', colorMap(1,:)*0.65, 'markersize', 4);
legHdl(2) = plot(lowAttenTrials.meanErgodicityRaw, lowAttenTrials.meanEstErrRaw, 'o', 'color', colorMap(2,:),...
    'markerfacecolor', colorMap(2,:)*0.65, 'markersize', 4);
plot(medAttenTrials.meanErgodicityRaw, medAttenTrials.meanEstErrRaw, 'o', 'color', colorMap(3,:),...
    'markerfacecolor', colorMap(3,:)*0.65, 'markersize', 4);
plot(highAttenTrials.meanErgodicityRaw, highAttenTrials.meanEstErrRaw, 'o', 'color', colorMap(4,:),...
    'markerfacecolor', colorMap(4,:)*0.65, 'markersize', 4);

xlabel('Mean Distance from Ergodicity', 'FontSize', 22);
ylabel('Mean Relative Tracking Error', 'FontSize', 22);
legend(legHdl, 'No Attenuation', 'With Attenuation', 'Location', 'southeast');

% Prettify figure
opt = [];
opt.BoxDim = [8,5];
opt.YLim = [40, 80];
opt.YTick = 40:10:80;
opt.XLim = [0.20, 0.45];
opt.XTick = 0.20:0.05:0.45;
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.ShowBox = 'off';
opt.FontName = 'Helvetica';
opt.IgnoreLines = 1;
setAxesProp(opt, gca);
set(gca,'YTickLabel', strcat(num2str((opt.YTick)'),'%'));
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.3, 0.25];
set(gca, 'Position', axesPosition);
% legend('off');
print(GEN_SAVE_PATH('fig5C.pdf'),'-dpdf');

%% Boxplot - Ergodic Metric vs. Attenuation
figure(2); clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [12 8]);
notBoxPlot(meanErgodicity(gAtten ~= 0), gAtten(gAtten ~= 0), ...
    'jitter', 0.2, 'jitterScale', 5, 'alpha', 0.5), hold on;
line([-5, max(gAtten)+1], [refTrial.meanErgodicityMean,refTrial.meanErgodicityMean], ...
    'LineStyle', '--', 'LineWidth', 2, 'Color', 'k');
xlabel('Wiggle Attenuation (dB)');
ylabel('Mean Distance from Ergodicity');

% Prettify figure
opt = [];
opt.BoxDim = [8,5];
opt.YLim = [0.20, 0.45];
opt.YTick = 0.20:0.05:0.45;
opt.XLim = [min(gAtten), max(gAtten)+5];
opt.XTick = [0, 25:25:150];
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.ShowBox = 'off';
opt.FontName = 'Helvetica';
opt.IgnoreLines = 1;
setAxesProp(opt, gca);
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.3, 0.25];
set(gca, 'Position', axesPosition);
print(GEN_SAVE_PATH('fig5B.pdf'),'-dpdf');

%% Boxplot - Mean Tracking Error vs. Attenuation
figure(3); clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [12 8]);
notBoxPlot(meanEstErrorMeanBelief(gAtten ~= 0), gAtten(gAtten ~= 0), ...
    'jitter', 0.2, 'jitterScale', 5, 'alpha', 0.5), hold on;
line([-5, max(gAtten)+1], [refTrial.meanEstErrMean,refTrial.meanEstErrMean], ...
    'LineStyle', '--', 'LineWidth', 2, 'Color', 'k');

xlabel('Wiggle Attenuation (dB)');
ylabel('Mean Relative Tracking Error');

% Prettify figure
opt = [];
opt.BoxDim = [8,5];
opt.YLim = [47.5, 80];
opt.XLim = [min(gAtten), max(gAtten)+5];
opt.YTick = 50:10:80;
opt.XTick = [0, 25:25:150];
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.ShowBox = 'off';
opt.FontName = 'Helvetica';
opt.IgnoreLines = 1;
setAxesProp(opt, gca);
set(gca,'YTickLabel', strcat(num2str((opt.YTick)'),'%'));
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.3, 0.25];
set(gca, 'Position', axesPosition);
print(GEN_SAVE_PATH('fig5A.pdf'),'-dpdf');
fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));
