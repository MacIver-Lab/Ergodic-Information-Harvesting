function makeFig6Plot(dataPath, savePath)
close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, '/wiggle_attenuation_sim', fname);
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile('./PublishedData/', 'animal_behavior_data', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
% Lock random number generator seed
rng(0);

%% Kinematic Energy of Electric Fish Behavior
load(GEN_BEHAVIOR_DATA_PATH('ElectricFish/Eigenmannia-sp-StrongSignal.mat'), 'strongSigData');
load(GEN_BEHAVIOR_DATA_PATH('ElectricFish/Eigenmannia-sp-WeakSignal.mat'), 'weakSigData');

% Fit fish behavioral trajectory into normalized simulation workspace
relativeEffortS = zeros(1, size(strongSigData,1));
relativeEffortW = zeros(1, size(weakSigData,1));

for i = 1:size(strongSigData,1)
    [relativeEffortS(i), perfDataS(i)] = calcForce(strongSigData{i,1}, strongSigData{i,2});
end
for i = 1:size(weakSigData,1)
    [relativeEffortW(i), perfDataW(i)] = calcForce(weakSigData{i,1}, weakSigData{i,2});
end

fishPowS_SEM = std([perfDataS.fishPow]) / sqrt(length([perfDataS.fishPow]));
fishPowW_SEM = std([perfDataW.fishPow]) / sqrt(length([perfDataW.fishPow]));

ref = [[perfDataS.refPow], [perfDataW.refPow]];
refPow_SEM = std(ref) / sqrt(length(ref));
fprintf('Weak signal     : %.3f +/- %.3f mW\n', mean([perfDataW.fishPow]), fishPowW_SEM);
fprintf('Strong signal   : %.3f +/- %.3f mW\n', mean([perfDataS.fishPow]), fishPowS_SEM);
fprintf('Reference signal: %.3f +/- %.3f mW\n', mean([ref]), refPow_SEM);

% Make Plot
figure(1); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
energyData = [relativeEffortS,relativeEffortW];
energyLabel = [ones(1,length(relativeEffortS)), 2*ones(1,length(relativeEffortW))];
notBoxPlot(energyData, energyLabel,...
    'jitter', 0.04, 'alpha', 0.5);
opt = [];
opt.BoxDim = [8, 5]*0.5;
opt.YLabel = 'Relative Energy'; % ylabel
opt.XLim = [0.5, 2.5];
opt.YLim = [0, 30];
opt.YTick = [0, 10:10:30];
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
ytick = get(gca, 'YTick');
set(gca,'YTickLabel', strcat(num2str((ytick)'),'x'));
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.4, 0.6];
set(gca, 'Position', axesPosition);

% Test statistics
p = kruskalwallis(energyData, energyLabel, 'off');
p = ceil(p * 1e3) / 1e3;
fprintf('Electrc Fish (Behavior)  : Kruskal-Wallis one-way ANOVA test, p < %.3f\n', p);

%% EIH relative tracking error vs. relative energy
load(GEN_DATA_PATH('PerfData.mat'));
targetTrajAmp = 0.2; % Amplitude of the target trajectory
gAtten = PerfData.gAtten;
% meanEstErrorMeanBelief 
% Unit: percentage tracking error with respect to target trajectory
%       amplitude
meanEstErrorMeanBelief = 100.0 * PerfData.meanEstErrorMeanBelief / targetTrajAmp;
% Categorize data fields according to wiggle attenuation gain
AttenConditions = unique(gAtten);
nAttenConditions = length(AttenConditions);
for idx = 1:nAttenConditions
    attenDataSet(idx).gAtten = AttenConditions(idx);
    attenDataSet(idx).meanEstError = meanEstErrorMeanBelief(gAtten == AttenConditions(idx));
    attenDataSet(idx).reEffort = PerfData.reEffort(gAtten == AttenConditions(idx));
end
% Threshold for low, medium, and high attenuation trials
%     Low    = [ 5 -  15] dB
%     Medium = (15 -  80] dB
%     High   = (80 - 150] dB
AttenThreshLow = 15;
AttenThreshHigh = 80;
% Exclude the reference
AttenConditions(1) = NaN;
refTrial.meanEstErrRaw = [attenDataSet(1).meanEstError];
refTrial.reEffortRaw = [attenDataSet(1).reEffort];
refTrial.meanEstErrMean = mean([attenDataSet(1).meanEstError]);
refTrial.meanEstErrStd = std([attenDataSet(1).meanEstError]);
refTrial.reEffortMean = mean([attenDataSet(1).reEffort]);
refTrial.reEffortStd = std([attenDataSet(1).reEffort]);

% Low Attenuation Trials
% Statistics
indices = AttenConditions <= AttenThreshLow;
lowAttenTrials.meanEstErrRaw = [attenDataSet(indices).meanEstError];
lowAttenTrials.reEffortRaw = [attenDataSet(indices).reEffort];
lowAttenTrials.meanEstErrMean = mean([attenDataSet(indices).meanEstError]);
lowAttenTrials.meanEstErrStd = std([attenDataSet(indices).meanEstError]);
lowAttenTrials.reEffortMean = mean([attenDataSet(indices).reEffort]);
lowAttenTrials.reEffortStd = std([attenDataSet(indices).reEffort]);

% Medium Attenuation Trials
% Statistics
indices = AttenConditions>AttenThreshLow & AttenConditions<=AttenThreshHigh;
medAttenTrials.meanEstErrRaw = [attenDataSet(indices).meanEstError];
medAttenTrials.reEffortRaw = [attenDataSet(indices).reEffort];
medAttenTrials.meanEstErrMean = mean([attenDataSet(indices).meanEstError]);
medAttenTrials.meanEstErrStd = std([attenDataSet(indices).meanEstError]);
medAttenTrials.reEffortMean = mean([attenDataSet(indices).reEffort]);
medAttenTrials.reEffortStd = std([attenDataSet(indices).reEffort]);

% High Attenuation Trials
% Statistics
indices = AttenConditions > AttenThreshHigh;
highAttenTrials.meanEstErrRaw = [attenDataSet(indices).meanEstError];
highAttenTrials.reEffortRaw = [attenDataSet(indices).reEffort];
highAttenTrials.meanEstErrMean = mean([attenDataSet(indices).meanEstError]);
highAttenTrials.meanEstErrStd = std([attenDataSet(indices).meanEstError]);
highAttenTrials.reEffortMean = mean([attenDataSet(indices).reEffort]);
highAttenTrials.reEffortStd = std([attenDataSet(indices).reEffort]);

% Show Plot
colorMap = lines(4);
colorMap(1, :) = [0, 0, 0];
colorMap(3, :) = colorMap(2, :);
colorMap(4, :) = colorMap(2, :);
axes; hold on;
% Plot raw data
plot(refTrial.reEffortRaw, refTrial.meanEstErrRaw, 'o', 'color', colorMap(1,:),...
    'markerfacecolor', colorMap(1,:)*0.65, 'markersize', 4);
legHdl(2) = plot(lowAttenTrials.reEffortRaw, lowAttenTrials.meanEstErrRaw, 'o', 'color', colorMap(2,:),...
    'markerfacecolor', colorMap(2,:)*0.65, 'markersize', 4);
plot(medAttenTrials.reEffortRaw, medAttenTrials.meanEstErrRaw, 'o', 'color', colorMap(3,:),...
    'markerfacecolor', colorMap(3,:)*0.65, 'markersize', 4);
plot(highAttenTrials.reEffortRaw, highAttenTrials.meanEstErrRaw, 'o', 'color', colorMap(4,:),...
    'markerfacecolor', colorMap(4,:)*0.65, 'markersize', 4);

xlabel('Relative Energy', 'FontSize', 22);
ylabel('Relative Tracking Error', 'FontSize', 22);

% Reference line
legHdl(1) = line('XData', [1, mean(refTrial.reEffortRaw)],...
    'YData', [mean(refTrial.meanEstErrRaw), mean(refTrial.meanEstErrRaw)],...
    'LineStyle', '--', 'LineWidth', 2, 'Color', colorMap(1,:));
line('XData', [mean(refTrial.reEffortRaw), mean(refTrial.reEffortRaw)],...
    'YData', [mean(refTrial.meanEstErrRaw), 80],...
    'LineStyle', '--', 'LineWidth', 2, 'Color', colorMap(1,:));

% legend
legend(legHdl, 'No Attenuation', 'With Attenuation', 'Location', 'best');

% Prettify figure
opt = [];
opt.BoxDim = [8,5]*0.5;
opt.YLim = [45, 80];
opt.YTick = 50:10:80;
opt.XLim = [1, 35];
opt.XTick = [1, 10:10:35];
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.ShowBox = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 14;
opt.IgnoreLines = 1;
setAxesProp(opt, gca);
set(gca,'XTickLabel', strcat(num2str((opt.XTick)'),'x'));
set(gca,'YTickLabel', strcat(num2str((opt.YTick)'),'%'));
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.4, 0.2];
set(gca, 'Position', axesPosition);
print(gcf,'-dpdf',GEN_SAVE_PATH('fig6.pdf'));


function [fTraj, rTraj] = normalizeFishTraj(fTraj, rTraj)
%% Fit fish behavioral trajectory into normalized simulation workspace
refSimAmp = 0.2;  % Amplitude of simulated target trajectory
refSimMean = 0.5; % Center of simulated target trajectory
peakH = findpeaks(rTraj, ...
    'SortStr', 'descend', ...
    'MinPeakHeight', max(rTraj)-5);
peakL = abs(findpeaks(-rTraj, ...
    'SortStr', 'descend', ...
    'MinPeakHeight', max(-rTraj)-5));
refAmpPixel = (mean(peakH) - mean(peakL)) / 2.0;
% Remove offset
rTraj = rTraj - (refAmpPixel + mean(peakL));
fTraj = fTraj - mean(fTraj);
% Scale and offset to match simulated amplitude
rTraj = rTraj / refAmpPixel * refSimAmp + refSimMean;
fTraj = fTraj / refAmpPixel * refSimAmp + refSimMean;

clf; hold on;
plot(rTraj);
plot(fTraj);

function [relativeEffort, perfData] = calcForce(fTraj, rTraj)
%% Estimate the amount of force required for the acceleration
%  Using simple F = ma model
USE_DRAG = 0;
IGNORE_NEG_EFFORT = 0;

% Video has a 46.8 pixel/cm (DPC) pixel density in the longitudinal axis
% The refuge measures 15 cm long
videoDPC = 46.8;
refugeLenCM = 15.0;
fishLenCM = 14.36;
% Video FPS
videoFPS = 60.0; % [frame/second]

% Fish body weight
% This is estimated by using black ghost knifefish as a reference
% According to MacI10a Energy-Information Trade-Offs between Movement and Sensing
% black ghost knifefish weights 23 grams with a body length of 19 cm
% Therefore, we can estimate eigenmannia's weight considering its body
% length is 14.36 cm by assuming it has the similar density, note that this
% will likely to slightly overestimate the weight since eigenmannia's body
% is slimmer than black ghost, especially near the caudal part. Hence, we 
% discounted the overall weight by a ratio of 0.8 to take that into
% consideration.
% Added mass is taken from Clai08a
addedMass = 6.04;
fishWeightGram = 23 + addedMass;
% fishWeightGram = 14.36 * (23.0 / 19.0) + addedMass;
fishWeightGram = fishWeightGram * 1.0;

% Drag model
rho = 1.0;        % Density of fluit, water = 1.0 [g/cm^3]
area = 30;        % unit [cm^2]
dragCoef = 0.04;  % Drag coefficient, streamlined body = 0.04, dimensionless
dragConstant = 0.71; % Drag constant, approximated by reverse computing with the 
                     % fact that black ghost knifefish has an approximated 2.0 mN
                     % drag when swimming at 15cm/sec velocity which lead to 
                     % a drag constant of around 0.89, take a 0.8 discount 
                     % and we got 0.80
dragForce = @(v) dragConstant * (v.^2); % unit [g*cm/sec^2]

%% Filter raw trajectory
% Position unit [pixel]
stopFreq = 1.6; % [Hz]
fTraj = LPF(fTraj, videoFPS, stopFreq);
rTraj = LPF(rTraj, videoFPS, stopFreq);
% flatten
fTraj = fTraj(:);
rTraj = rTraj(:);

%% Approximate velocity through time delta
% Unit [cm/sec]
fVel = diff(fTraj) / videoDPC * videoFPS;
rVel = diff(rTraj) / videoDPC * videoFPS;
% fprintf('Maximum fish velocity       = %.4f cm/sec\n', max(abs(fVel)));

%% Approximate acceleration through time delta
% Unit [cm/sec^2]
fAcc = [0; diff(fVel)] * videoFPS;
rAcc = [0; diff(rVel)] * videoFPS;
% fprintf('Maximum fish acceleration   = %.4f cm/sec^2\n', max(abs(fAcc)));

%% Estimate mechanical force for locomotion
% Force unit [g*cm/sec^2] = [1e-5 Newton]
if USE_DRAG
    fForce = zeros(size(fAcc), 'double');
    bwdSwimIdx = find(fVel < 0);
    fwdSwimIdx = find(fVel >= 0);
    fForce(fwdSwimIdx) = fishWeightGram * fAcc(fwdSwimIdx) + dragForce(fVel(fwdSwimIdx));
    fForce(bwdSwimIdx) = fishWeightGram * fAcc(bwdSwimIdx) - dragForce(fVel(bwdSwimIdx));
    % Reference trajectory
    rForce = zeros(size(rAcc), 'double');
    bwdSwimIdx = find(rVel < 0);
    fwdSwimIdx = find(rVel >= 0);
    rForce(fwdSwimIdx) = fishWeightGram * rAcc(fwdSwimIdx) + dragForce(rVel(fwdSwimIdx));
    rForce(bwdSwimIdx) = fishWeightGram * rAcc(bwdSwimIdx) - dragForce(rVel(bwdSwimIdx));
else
    fForce = fishWeightGram * fAcc;
    rForce = fishWeightGram * rAcc;
end

%% Estimate Mechanical Effort
% Compute instant power, unit [g*cm^2/sec^3] = [1e-7 W]
if IGNORE_NEG_EFFORT
    fPow = fForce .* fVel;
    rPow = rForce .* rVel;
    fPow(fPow < 0) = 0;
    rPow(rPow < 0) = 0;
else
    fPow = abs(fForce .* fVel);
    rPow = abs(rForce .* rVel);
end

% Integrate to get final mechanical effort estimation
% unit [g*cm^2/sec^2] = [1e-7 J]
dt = 1 / videoFPS;
time = dt:dt:(length(fForce)/videoFPS);
fEffort = trapz(time, fPow);
rEffort = trapz(time, rPow);

% Convert unit to micro-J
% unit [g*cm^2/sec^2] = [1e-7 J]
fEffort = 1e-1 * fEffort;
rEffort = 1e-1 * rEffort;

relativeEffort = fEffort / rEffort;

perfData.fEffort = fEffort;
perfData.rEffort = rEffort;
perfData.fishPow = 1e-4 * max(fPow);
perfData.refPow = 1e-4 * max(rPow);
perfData.fishF = max(abs(fForce));
perfData.fishVel = max(abs(fVel));
perfData.fishAccel = max(abs(fAcc));
perfData.refVel = max(abs(rVel));
perfData.refAccel = max(abs(rAcc));