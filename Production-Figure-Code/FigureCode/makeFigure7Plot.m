function makeFigure7Plot(dataPath, savePath)

close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, 'behavior_sim', fname);
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile('./PublishedData/', 'animal_behavior_data', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
% Lock random number generator seed
rng(0);

%% Electric Fish Behavioral Trials
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
    'PaperSize', [30, 30]);
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
axesPosition(1:2) = [0.5, 0.5];
set(gca, 'Position', axesPosition);
title('Electric Fish (Behavior)', 'FontSize', 12)
% Test statistics (non-parametric ANOVA)
p = kruskalwallis(magList, magLabel, 'off');
fprintf('Statistics (electric fish behavior): Kruskal-wallis test - p = %.4f (n = %d)\n', ...
    p, length(magList));

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
axesPosition(1:2) = [0.1, 0.5];
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
axesPosition(1:2) = [0.7, 0.5];
set(gca, 'Position', axesPosition);
title('Electric Fish (Simulation)', 'FontSize', 12)
% Test statistics (non-parametric ANOVA)
p = kruskalwallis(magList, magLabel, 'off');
fprintf('Statistics (electric fish simulation): Kruskal-wallis test - p = %.4f (n = %d)\n', ...
    p, length(magList));

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
axesPosition(1:2) = [0.3, 0.5];
set(gca, 'Position', axesPosition);
title('Electric Fish (Simulation)', 'FontSize', 12)
patch(...
    'Vertices', [0.25, 0; 1, 0; 1, 250; 0.25, 250], ...
    'Faces', [1, 2, 3, 4], ...
    'FaceColor', 'k', 'FaceAlpha', 0.1, ...
    'EdgeAlpha', 0)
print(gcf,'-dpdf',GEN_SAVE_PATH('fig7A-D.pdf'));

%% Mole behavioral data
cmap = [56, 180, 74; ...
    161, 30, 34] / 255;
barColors = [[75, 111, 182]/255; cmap];
% Behavior
fNames = dir(GEN_BEHAVIOR_DATA_PATH('Mole/fig*.mat'));
idx_lSNR = 1;
idx_hSNR = 1;
findAngle = @(a, b) acosd(min(1,max(-1, a(:).' * b(:) / norm(a) / norm(b) )));
for i = 1:length(fNames)
    moleData = load([fNames(i).folder, '/', fNames(i).name]);
    if ~contains(fNames(i).name, 'b1')
        % all figures other than b1* have inverted input
        % due to digitization
        moleData.molePath = moleData.molePath(end:-1:1, :);
        moleData.refPath = moleData.refPath(end:-1:1, :);
    end
    refVec = moleData.refPath(2, :) - moleData.refPath(1, :);
    refAngle = findAngle(refVec/norm(refVec), [0, 1]);
    if refVec(1) > 0
        refAngle = -refAngle;
    end
    rotPath = rotate2D(moleData.molePath, refAngle);
    rotRef = rotate2D(moleData.refPath, refAngle);
    % Process data
    if ~contains(fNames(i).name,'normal')
        %%%%%%- One-side nostril block - Low SNR -%%%%%%
        % FFT
        [fftTicks, MoleFreqMagRawX_lSNR(:, idx_lSNR)] = ...
            decomposeFourierMag(rotPath(:, 1), [], 1/15, 2);
        [~, MoleFreqMagRawY_lSNR(:, idx_lSNR)] = ...
            decomposeFourierMag(rotPath(:, 2), [], 1/15, 2);
        % Compute averaged spectrum magnitude
        MoleFreqMagX_lSNR(idx_lSNR) = mean(MoleFreqMagRawX_lSNR(:, idx_lSNR));
        MoleFreqMagY_lSNR(idx_lSNR) = mean(MoleFreqMagRawY_lSNR(:, idx_lSNR));
        idx_lSNR = idx_lSNR + 1;
    else
        %%%%%%- Normal - High SNR -%%%%%%
        % FFT
        [fftTicks, MoleFreqMagRawX_hSNR(:, idx_hSNR)] = ...
            decomposeFourierMag(rotPath(:, 1), [], 1/15, 2);
        [~, MoleFreqMagRawY_hSNR(:, idx_hSNR)] = ...
            decomposeFourierMag(rotPath(:, 2), [], 1/15, 2);
        % Compute averaged spectrum magnitude
        MoleFreqMagX_hSNR(idx_hSNR) = mean(MoleFreqMagRawX_hSNR(:, idx_hSNR));
        MoleFreqMagY_hSNR(idx_hSNR) = mean(MoleFreqMagRawY_hSNR(:, idx_hSNR));
        idx_hSNR = idx_hSNR + 1;
    end
end
% Lateral
magListX = [...
    MoleFreqMagX_hSNR, ...
    MoleFreqMagX_lSNR];
magListX = (magListX - mean(MoleFreqMagX_hSNR)) / max(magListX);
magLabelX = [...
    1*ones(1, length(MoleFreqMagX_hSNR)), ...
    2*ones(1, length(MoleFreqMagX_lSNR))];
%%%% Main plot %%%%
figure(2); clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [30, 30]);
%%% Spectrum (Lateral) %%%
axes; hold on;
nrmFactor = max([MoleFreqMagRawX_lSNR(:, 4); MoleFreqMagRawX_hSNR(:, 4)]);
plot(fftTicks, MoleFreqMagRawX_lSNR(:, 4)/nrmFactor, '-', ...
    'LineWidth', 2, 'Color', [0.5, 0, 0, 0.9]);
plot(fftTicks, MoleFreqMagRawX_hSNR(:, 4)/nrmFactor, '-', ...
    'LineWidth', 2, 'Color', [0, 0.5, 0, 0.9]);
set(gca, 'YTick', [0, 1]);
set(gca, 'XTick', [0, 2]);
set(gca, 'FontName', 'Helvetica');
xlabel('Frequency');
ylabel('Fourier Magnitude');
hLine = findobj(gca,'Type','line');
patch(...
    'Vertices', [0, 0; 2, 0; 2, 1.7e5; 0, 1.7e5], ...
    'Faces', [1, 2, 3, 4], ...
    'FaceColor', 'k', 'FaceAlpha', 0.1, ...
    'EdgeAlpha', 0)
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XLim = [0, 2];
opt.YLim = [0, 1];
opt.YTick = [0, 1];
opt.FontName = 'Helvetica';
opt.FontSize = 12;
setAxesProp(opt);
hLine(1).LineWidth = 2;
hLine(1).Color = [0, 0.5, 0, 0.9];
hLine(2).LineWidth = 2;
hLine(2).Color = [0.5, 0, 0, 0.9];
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.1, 0.5];
set(gca, 'Position', axesPosition);
title('Mole (Behavior)', 'FontSize', 12)
%%% Spectrum Statistics %%%
hAxe = axes; 
notBoxPlot(magListX, magLabelX);
xlim([0.5, 2.5]);
ylim([-0.3, 1]);
set(gca, 'YTick', [0, 1]);
set(gca, 'XTick', [1, 2]);
set(gca, 'XTickLabel', ...
    {'\color[rgb]{0.1961,0.7059,0.2902}Strong Signal', ...
    '\color[rgb]{0.6353,0,0}Weak Signal'});
set(gca, 'FontName', 'Helvetica');
ylabel('Fourier Magnitude');
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 12;
setAxesProp(opt);
legend('off');
hLine = findobj(gca,'Type','line');
hLine(1).LineWidth = 4;
hLine(1).Color = cmap(2, :);
hLine(2).LineWidth = 4;
hLine(2).Color = cmap(1, :);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.5, 0.5];
set(gca, 'Position', axesPosition);
title('Mole (Behavior)', 'FontSize', 12)
% Test statistics (non-parametric ANOVA)
pX = kruskalwallis(magListX, magLabelX, 'off');
fprintf('Statistics (mole behavior): Kruskal-wallis test - p = %.4f (n = %d)\n', ...
    pX, length(magListX));


% EIH
EH_lSNR_files = dir(GEN_DATA_PATH('EIH-Mole-WeakSignal*.mat'));
EH_hSNR_files = dir(GEN_DATA_PATH('EIH-Mole-StrongSignal*.mat'));
for i = 1:length(EH_lSNR_files)
    Mole_EIH_lSNR = load(GEN_DATA_PATH(EH_lSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    Mole_EIH_hSNR = load(GEN_DATA_PATH(EH_hSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    % crop the segment before sensor crosses over the target
    % the first time
    crossIdx_lSNR = find(Mole_EIH_lSNR.sTrajList < Mole_EIH_lSNR.oTrajList(1), 1);
    crossIdx_hSNR = find(Mole_EIH_hSNR.sTrajList > Mole_EIH_hSNR.oTrajList(1), 1);
    % FFT
    [fftTicks, MoleEIHMag_lSNR(:, i)] = ...
        decomposeFourierMag(Mole_EIH_lSNR.sTrajList(crossIdx_lSNR:end), [], Mole_EIH_lSNR.dt);
    [fftTicks, MoleEIHMag_hSNR(:, i)] = ...
        decomposeFourierMag(Mole_EIH_hSNR.sTrajList(crossIdx_hSNR:end), [], Mole_EIH_hSNR.dt);
    % Compute averaged spectrum magnitude
    MoleFreqMagEIH_lSNR(i) = mean(MoleEIHMag_lSNR(:, i));
    MoleFreqMagEIH_hSNR(i) = mean(MoleEIHMag_hSNR(:, i));
end
magListEIH = [...
    MoleFreqMagEIH_hSNR, ...
    MoleFreqMagEIH_lSNR];
magListEIH = (magListEIH - mean(MoleFreqMagEIH_hSNR)) / max(magListEIH);
magLabelEIH = [...
    1*ones(1, length(MoleFreqMagEIH_hSNR)), ...
    2*ones(1, length(MoleFreqMagEIH_lSNR))];
%%% Spectrum Statistics %%%
axes; hold on;
notBoxPlot(magListEIH, magLabelEIH);
xlim([0.5, 2.5])
ylim([-0.3, 1])
set(gca, 'YTick', [0, 1]);
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
axesPosition(1:2) = [0.7, 0.5];
set(gca, 'Position', axesPosition);
title('Mole (Simulation)', 'FontSize', 12)
%%% Spectrum %%%
axes; hold on;
nrmFactor = max([MoleEIHMag_lSNR(:, 2);MoleEIHMag_hSNR(:, 2)]);
plot(fftTicks, MoleEIHMag_lSNR(:, 2)/nrmFactor, '-', ...
    'LineWidth', 2, 'Color', [0.5, 0, 0, 0.9]);
plot(fftTicks, MoleEIHMag_hSNR(:, 2)/nrmFactor, '-', ...
    'LineWidth', 2, 'Color', [0, 0.5, 0, 0.9]);
set(gca, 'YTick', [0, 1]);
set(gca, 'XTick', [0, 0.5, 1]);
set(gca, 'FontName', 'Helvetica');
xlabel('Frequency');
ylabel('Fourier Magnitude');
hLine = findobj(gca,'Type','line');
patch(...
    'Vertices', [0, 0; 1, 0; 1, 1.7e5; 0, 1.7e5], ...
    'Faces', [1, 2, 3, 4], ...
    'FaceColor', 'k', 'FaceAlpha', 0.1, ...
    'EdgeAlpha', 0)
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.ShowBox = 'off';
opt.YLim = [0, 1];
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 12;
setAxesProp(opt);
hLine(1).LineWidth = 2;
hLine(1).Color = [0, 0.5, 0, 0.9];
hLine(2).LineWidth = 2;
hLine(2).Color = [0.5, 0, 0, 0.9];
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.3, 0.5];
set(gca, 'Position', axesPosition);
title('Mole (Simulation)', 'FontSize', 12)
p = kruskalwallis(magListEIH, magLabelEIH, 'off');
fprintf('Statistics (mole simulation): Kruskal-wallis test - p = %.4f (n = %d)\n', ...
    p, length(magListEIH));

print(gcf,'-dpdf',GEN_SAVE_PATH('fig7E-H.pdf'));

%% Cockroach
cmap = [56, 180, 74; ...
    161, 30, 34] / 255;
barColors = [[75, 111, 182]/255; cmap];
cockroachData = load(GEN_BEHAVIOR_DATA_PATH('Cockroach/cockroach_data.mat'));
getCockroachFilename = @(x) dir(GEN_BEHAVIOR_DATA_PATH(...
    ['Cockroach/raw_tracked_data/', x, '*_side_c.csv']));
for i = 1:length(cockroachData.trial_c12.trial)
    file = getCockroachFilename(cockroachData.trial_c12.trial{i});
    fprintf('Parsing file %s...\n', file.name);
    rawData = csvread([file.folder, '/', file.name]);
    rawData = parseCockroachData(rawData);
    % FFT
    [fftTicks, fftMagX_lSNR] = decomposeFourierMag(rawData.head(:, 1), [], 1/15, 2);
    [~, fftMagY_lSNR] = decomposeFourierMag(rawData.head(:, 2), [], 1/15, 2);
    % Keep a copy of FFT for showing illustrational frequency spectrum
    pltMagX_lSNR(:, i) = fftMagX_lSNR;
    pltMagY_lSNR(:, i) = fftMagY_lSNR;
    % Compute averaged spectrum magnitude
    CockroachFreqMagX_lSNR(i) = mean(fftMagX_lSNR);
    CockroachFreqMagY_lSNR(i) = mean(fftMagY_lSNR);
end
for i = 1:length(cockroachData.trial_c4.trial)
    file = getCockroachFilename(cockroachData.trial_c4.trial{i});
    fprintf('Parsing file %s...\n', file.name);
    rawData = csvread([file.folder, '/', file.name]);
    rawData = parseCockroachData(rawData);
    % FFT
    [fftTicks, fftMagX_hSNR] = decomposeFourierMag(rawData.head(:, 1), [], 1/15, 2);
    [~, fftMagY_hSNR] = decomposeFourierMag(rawData.head(:, 2), [], 1/15, 2);
    % Keep a copy of FFT for showing illustrational frequency spectrum
    pltMagX_hSNR(:, i) = fftMagX_hSNR;
    pltMagY_hSNR(:, i) = fftMagY_hSNR;
    % Compute averaged spectrum magnitude
    CockroachFreqMagX_hSNR(i) = mean(fftMagX_hSNR);
    CockroachFreqMagY_hSNR(i) = mean(fftMagY_hSNR);
end
magListX = [...
    CockroachFreqMagX_hSNR, ...
    CockroachFreqMagX_lSNR];
magLabelX = [...
    0.85*ones(1, length(CockroachFreqMagX_hSNR)), ...
    1.15*ones(1, length(CockroachFreqMagX_lSNR))];
magListY = [...
    CockroachFreqMagY_hSNR, ...
    CockroachFreqMagY_lSNR];
magLabelY = [...
    1*ones(1, length(CockroachFreqMagY_hSNR)), ...
    2*ones(1, length(CockroachFreqMagY_lSNR))];

%%%% Main plot %%%%
figure(3); clf; 
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [30, 30]);
%%% Illustrational frequency spectrum %%%
axes; hold on;
nrmFactor = max([pltMagY_lSNR(:, 7);pltMagY_hSNR(:, 7)]);
plot(fftTicks, pltMagY_lSNR(:, 7)/nrmFactor, '-', ...
    'LineWidth', 2, 'Color', [0.5, 0, 0, 0.9]);
plot(fftTicks, pltMagY_hSNR(:, 7)/nrmFactor, '-', ...
    'LineWidth', 2, 'Color', [0, 0.5, 0, 0.9]);
set(gca, 'YTick', [0, 1]);
set(gca, 'XTick', [0, 2]);
set(gca, 'FontName', 'Helvetica');
ylim([0, 1]);
xlabel('Frequency');
ylabel('Fourier Magnitude');
hLine = findobj(gca,'Type','line');
yRange = ylim;
patch(...
    'Vertices', [0, 0; 2, 0; 2, yRange(2); 0, yRange(2)], ...
    'Faces', [1, 2, 3, 4], ...
    'FaceColor', 'k', 'FaceAlpha', 0.1, ...
    'EdgeAlpha', 0);
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 14;
setAxesProp(opt);
hLine(1).LineWidth = 2;
hLine(1).Color = [0, 0.5, 0, 0.9];
hLine(2).LineWidth = 2;
hLine(2).Color = [0.5, 0, 0, 0.9];
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.1, 0.5];
set(gca, 'Position', axesPosition);
title('Cockroach Behavior', 'FontSize', 12)
%%% Statistics %%%
axes;
magListY = (magListY - mean(CockroachFreqMagY_hSNR)) / max(magListY);
notBoxPlot(magListY, magLabelY);
xlim([0.5, 2.5])
ylim([-0.3, 1])
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
opt.YTick = [0, 1];
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 14;
setAxesProp(opt);
legend('off');
hLine(1).LineWidth = 4;
hLine(1).Color = cmap(2, :);
hLine(2).LineWidth = 4;
hLine(2).Color = cmap(1, :);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.5, 0.5];
set(gca, 'Position', axesPosition);
title('Cockroach Behavior', 'FontSize', 12)
hAxe = gca;
% Test statistics (non-parametric ANOVA)
pX = kruskalwallis(magListX, magLabelX, 'off');
pY = kruskalwallis(magListY, magLabelY, 'off');
fprintf('Statistics (cockroach behavior): Kruskal-wallis test - p = %.4f (n = %d)\n', ...
    pY, length(magListY));


% EIH
EH_lSNR_files = dir(GEN_DATA_PATH('EIH-Cockroach-WeakSignal*.mat'));
EH_hSNR_files = dir(GEN_DATA_PATH('EIH-Cockroach-StrongSignal*.mat'));
for i = 1:length(EH_lSNR_files)
    Cockroach_EIH_lSNR = load(GEN_DATA_PATH(EH_lSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');
    Cockroach_EIH_hSNR = load(GEN_DATA_PATH(EH_hSNR_files(i).name), 'oTrajList', 'sTrajList', 'dt');   
    % FFT
    [Cockroach_EIH_hSNR.fftTicks, Cockroach_EIH_hSNR.fftMag] = ...
        decomposeFourierMag(Cockroach_EIH_hSNR.sTrajList, [], Cockroach_EIH_hSNR.dt);
    [Cockroach_EIH_lSNR.fftTicks, Cockroach_EIH_lSNR.fftMag] = ...
        decomposeFourierMag(Cockroach_EIH_lSNR.sTrajList, [], Cockroach_EIH_lSNR.dt);
    % Compute averaged spectrum magnitude
    CockroachFreqMagEIH_hSNR(i) = mean(Cockroach_EIH_hSNR.fftMag);
    CockroachFreqMagEIH_lSNR(i) = mean(Cockroach_EIH_lSNR.fftMag);
end
magListEIH = [...
    CockroachFreqMagEIH_hSNR, ...
    CockroachFreqMagEIH_lSNR];
magLabelEIH = [...
    1*ones(1, length(CockroachFreqMagEIH_hSNR)), ...
    2*ones(1, length(CockroachFreqMagEIH_lSNR))];
figure(3); 
%%% Spectrum %%%
axes; hold on;
nrmFactor = max([Cockroach_EIH_lSNR.fftMag, Cockroach_EIH_hSNR.fftMag]);
plot(Cockroach_EIH_lSNR.fftTicks, Cockroach_EIH_lSNR.fftMag/nrmFactor, '-', ...
    'LineWidth', 2, 'Color', [0.5, 0, 0, 0.9]);
plot(Cockroach_EIH_hSNR.fftTicks, Cockroach_EIH_hSNR.fftMag/nrmFactor, '-', ...
    'LineWidth', 2, 'Color', [0, 0.5, 0, 0.9]);
set(gca, 'YTick', [0, 1]);
set(gca, 'XTick', [0, 1]);
set(gca, 'FontName', 'Helvetica');
xlabel('Frequency');
ylabel('Fourier Magnitude');
hLine = findobj(gca,'Type','line');
ylim([0, 1]);
yRange = ylim;
patch(...
    'Vertices', [0, 0; 1, 0; 1, yRange(2); 0, yRange(2)], ...
    'Faces', [1, 2, 3, 4], ...
    'FaceColor', 'k', 'FaceAlpha', 0.1, ...
    'EdgeAlpha', 0);
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.YLim = yRange;
opt.FontName = 'Helvetica';
opt.FontSize = 14;
setAxesProp(opt);
hLine(1).LineWidth = 2;
hLine(1).Color = [0, 0.5, 0, 0.9];
hLine(2).LineWidth = 2;
hLine(2).Color = [0.5, 0, 0, 0.9];
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.3, 0.5];
set(gca, 'Position', axesPosition);
title('Cockroach Simulation', 'FontSize', 12)
%%% Statistics %%%
axes;
magListEIH = (magListEIH - mean(CockroachFreqMagEIH_hSNR)) / max(magListEIH);
notBoxPlot(magListEIH, magLabelEIH);
xlim([0.5, 2.5])
ylim([-0.3, 1])
set(gca, 'YTick', [0, 1]);
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
opt.FontName = 'Helvetica';
opt.FontSize = 14;
setAxesProp(opt);
legend('off');
hLine(1).LineWidth = 4;
hLine(1).Color = cmap(2, :);
hLine(2).LineWidth = 4;
hLine(2).Color = cmap(1, :);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.7, 0.5];
set(gca, 'Position', axesPosition);
title('Cockroach Simulation (Lateral Position)', 'FontSize', 12)
p = kruskalwallis(magListEIH, magLabelEIH, 'off');
fprintf('Statistics (cockroach simulation): Kruskal-wallis test - p = %.4f (n = %d)\n', ...
    p, length(magListEIH));

print(gcf,'-dpdf',GEN_SAVE_PATH('fig7I-L.pdf'));


%% Hawkmoth
FreqResolution = 0.001;
CutOffFreq = 3.5;
LPFHighCutFreq = 6.0;
NPeaks = 18;
hsnrFiles = dir(GEN_DATA_PATH('EIH-Moth-StrongSignal-*.mat'));
lsnrFiles = dir(GEN_DATA_PATH('EIH-Moth-WeakSignal-*.mat'));
hsnrBodeGain = zeros(length(hsnrFiles), NPeaks, 'double');
lsnrBodeGain = zeros(length(lsnrFiles), NPeaks, 'double');
locs = [53, 73, 126, 177, 276, 325, 427, 475, 576, 725, 926, ...
    1075, 1325, 1526, 1975, 2225, 2825, 3426];
calcDist = @(x) sum(abs(diff(x)));
for i = 1:length(hsnrFiles)
    load([hsnrFiles(i).folder, '/', hsnrFiles(i).name], ...
        'dt', 'oTrajList', 'sTrajList');      
    reHSNR_EIH(i) = calcDist(sTrajList) / calcDist(oTrajList);
    Fs = 1 / dt;
    nFFTSamples = round(Fs / FreqResolution);
    CutOffIdx = round(CutOffFreq / FreqResolution);
    freqIdx = FreqResolution:FreqResolution:CutOffFreq;
    oPos = oTrajList;
    sPos = sTrajList;
    oPos = LPF(oPos, Fs, LPFHighCutFreq);
    sPos = LPF(sPos, Fs, LPFHighCutFreq);
    
    oFFT = fft(oPos-mean(oPos), nFFTSamples);
    sFFT = fft(sPos-mean(sPos), nFFTSamples);
    oFFT = abs(oFFT(1:CutOffIdx));
    sFFT = abs(sFFT(1:CutOffIdx));
    oLocs = findLocalPeaks(oFFT, locs);
    sFFTLocs = findLocalPeaks(sFFT, oLocs);
    hsnrBodeGain(i, :) = sFFT(sFFTLocs) ./ oFFT(oLocs);
    fprintf('file %s processed\n', hsnrFiles(i).name);
end
for i = 1:length(lsnrFiles)
    load([lsnrFiles(i).folder, '/', lsnrFiles(i).name], ...
        'dt', 'oTrajList', 'sTrajList');        
    reLSNR_EIH(i) = calcDist(sTrajList) / calcDist(oTrajList);
    Fs = 1 / dt;
    nFFTSamples = round(Fs / FreqResolution);
    CutOffIdx = round(CutOffFreq / FreqResolution);
    freqIdx = FreqResolution:FreqResolution:CutOffFreq;
    oPos = oTrajList;
    sPos = sTrajList;
    oPos = LPF(oPos, Fs, LPFHighCutFreq);
    sPos = LPF(sPos, Fs, LPFHighCutFreq);
    
    oFFT = fft(oPos-mean(oPos), nFFTSamples);
    sFFT = fft(sPos-mean(sPos), nFFTSamples);
    oFFT = abs(oFFT(1:CutOffIdx));
    sFFT = abs(sFFT(1:CutOffIdx));
    oLocs = findLocalPeaks(oFFT, locs);
    sFFTLocs = findLocalPeaks(sFFT, oLocs);
    lsnrBodeGain(i, :) = sFFT(sFFTLocs) ./ oFFT(oLocs);
    fprintf('file %s processed\n', lsnrFiles(i).name);
end
meanHSNR = mean(hsnrBodeGain);
meanLSNR = mean(lsnrBodeGain);
semHSNR = 1.96 * std(hsnrBodeGain) ./ sqrt(size(hsnrBodeGain,1));
semLSNR = 1.96 * std(lsnrBodeGain) ./ sqrt(size(lsnrBodeGain,1));
medFreqGainHSNR = mean(hsnrBodeGain(:, 8:11), 2, 'omitnan');
medFreqGainLSNR = mean(lsnrBodeGain(:, 8:11), 2, 'omitnan');
lowFreqGainHSNR = mean(hsnrBodeGain(:, 1:7), 2, 'omitnan');
lowFreqGainLSNR = mean(lsnrBodeGain(:, 1:7), 2, 'omitnan');
medFreqGainData = [medFreqGainHSNR; medFreqGainLSNR];
medFreqGainLabel = [ones(length(medFreqGainHSNR),1); 2*ones(length(medFreqGainLSNR),1)];
cmap = [56, 180, 74; ...
    161, 30, 34] / 255;

figure(4); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [30, 30]);
freq = [0.2;0.3;0.5;0.7;1.100;1.300;1.700;1.900;2.300;2.900;3.700;4.300;5.300;6.100;7.900;8.900;11.30;13.7];
boundedline(freq, [meanHSNR; meanLSNR], [semLSNR',semHSNR'], ...
    'alpha', 'transparency', 0.4, 'cmap', cmap);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
hLines = findobj(gca, 'Type', 'line');
yRange = ylim;
patch(...
    'Vertices', [freq(8)-0.1, 1e-5; freq(11)+0.25, 1e-5; freq(11)+0.25, 3.5; freq(8)-0.1, 3.5], ...
    'Faces', [1, 2, 3, 4], ...
    'FaceColor', 'k', 'FaceAlpha', 0.1, ...
    'EdgeAlpha', 0);
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.XLabel = 'Frequency'; % xlabel
opt.YLabel = 'Gain'; % ylabel
opt.XLim = [min(freq'), max(freq')];
opt.YLim = [0.10, 3.3];
opt.YTick = [0.2, 1];
opt.ShowBox = 'off';
opt.XMinorTick = 'on';
opt.YMinorTick = 'on';
opt.FontName = 'Helvetica';
opt.FontSize = 16;
% opt.Colors = cmap;
setAxesProp(opt);
hLines(1).Marker = '.';
hLines(2).Marker = '.';
hLines(1).MarkerSize = 16;
hLines(2).MarkerSize = 16;
hLines(1).LineWidth = 2;
hLines(2).LineWidth = 2;
hLines(1).Color = cmap(2, :);
hLines(2).Color = cmap(1, :);
legend('off');
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', {'0.2', '1'})
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.3, 0.5];
set(gca, 'Position', axesPosition);
title('Simulation', 'FontSize', 16)

%%BEHAVIOR DATA%%
dat = load(GEN_BEHAVIOR_DATA_PATH('Moth/MothBodeData_Macroglossum.mat'));
dat = dat.AnalysisMacroglossum;
for i = 1:3
    g{i} = dat(dat.Group == i, 'Gain').Gain;
    g{i} = reshape(g{i}, length(freq), length(g{i})/length(freq));
    meanG{i} = mean(g{i}, 2, 'omitnan');
    semG{i} = 1.96 * std(g{i}, 0, 2, 'omitnan') ./ sqrt(size(g{i},2));
    lowGain{i} = mean(g{i}(1:7, :), 'omitnan')';
    medGain{i} = mean(g{i}(8:11, :), 'omitnan')';
end
axes; hold on;
boundedline(freq, meanG{1}, semG{1}, ...
    'alpha', 'transparency', 0.4, 'cmap', cmap(1,:));
boundedline(freq, meanG{3}, semG{3}, ...
    'alpha', 'transparency', 0.4, 'cmap', cmap(2,:));
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
hLines = findobj(gca, 'Type', 'line');
yRange = ylim;
patch(...
    'Vertices', [freq(8)-0.1, 1e-5; freq(11)+0.25, 1e-5; freq(11)+0.25, 3.5; freq(8)-0.1, 3.5], ...
    'Faces', [1, 2, 3, 4], ...
    'FaceColor', 'k', 'FaceAlpha', 0.1, ...
    'EdgeAlpha', 0);
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.XLabel = 'Frequency'; % xlabel
opt.YLabel = 'Gain'; % ylabel
opt.XLim = [min(freq'), max(freq')];
opt.YLim = [0.10, 3.3];
opt.YTick = [0.2, 1];
opt.ShowBox = 'off';
opt.XMinorTick = 'on';
opt.YMinorTick = 'on';
opt.FontName = 'Helvetica';
opt.FontSize = 16;
% opt.Colors = cmap;
setAxesProp(opt);
hLines(1).Marker = '.';
hLines(2).Marker = '.';
hLines(1).MarkerSize = 16;
hLines(2).MarkerSize = 16;
hLines(1).LineWidth = 2;
hLines(2).LineWidth = 2;
hLines(1).Color = cmap(2, :);
hLines(2).Color = cmap(1, :);
legend('off');
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', {'0.2', '1'})
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.1, 0.5];
set(gca, 'Position', axesPosition);
title('Behavior', 'FontSize', 16)
% Simulation
axes; hold on;
notBoxPlot(medFreqGainData, medFreqGainLabel, ...
    'jitter', 0.02, 'alpha', 0.125);
xlim([0.5, 2.5]);
ylim([1, 3.7]);
set(gca, 'YTick', 1:3);
set(gca, 'XTick', [1, 2]);
set(gca, 'XTickLabel', ...
    {'\color[rgb]{0.1961,0.7059,0.2902}Strong Signal', ...
    '\color[rgb]{0.6353,0,0}Weak Signal'});
set(gca, 'FontName', 'Helvetica');
ylabel('Gain');
hLine = findobj(gca,'Type','line');
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 16;
setAxesProp(opt);
legend('off');
hLine(1).LineWidth = 4;
hLine(1).Color = cmap(2, :);
hLine(2).LineWidth = 4;
hLine(2).Color = cmap(1, :);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.7, 0.5];
set(gca, 'Position', axesPosition);
title('Simulation', 'FontSize', 16)
% Test statistics (non-parametric ANOVA)
p = kruskalwallis(medFreqGainData, medFreqGainLabel, 'off');
fprintf('Statistics (moth simulation): Kruskal-wallis test - p = %.4f (n = %d)\n', ...
    p, length(medFreqGainData));

% Behavior
medGainData = [medGain{1}; medGain{3}];
medGainLabel = [ones(length(medGain{1}),1); 2*ones(length(medGain{3}),1)];
axes; hold on;
notBoxPlot(medGainData, medGainLabel, 'jitter', 0.05, 'alpha', 0.3);
xlim([0.5, 2.5]);
ylim([1, 3.7]);
set(gca, 'YTick', 1:3);
set(gca, 'XTick', [1, 2]);
set(gca, 'XTickLabel', ...
    {'\color[rgb]{0.1961,0.7059,0.2902}Strong Signal', ...
    '\color[rgb]{0.6353,0,0}Weak Signal'});
set(gca, 'FontName', 'Helvetica');
ylabel('Gain');
hLine = findobj(gca,'Type','line');
opt = [];
opt.BoxDim = [8, 5]*0.35;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 16;
setAxesProp(opt);
legend('off');
hLine(1).LineWidth = 4;
hLine(1).Color = cmap(2, :);
hLine(2).LineWidth = 4;
hLine(2).Color = cmap(1, :);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.5, 0.5];
set(gca, 'Position', axesPosition);
title('Behavior', 'FontSize', 16)
print(gcf,'-dpdf',GEN_SAVE_PATH('fig7M-P.pdf'));
% Test statistics (non-parametric ANOVA)
p = kruskalwallis(medGainData, medGainLabel, 'off');
fprintf('Statistics (moth behavior): Kruskal-wallis test - p = %.4f (n = %d)\n', ...
    p, length(medGainData));
fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));

function out = parseCockroachData(dat)
out.head = dat(:, 1:2);
out.tail = dat(:, 3:4);
out.src = dat(1, 5:6);


function x = rotate2D(x, theta)
offset = x(1, :);
x = x - offset;
rot2x = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
for i = 1:size(x, 1)
    x(i, :) = x(i, :) * rot2x;
end
x = x + offset;

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