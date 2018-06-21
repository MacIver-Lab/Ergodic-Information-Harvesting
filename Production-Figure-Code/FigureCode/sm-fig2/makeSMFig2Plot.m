function makeSMFig2Plot(dataPath, savePath)
%% Plot supplement figure 2
% Chen Chen

close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
GEN_DATA_PATH = @(fname) fullfile(dataPath, fname);
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile(pwd, 'FigureCode', 'sm-fig2', 'BehaviorData', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);

barColor = [72, 110, 181;...
            50, 180, 74; ...
            162, 0, 0] / 255;
dFreq = 0.005;
maxFreq = 2;
fftEndIdx = floor(maxFreq/dFreq) + 1;
freqTicks = dFreq:dFreq:maxFreq;
simTrajHighCutFreq = 2.10;
nFreqSamps = @(t) ceil((dFreq*t)^-1);

%% Electric Fish Simulation
% Load data
EH_lSNR = load(GEN_DATA_PATH('fig2-ErgodicHarvest-ElectricFish-SNR-30.mat'), ...
    'oTrajList', 'sTrajList', 'dt');
EH_hSNR = load(GEN_DATA_PATH('fig2-ErgodicHarvest-ElectricFish-SNR-60.mat'), ...
    'oTrajList', 'sTrajList', 'dt');

% Load fish behavioral data
fish.hSNR = load(GEN_BEHAVIOR_DATA_PATH('ElectricFish-StrongSignal-Sine.mat'));
fish.lSNR = load(GEN_BEHAVIOR_DATA_PATH('ElectricFish-WeakSignal-Sine.mat'));

% Filter trajectory
EH_lSNR.sTrajList = LPF(EH_lSNR.sTrajList, 1/EH_lSNR.dt, simTrajHighCutFreq);
EH_hSNR.sTrajList = LPF(EH_hSNR.sTrajList, 1/EH_hSNR.dt, simTrajHighCutFreq);

% FFT
%   sTraj - Sensor trajectory
%   tTraj - Target trajectory
fft_EH_hSNR_sTraj = fft(EH_hSNR.sTrajList-mean(EH_hSNR.sTrajList), nFreqSamps(EH_hSNR.dt));
fft_EH_hSNR_tTraj = fft(EH_hSNR.oTrajList-mean(EH_hSNR.oTrajList), nFreqSamps(EH_hSNR.dt));
fft_EH_hSNR_sTraj = abs(fft_EH_hSNR_sTraj(2:fftEndIdx));
fft_EH_hSNR_tTraj = abs(fft_EH_hSNR_tTraj(2:fftEndIdx));
fft_EH_lSNR_sTraj = fft(EH_lSNR.sTrajList-mean(EH_lSNR.sTrajList), nFreqSamps(EH_lSNR.dt));
fft_EH_lSNR_tTraj = fft(EH_lSNR.oTrajList-mean(EH_lSNR.oTrajList), nFreqSamps(EH_lSNR.dt));
fft_EH_lSNR_sTraj = abs(fft_EH_lSNR_sTraj(2:fftEndIdx));
fft_EH_lSNR_tTraj = abs(fft_EH_lSNR_tTraj(2:fftEndIdx));
fft_fish_lSNR_sTraj = fft(fish.lSNR.fishTraj-mean(fish.lSNR.fishTraj), nFreqSamps(1/fish.lSNR.FPS));
fft_fish_lSNR_tTraj = fft(fish.lSNR.refugeTraj-mean(fish.lSNR.refugeTraj), nFreqSamps(1/fish.lSNR.FPS));
fft_fish_lSNR_sTraj = abs(fft_fish_lSNR_sTraj(2:fftEndIdx));
fft_fish_lSNR_tTraj = abs(fft_fish_lSNR_tTraj(2:fftEndIdx));
fft_fish_hSNR_sTraj = fft(fish.hSNR.fishTraj-mean(fish.hSNR.fishTraj), nFreqSamps(1/fish.hSNR.FPS));
fft_fish_hSNR_tTraj = fft(fish.hSNR.refugeTraj-mean(fish.hSNR.refugeTraj), nFreqSamps(1/fish.hSNR.FPS));
fft_fish_hSNR_sTraj = abs(fft_fish_hSNR_sTraj(2:fftEndIdx));
fft_fish_hSNR_tTraj = abs(fft_fish_hSNR_tTraj(2:fftEndIdx));

%--------- Fish Sinusoidal Tracking ---------%
figure(1);clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
plot(freqTicks, ...
    fft_fish_lSNR_tTraj/max(fft_fish_lSNR_tTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(freqTicks, ...
    fft_fish_hSNR_sTraj/max(fft_fish_hSNR_sTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
plot(freqTicks, ...
    fft_fish_lSNR_sTraj/max(fft_fish_lSNR_sTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Frequency');
ylabel('FFT Magnitude');
opt = [];
opt.BoxDim = [8,5] * 0.4;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, 1];
opt.FontSize = 12;
opt.FontName = 'Helvetica';
opt.Colors = barColor;
setPlotProp(opt);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.3, 0.3];
set(gca, 'Position', axesPosition);

% Ergodic Harvesting
axes; hold on;
plot(freqTicks, ...
    fft_EH_hSNR_tTraj/max(fft_EH_hSNR_tTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(freqTicks, ...
    fft_EH_hSNR_sTraj/max(fft_EH_hSNR_sTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
plot(freqTicks, ...
    fft_EH_lSNR_sTraj/max(fft_EH_lSNR_sTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Frequency');
ylabel('Fourier Magnitude');
opt = [];
opt.BoxDim = [8,5] * 0.4;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, 0.5];
opt.FontSize = 12;
opt.FontName = 'Helvetica';
opt.Colors = barColor;
setPlotProp(opt);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.3, 0.6];
set(gca, 'Position', axesPosition);


%% Rat Odor Tracking
% Load Data
lSNR = load(GEN_DATA_PATH('fig2-ErgodicHarvest-Rat-WeakSignal.mat'));
hSNR = load(GEN_DATA_PATH('fig2-ErgodicHarvest-Rat-StrongSignal.mat'));

khan = load(GEN_BEHAVIOR_DATA_PATH('Khan12a_fig2.mat'));
khan.lSNR.sTraj = khan.fig2b_nose;
khan.lSNR.oTraj = khan.fig2b_trail;
khan.hSNR.sTraj = khan.fig2a_nose;
khan.hSNR.oTraj = khan.fig2a_trail;

% Adjust time horizon to fit into the actual data length (provided in Khan12a)
lSNR.dt = 7.8947 / length(lSNR.sTrajList);
hSNR.dt = 6.8421 / length(hSNR.sTrajList);
khan.lSNR.dt = 7.8947 / length(khan.lSNR.sTraj);
khan.hSNR.dt = 6.8421 / length(khan.hSNR.sTraj);

% FFT
%   sTraj - Sensor trajectory
%   tTraj - Target trajectory
fft_EH_hSNR_sTraj = fft(hSNR.sTrajList-mean(hSNR.sTrajList), nFreqSamps(hSNR.dt));
fft_EH_hSNR_tTraj = fft(hSNR.oTrajList-mean(hSNR.oTrajList), nFreqSamps(hSNR.dt));
fft_EH_hSNR_sTraj = abs(fft_EH_hSNR_sTraj(2:fftEndIdx));
fft_EH_hSNR_tTraj = abs(fft_EH_hSNR_tTraj(2:fftEndIdx));
fft_EH_lSNR_sTraj = fft(lSNR.sTrajList-mean(lSNR.sTrajList), nFreqSamps(lSNR.dt));
fft_EH_lSNR_tTraj = fft(lSNR.oTrajList-mean(lSNR.oTrajList), nFreqSamps(lSNR.dt));
fft_EH_lSNR_sTraj = abs(fft_EH_lSNR_sTraj(2:fftEndIdx));
fft_EH_lSNR_tTraj = abs(fft_EH_lSNR_tTraj(2:fftEndIdx));
fft_rat_lSNR_sTraj = fft(khan.lSNR.sTraj-mean(khan.lSNR.sTraj), nFreqSamps(khan.lSNR.dt));
fft_rat_lSNR_tTraj = fft(khan.lSNR.oTraj-mean(khan.lSNR.oTraj), nFreqSamps(khan.lSNR.dt));
fft_rat_lSNR_sTraj = abs(fft_rat_lSNR_sTraj(2:fftEndIdx));
fft_rat_lSNR_tTraj = abs(fft_rat_lSNR_tTraj(2:fftEndIdx));
fft_rat_hSNR_sTraj = fft(khan.hSNR.sTraj-mean(khan.hSNR.sTraj), nFreqSamps(khan.hSNR.dt));
fft_rat_hSNR_tTraj = fft(khan.hSNR.oTraj-mean(khan.hSNR.oTraj), nFreqSamps(khan.hSNR.dt));
fft_rat_hSNR_sTraj = abs(fft_rat_hSNR_sTraj(2:fftEndIdx));
fft_rat_hSNR_tTraj = abs(fft_rat_hSNR_tTraj(2:fftEndIdx));

%--------- Rat Odor Tracking data from Khan12a ---------%
% Strong Signal Trajectory plot
axes; hold on;
plot(freqTicks, ...
    fft_rat_lSNR_tTraj/max(fft_rat_lSNR_tTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(freqTicks, ...
    fft_rat_hSNR_sTraj/max(fft_rat_hSNR_sTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
plot(freqTicks, ...
    fft_rat_lSNR_sTraj/max(fft_rat_lSNR_sTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Frequency');
ylabel('FFT Magnitude');
opt = [];
opt.BoxDim = [8,5] * 0.4;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, 1.8];
opt.FontSize = 12;
opt.FontName = 'Helvetica';
opt.Colors = barColor;
setPlotProp(opt);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.5, 0.6];
set(gca, 'Position', axesPosition);

% Ergodic Harvesting
axes; hold on;
plot(freqTicks, ...
    fft_EH_hSNR_tTraj/max(fft_EH_hSNR_tTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(freqTicks, ...
    fft_EH_hSNR_sTraj/max(fft_EH_hSNR_sTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
plot(freqTicks, ...
    fft_EH_lSNR_sTraj/max(fft_EH_lSNR_sTraj), ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Frequency');
ylabel('Fourier Magnitude');
opt = [];
opt.BoxDim = [8,5] * 0.4;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, 1.8];
opt.FontSize = 12;
opt.FontName = 'Helvetica';
opt.Colors = barColor;
setPlotProp(opt);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.5, 0.3];
set(gca, 'Position', axesPosition);

%% Print to file
print(GEN_SAVE_PATH('sm-fig2a-FFT.pdf'),'-dpdf');
fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));