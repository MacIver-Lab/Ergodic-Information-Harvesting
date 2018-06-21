function makeSMFig3Plots(dataPath, savePath)
%% Plot supplement figure 3
% Chen Chen

warning('off', 'MATLAB:print:FigureTooLargeForPage');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
GEN_DATA_PATH = @(fname) fullfile(dataPath, fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
%% Load Source Trajectory

load(GEN_DATA_PATH('fig2-ErgodicHarvest-ElectricFish-SNR-30.mat'), 'sTrajList', 'dt');
sPosRaw = sTrajList;

%% Attenuate the Wiggle Frequencies
freqWin = [0.1, 1.5];
attenSamp = [50, 150];
for i = 1:length(attenSamp)
    [so, am] = AdaptiveWiggleAttenuator(sPosRaw, dt, freqWin, attenSamp(i), 0, 0);
    AttenuateMetrics{i} = am;
    sPosOut{i} = so;
    sPos = sPosRaw;
end

%% Plot
figure(1); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [20 10]);
lineColors = [247 148 29; ...
              11  146 70; ...
              43  57  144] / 255.0;
hdl = line([0, freqWin(2)], [0, 0], ...
    'LineWidth', 2, 'Color', lineColors(1, :)); hold on;
hLine = line([freqWin(1), freqWin(1)], [-150, 0], ...
    'LineStyle', '--', 'LineWidth', 2, ...
    'Color', [236 29 39]/255.0);
% Second segment - attenuated freq band
for i = 1:length(attenSamp)
    hdl = [hdl, plot(AttenuateMetrics{i}.magResponseFreqIdx(AttenuateMetrics{i}.cutoffIdx:end), ...
        AttenuateMetrics{i}.magResponsedB(AttenuateMetrics{i}.cutoffIdx:end), ...
        'LineWidth', 2, ...
        'Color', lineColors(i+1, :))];
end

opt = [];
opt.BoxDim = [8,5];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XLabel = 'Frequency (Hz)';
opt.YLabel = 'Magnitude (dB)';
opt.XLim = [0, freqWin(2)];
opt.XTick = [0, freqWin(1), 0.5:0.5:2];
opt.YLim = [-150, 2];
opt.YTick = -150:50:0;
opt.FontName = 'Helvetica';
setPlotProp(opt);
hLine.LineStyle = '--';
ytickangle(45);
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.3, 0.3];
set(gca, 'Position', axesPosition);
print(GEN_SAVE_PATH('sm-fig3a-MagResponse.pdf'), '-dpdf');

figure(2); clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [20 10]);
simTimestamp = (1:length(sPos)) * dt;
plot(simTimestamp, sPos, 'LineWidth',2, 'Color', lineColors(1, :)); hold on;
for i = 1:length(attenSamp)
    plot(simTimestamp, sPosOut{i}, 'LineWidth', 2, 'Color', lineColors(i+1, :));
end
xlabel('Time');
ylabel('Position');
set(gca,'XLim',[0,50]);
set(gca,'YLim',[0.1,0.9]);
opt = [];
opt.BoxDim = [8,5];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XLabel = 'Time (Seconds)';
opt.YLabel = 'Position';
opt.YTick = [];
opt.FontName = 'Helvetica';
setPlotProp(opt);
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.3, 0.3];
set(gca, 'Position', axesPosition);
print(GEN_SAVE_PATH('sm-fig3b-Traj.pdf'), '-dpdf');

figure(3); clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [20 10]);
plot(AttenuateMetrics{1}.FreqResolution:AttenuateMetrics{1}.FreqResolution:freqWin(2), ...
    abs(AttenuateMetrics{1}.sPosFFT(2:freqWin(2)/AttenuateMetrics{1}.FreqResolution+1)),...
    'LineWidth', 2, 'Color', lineColors(1, :)); hold on;
for i = 1:length(attenSamp)
    plot(AttenuateMetrics{i}.FreqResolution:AttenuateMetrics{i}.FreqResolution:freqWin(2), ...
        abs(AttenuateMetrics{i}.sPosFiltFFT(2:freqWin(2)/AttenuateMetrics{i}.FreqResolution+1)),...
        'LineWidth', 2, 'Color', lineColors(i+1, :));
end
hLine = line([freqWin(1), freqWin(1)], [0, 180], ...
    'LineStyle', '--', 'LineWidth', 2, 'Color', [236 29 39]/255.0);
xlabel('Frequency (Hz)');
set(gca,'YTick',[]);
set(gca,'XTick',[0, freqWin(1), 0.5:0.5:freqWin(2)]);
ylabel('Normalized Gain');
opt = [];
opt.BoxDim = [8,5];
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XLabel = 'Frequency (Hz)';
opt.YLabel = 'Magnitude';
opt.YLim = [0, 175];
opt.XLim = [0, 1];
opt.FontName = 'Helvetica';
setPlotProp(opt);
hLine.LineStyle = '--';
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.3, 0.3];
set(gca, 'Position', axesPosition);
print(GEN_SAVE_PATH('sm-fig3c-FFT.pdf'), '-dpdf');

fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));