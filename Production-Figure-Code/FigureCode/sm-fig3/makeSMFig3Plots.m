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
% First segment - pass band
hdl = plot(AttenuateMetrics{1}.magResponseFreqIdx(1:AttenuateMetrics{1}.cutoffIdx), ...
    AttenuateMetrics{1}.magResponsedB(1:AttenuateMetrics{1}.cutoffIdx), 'LineWidth', 2); hold on;
hLine = line([freqWin(1), freqWin(1)], [-150, 0], ...
    'LineStyle', '--', 'LineWidth', 2);
% Second segment - attenuated freq band
for i = 1:length(attenSamp)
    hdl = [hdl, plot(AttenuateMetrics{i}.magResponseFreqIdx(AttenuateMetrics{i}.cutoffIdx:end), ...
        AttenuateMetrics{i}.magResponsedB(AttenuateMetrics{i}.cutoffIdx:end), 'LineWidth', 2)];
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
legend(hdl, 'Original', sprintf('Filtered (Attenuation = %ddB)', attenSamp(1)), ...
    sprintf('Filtered (Attenuation = %ddB)', attenSamp(2)),...
    'Location', 'Best');
print(GEN_SAVE_PATH('sm-fig3a-MagResponse.pdf'),'-dpdf');

figure(2); clf;
simTimestamp = (1:length(sPos)) * dt;
plot(simTimestamp, sPos, 'LineWidth',2); hold on;
for i = 1:length(attenSamp)
    plot(simTimestamp, sPosOut{i}, 'LineWidth',2);
end
legend('Original', sprintf('Filtered (Attenuation = %ddB)', attenSamp(1)), ...
    sprintf('Filtered (Attenuation = %ddB)', attenSamp(2)), ...
    'Location', 'Best');
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
print(GEN_SAVE_PATH('sm-fig3b-Traj.pdf'),'-dpdf');

figure(3); clf;
plot(AttenuateMetrics{1}.FreqResolution:AttenuateMetrics{1}.FreqResolution:freqWin(2), ...
    abs(AttenuateMetrics{1}.sPosFFT(2:freqWin(2)/AttenuateMetrics{1}.FreqResolution+1)),...
    'LineWidth',2); hold on;
for i = 1:length(attenSamp)
    plot(AttenuateMetrics{i}.FreqResolution:AttenuateMetrics{i}.FreqResolution:freqWin(2), ...
        abs(AttenuateMetrics{i}.sPosFiltFFT(2:freqWin(2)/AttenuateMetrics{i}.FreqResolution+1)),...
        'LineWidth',2);
end
hLine = line([freqWin(1), freqWin(1)], [0, 180], ...
    'LineStyle', '--', 'LineWidth', 2);
legend('Original', sprintf('Filtered (Attenuation = %ddB)', attenSamp(1)), ...
    sprintf('Filtered (Attenuation = %ddB)', attenSamp(2)), ...
    'Location', 'Best');
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
print(GEN_SAVE_PATH('sm-fig3c-FFT.pdf'),'-dpdf');