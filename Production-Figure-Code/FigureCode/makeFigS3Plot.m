function makeFigS3Plot(dataPath, savePath)
%% Analyze sensitivity simulations on R/Lambda ratio
% R/Lambda ratio controls the trade-off between energy and ergodicity
% Chen Chen

close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
GEN_DATA_PATH = @(fname) fullfile(dataPath, fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);

%% Load data
wsFiles = dir(GEN_DATA_PATH('sensitivity_sim/EIH-SA-WeakSignal-*.mat'));
% lambda function for computing cumulative distance traveled
cumDist = @(x) sum(abs(diff(x)));
for i = 1:length(wsFiles)
    ws = load([wsFiles(i).folder, '/', wsFiles(i).name], ...
        'sTrajList', 'oTrajList', 'wControl', 'dt');
    % filter trajectory
    wsTrajTarget = LPF(ws.oTrajList(500:end), 1 / ws.dt, 2.10);
    wsTrajSensor = LPF(ws.sTrajList(500:end), 1 / ws.dt, 2.10);
    
    % Relative exploration
    wsData(i).RE = cumDist(wsTrajSensor) / cumDist(wsTrajTarget);
    % control cost -> R
    wsData(i).R = double(ws.wControl);
    fprintf('Proessing trial %3d (%.1f%%)...\n', i, 100*i/length(wsFiles));
end

%% Make plots
wsDataR = [wsData.R];
wsDataRE = [wsData.RE];
figure(1); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [16 8]);
notBoxPlot(wsDataRE, wsDataR/5, ...
    'jitter', 0, 'jitterScale', 0.1, ...
    'plotRawData', false, 'plotColor', 'b', 'blendPatchColor', true, ...
    'scaledJitter', true);
% mark reference
line([10/5, 10/5], [0, 4], ...
    'LineStyle', '--', 'LineWidth', 2, 'Color', [0.65, 0, 0]);
set(gca, 'XScale', 'log');
opt = [];
opt.BoxDim = [8, 5];
opt.ShowBox = 'off';
opt.XTick = 10.^(-1:2);
opt.XLim = [10^-1.1, 10^2.1];
opt.YLim = [0, 4];
opt.XMinorTick = 'on';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 12;
opt.XLabel = 'R/\lambda';
opt.YLabel = 'Relative exploration';
opt.IgnoreLines = true;
setAxesProp(opt);
legend('off');
print(gcf,'-dpdf',GEN_SAVE_PATH('figS3.pdf'));