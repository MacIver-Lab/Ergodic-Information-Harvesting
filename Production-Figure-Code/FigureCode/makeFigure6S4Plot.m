function makeFigure6S4Plot(dataPath, savePath)
%% Plot figure 6---figure supplement 4
% Analyze sensitivity simulations on R/Lambda ratio
% R/Lambda ratio controls the trade-off between energy and ergodicity

warning('off', 'MATLAB:print:FigureTooLargeForPage');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:Axes:NegativeDataInLogAxis');
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
lambda = 5.0;
wsDataRoverLambda = [wsData.R] / lambda;
wsDataRE = [wsData.RE];
figure(1); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [30, 30]);
notBoxPlot(wsDataRE, wsDataRoverLambda, ...
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
opt.YTick = 0:4;
opt.XMinorTick = 'on';
opt.YMinorTick = 'off';
opt.FontName = 'Helvetica';
opt.FontSize = 12;
opt.XLabel = 'R/\lambda';
opt.YLabel = 'Relative exploration';
opt.IgnoreLines = true;
setAxesProp(opt);
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.4, 0.4];
set(gca, 'Position', axesPosition);
print(gcf,'-dpdf',GEN_SAVE_PATH('fig6s4.pdf'));

avgRE_RoverLambda_1 = mean(wsDataRE( abs(wsDataRoverLambda - 1.0)  < 1e-3));
avgRE_RoverLambda_10 = mean(wsDataRE(abs(wsDataRoverLambda - 10.0) < 1e-3));
fprintf(...
    'Average relative exploration for R over lambda == 1 is  %.2f\n', ...
    avgRE_RoverLambda_1);
fprintf(...
    'Average relative exploration for R over lambda == 10 is %.2f\n', ...
    avgRE_RoverLambda_10);