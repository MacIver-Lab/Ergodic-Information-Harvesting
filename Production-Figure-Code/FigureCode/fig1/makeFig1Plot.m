function makeFig1Plot(dataPath, savePath)
%% Plot simulated trajectory in figure 1 d-e
% Chen Chen

warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
%% Load data
EH_lSNR = load(GEN_DATA_PATH('fig1-ErgodicHarvest-SNR-10.mat'), 'oTrajList', 'sTrajList');
EH_hSNR = load(GEN_DATA_PATH('fig1-ErgodicHarvest-SNR-60.mat'), 'oTrajList', 'sTrajList');
IT_lSNR = load(GEN_DATA_PATH('fig1-Infotaxis-SNR-10.mat'), 'oTrajList', 'sTrajList');
IT_hSNR = load(GEN_DATA_PATH('fig1-Infotaxis-SNR-60.mat'), 'oTrajList', 'sTrajList');

%% Make plot
% Panel d - low light condition
figure(1); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
plot(EH_lSNR.oTrajList, ...
    'LineWidth', 2);
plot(EH_lSNR.sTrajList, ...
    'LineWidth', 2);
plot(IT_lSNR.sTrajList, ...
    'LineWidth', 2);
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,5] * 0.6;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [0];
opt.YTick = [0, 1];
opt.XLim = [0, length(EH_lSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 16;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; 216, 58, 58; 73, 109, 181] / 255.0;
setPlotProp(opt);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.4, 0.2];
set(gca, 'Position', axesPosition);


% Panel e - high light condition
axes; hold on;
plot(EH_hSNR.oTrajList, ...
    'LineWidth', 2);
plot(EH_hSNR.sTrajList, ...
    'LineWidth', 2);
plot(IT_hSNR.sTrajList, ...
    'LineWidth', 2);
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,5] * 0.6;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [0];
opt.YTick = [0, 1];
opt.XLim = [0, length(EH_hSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 16;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; 216, 58, 58; 73, 109, 181] / 255.0;
setPlotProp(opt);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.55];
set(gca, 'Position', axesPosition);

print(GEN_SAVE_PATH('fig1-de.pdf'), '-dpdf');
fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));