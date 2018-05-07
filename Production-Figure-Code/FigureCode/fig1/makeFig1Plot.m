function makeFig1Plot(dataPath, savePath)
%% Plot simulated trajectory in figure 1 d-e
%
% Chen Chen
% 12/19/2017

warning('off', 'MATLAB:print:FigureTooLargeForPage');
%% Load data
GEN_DATA_PATH = @(fname) fullfile(dataPath, fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
EH_lSNR = load(GEN_DATA_PATH('fig1-ErgodicHarvest-SNR-10.mat'), 'oTrajList', 'sTrajList');
EH_hSNR = load(GEN_DATA_PATH('fig1-ErgodicHarvest-SNR-60.mat'), 'oTrajList', 'sTrajList');
IT_lSNR = load(GEN_DATA_PATH('fig1-Infotaxis-SNR-10.mat'), 'oTrajList', 'sTrajList');
IT_hSNR = load(GEN_DATA_PATH('fig1-Infotaxis-SNR-60.mat'), 'oTrajList', 'sTrajList');

%% Make plot
% Panel d - low light condition
figure(1); clf; hold on;
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
set(gca,  'Position', [4, 4.5, 4.8, 3]);

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
set(gca,  'Position', [4, 1, 4.8, 3]);

print(GEN_SAVE_PATH('fig1-de.pdf'), '-dpdf');
