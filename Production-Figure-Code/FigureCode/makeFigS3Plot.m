function makeFigS3Plot(dataPath, savePath, usePublished)
%% Plot supplement figure 3
% Chen Chen

close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
GEN_DATA_PATH = @(fname) fullfile(dataPath, fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);

%% Load data
% load('./Data/sm-fig1-EH_IF_Data.mat');
if usePublished
    load(GEN_DATA_PATH('snr_sweep_sim/SNR_SweepSim_Data.mat'));
else
    % Search for simulated data in the working directory
    [snrErg, snrInf, RE_Erg, RE_Inf] = ...
        FigS3ProcessData(GEN_DATA_PATH('snr_sweep_sim/'), './');
end

%% Plot result
figure(1); clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperUnits', 'inches', ...
    'PaperSize', [13 8]);
notBoxPlot(RE_Inf, snrInf, 'plotRawData', false, 'blendPatchColor', true)
hold on;
notBoxPlot(RE_Erg, snrErg, ...
    'plotColor', 'b', 'plotRawData', false, 'blendPatchColor', true);
xlabel('SNR');
ylabel('Relative Exploration');
% title('Relative Exploration vs. SNR');
baseLine = line([5, 56], [1, 1], 'LineStyle', '--', 'LineWidth', 2);
hPatch = findobj(gca,'Type','patch');
legend([hPatch(1), hPatch(end), baseLine], ...
    {'Ergodic Harvesting', 'Infotaxis', 'Baseline'});
opt = [];
opt.BoxDim = [8,5] * 0.5;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [10:10:50];
opt.YTick = [0:0.5:2.5];
opt.XLim = [9, 57];
opt.FontName = 'Helvetica';
opt.FontSize = 14;
opt.IgnoreLines = 1;
setAxesProp(opt, gca);
baseLine.LineStyle = '--';
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.5];
set(gca, 'Position', axesPosition);

% Compute correlation coefficient and its 95% confidence interval
[Rerg, ~, RLerg, RUerg] = corrcoef(double(snrErg), RE_Erg);
[Rinf, ~, RLinf, RUinf] = corrcoef(double(snrInf), RE_Inf);
axes; hold on;
errorbar(Rerg(2),1,RLerg(2)-Rerg(2),RUerg(2)-Rerg(2), ...
    'Horizontal', '.', 'LineWidth', 4, 'MarkerSize', 50, ...
    'Color', 'k', 'CapSize', 20); 
errorbar(Rinf(2),2,RLinf(2)-Rinf(2),RUinf(2)-Rinf(2), ...
    'Horizontal', '.', 'LineWidth', 4, 'MarkerSize', 50, ...
    'Color', 'k', 'CapSize', 20);
ylim([0.5, 2.5]);
set(gca, 'YTick', [1,2]);
set(gca, 'YTickLabel', {'Ergodic Harvesting', 'Infotaxis'});
xlabel('Correlation Coefficient');
opt = [];
opt.BoxDim = [8,5] * 0.5;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [-1:0.2:0.4];
opt.XLim = [-1, 0.4];
% opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.FontSize = 14;
opt.IgnoreLines = 1;
setAxesProp(opt, gca);
legend(gca, 'off');
ytickangle(90);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.6, 0.5];
set(gca, 'Position', axesPosition);
print(GEN_SAVE_PATH('figS3.pdf'), '-dpdf');
fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));