function makeSMFig1Plot(dataPath, savePath, usePrevDat)
%% Plot supplement figure 1
% Chen Chen

warning('off', 'MATLAB:print:FigureTooLargeForPage');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
GEN_DATA_PATH = @(fname) fullfile(dataPath, fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);

%% Load data
if usePrevDat
    load(GEN_DATA_PATH('sm-fig1-EH_IF_Data.mat'));
else
    % Search for simulated data in the working directory
    [snrErg, snrInf, RE_Erg, RE_Inf] = ...
        SMFig1ProcessData(dataPath, savePath);
end

%% Plot result
figure(1); clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperUnits', 'inches', ...
    'PaperSize', [13 8]);
notBoxPlot(RE_Inf, snrInf)
hold on;
notBoxPlot(RE_Erg, snrErg, 'plotColor', 'b');
xlabel('SNR');
ylabel('Relative Exploration');
% title('Relative Exploration vs. SNR');
baseLine = line([5, 56], [1, 1], 'LineStyle', '--', 'LineWidth', 2);
hPatch = findobj(gca,'Type','patch');
legend([hPatch(1), hPatch(end), baseLine], ...
    {'Ergodic Harvesting', 'Infotaxis', 'Baseline'});
opt = [];
opt.BoxDim = [8,5] * 0.6;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [10:10:60];
opt.YTick = [0:0.5:2.5];
opt.XLim = [9, 57];
opt.FontName = 'Helvetica';
setPlotProp(opt);
baseLine.LineStyle = '--';
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.25, 0.2];
set(gca, 'Position', axesPosition);
print(GEN_SAVE_PATH('sm-fig1-RelativeTrackingEffort.pdf'), '-dpdf');

% Compute correlation coefficient and its 95% confidence interval
[Rerg, ~, RLerg, RUerg] = corrcoef(double(snrErg), RE_Erg);
[Rinf, ~, RLinf, RUinf] = corrcoef(double(snrInf), RE_Inf);
figure(2); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperUnits', 'inches', ...
    'PaperSize', [13 8]);
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
opt.BoxDim = [8,5] * 0.6;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [-1:0.2:0.2];
opt.XLim = [-1, 0.2];
% opt.FontSize = 10;
opt.FontName = 'Helvetica';
setPlotProp(opt);
legend(gca, 'off');
ytickangle(90);
fig3Path = strrep(GEN_SAVE_PATH(''), 'sm-fig1', 'fig3');
mkdir(fig3Path);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.25, 0.2];
set(gca, 'Position', axesPosition);
print(fullfile(fig3Path,'fig3c-CoorelationCoefficient.pdf'), '-dpdf');
fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));