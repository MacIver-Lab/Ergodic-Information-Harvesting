function makeFigure6S1Plot(dataPath, savePath)
%% Plot figure 6---figure supplement 1
% Note that due to the complexity of this figure, each animal's inset will
% be plotted separately into individual PDF files
%
warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, 'behavior_sim', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
barColor = [72, 110, 181;...
    50, 180, 74; ...
    236, 29, 36] / 255;

% Whether nor not to plot EER band
% set to 1 to enable EER plot overlay
% note that due to the complexity of EER bands, it's can be fairly slow
% to plot the EER bands
global PLOT_EER_BAND
PLOT_EER_BAND = 1;

%% Electric Fish Simulation
% Load data
EH_lSNR = load(GEN_DATA_PATH('EIH-ElectricFish-WeakSignal-RandSeed-1.mat'), ...
    'oTrajList', 'sTrajList', 'dt', 'phi', 'pB');
EH_hSNR = load(GEN_DATA_PATH('EIH-ElectricFish-StrongSignal-RandSeed-1.mat'), ...
    'oTrajList', 'sTrajList', 'dt', 'phi', 'pB');
EH_lSNR.eidList = flattenResultList(EH_lSNR.pB(:,:,1:end-1))';
EH_hSNR.eidList = flattenResultList(EH_hSNR.pB(:,:,1:end-1))';
EH_lSNR.eidList = [ones(size(EH_lSNR.eidList,1), 1)/size(EH_lSNR.eidList,1), EH_lSNR.eidList(:, 1:end-1)];
EH_hSNR.eidList = [ones(size(EH_hSNR.eidList,1), 1)/size(EH_hSNR.eidList,1), EH_hSNR.eidList(:, 1:end-1)];

%--------- Fish Sinusoidal Tracking ---------%
% Strong Signal Trajectory plot
figure(1); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8], ...
    'Renderer', 'OpenGL');
%--------- Ergodic Harvesting Simulation Trajectory ---------%
% Strong Signal Trajectory plot
plot(EH_hSNR.oTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(EH_hSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.7;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(EH_hSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 16;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(2,:)];
setAxesProp(opt, gca);
legend(gca, 'off');
mPlotContinuousEID(EH_hSNR);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.5];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
plot(EH_lSNR.oTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(EH_lSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.7;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(EH_lSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 16;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(3,:)];
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(EH_lSNR);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.5, 0.5];
set(gca, 'Position', axesPosition);

% All set, now print the first section into PDF
print(GEN_SAVE_PATH('fig6s1_weakly_electric_fish.pdf'),'-dpdf');

%% Mole Odor Localization
% EIH
hSNR = load(GEN_DATA_PATH('EIH-Mole-StrongSignal-RandSeed-3.mat'));
lSNR = load(GEN_DATA_PATH('EIH-Mole-WeakSignal-RandSeed-3.mat'));
lSNR.eidList = flattenResultList(lSNR.pB(:,:,1:end))';
hSNR.eidList = flattenResultList(hSNR.pB(:,:,1:end))';
lSNR.eidList = [ones(size(lSNR.eidList,1), 1)/size(lSNR.eidList,1), lSNR.eidList(:, 1:end-1)];
hSNR.eidList = [ones(size(hSNR.eidList,1), 1)/size(hSNR.eidList,1), hSNR.eidList(:, 1:end-1)];
hSNR.sTrajList = hSNR.sTrajList(1:1200);
lSNR.sTrajList = lSNR.sTrajList(1:1200);
% Strong Signal Trajectory plot
figure(2); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8], ...
    'Renderer', 'OpenGL');
%--------- Ergodic Harvesting Simulation Trajectory ---------%
% Strong Signal Trajectory plot
hLine = line([0, length(hSNR.sTrajList)], [0.5, 0.5], ...
    'LineStyle', '--', ...
    'LineWidth', 2, ...
    'Color', 'k');
plot(hSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time');
ylabel('Lateral Position');
opt = [];
opt.BoxDim = [8,4] * 0.7;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(hSNR.sTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 16;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(2,:)];
setAxesProp(opt, gca);
legend(gca, 'off');
mPlotContinuousEID(hSNR);
hLine.LineStyle = '--';
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.5];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
hLine = line([0, length(lSNR.sTrajList)], [0.5, 0.5], ...
    'LineStyle', '--', ...
    'LineWidth', 2, ...
    'Color', 'k');
plot(lSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Time');
ylabel('Lateral Position');
opt = [];
opt.BoxDim = [8,4] * 0.7;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(lSNR.sTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 16;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(3,:)];
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(lSNR);
hLine.LineStyle = '--';
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.5, 0.5];
set(gca, 'Position', axesPosition);

% All set, now print the first section into PDF
print(GEN_SAVE_PATH('fig6s1_mole.pdf'),'-dpdf');

%% Cockroach odor source localization
% Plot
figure(3); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8], ...
    'Renderer', 'OpenGL');
%--------- Ergodic Harvesting Simulation Trajectory ---------%
lSNR = load(GEN_DATA_PATH('EIH-Cockroach-WeakSignal-RandSeed-2.mat'), ...
    'dt', 'sTrajList', 'oTrajList', 'phi', 'pB');
hSNR = load(GEN_DATA_PATH('EIH-Cockroach-StrongSignal-RandSeed-2.mat'), ...
    'dt', 'sTrajList', 'oTrajList', 'phi', 'pB');
trajLen = 1400;
hSNR.oTrajList = hSNR.oTrajList(1:trajLen);
hSNR.sTrajList = LPF(hSNR.sTrajList(1:trajLen), 1/hSNR.dt, 2);
lSNR.oTrajList = lSNR.oTrajList(1:trajLen);
lSNR.sTrajList = LPF(lSNR.sTrajList(1:trajLen), 1/lSNR.dt, 2);
hSNR.eidList = flattenResultList(hSNR.pB(:,:,1:end-1))';
lSNR.eidList = flattenResultList(lSNR.pB(:,:,1:end-1))';
hSNR.eidList = hSNR.eidList(:, 1:trajLen);
lSNR.eidList = lSNR.eidList(:, 1:trajLen);
lSNR.eidList = [ones(size(lSNR.eidList,1), 1)/size(lSNR.eidList,1), lSNR.eidList(:, 1:end-1)];
hSNR.eidList = [ones(size(hSNR.eidList,1), 1)/size(hSNR.eidList,1), hSNR.eidList(:, 1:end-1)];

% Plot
% Reference path
line([0, length(hSNR.oTrajList)], [0.5, 0.5], ...
    'LineStyle', '--', 'LineWidth', 4, 'Color', barColor(1, :));
% Strong signal
plot(hSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(2, :));
% Configure figure
xlabel('Time');
ylabel('Lateral Position');
opt = [];
opt.BoxDim = [8,4] * 0.7;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XLim = [0, length(hSNR.oTrajList)];
opt.YLim = [0.2, 0.8];
opt.XTick = [];
opt.YTick = [];
opt.FontSize = 16;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(2,:)];
setAxesProp(opt, gca);
legend(gca, 'off');
mPlotContinuousEID(hSNR);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.5];
set(gca, 'Position', axesPosition);
% Weak Signal
trajAxes = axes; hold on;
% Reference path
line([0, length(hSNR.oTrajList)], [0.5, 0.5], ...
    'LineStyle', '--', 'LineWidth', 4, 'Color', barColor(1, :));
% Weak signal
plot(lSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(3, :));
% Configure figure
xlabel('Time');
ylabel('Lateral Position');
opt = [];
opt.BoxDim = [8,4] * 0.7;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XLim = [0, length(hSNR.oTrajList)];
opt.YLim = [0.2, 0.8];
opt.XTick = [];
opt.YTick = [];
opt.FontSize = 16;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; barColor(3,:)];
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(lSNR);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.5, 0.5];
set(gca, 'Position', axesPosition);

% All set, now print the first section into PDF
print(GEN_SAVE_PATH('fig6s1_cockroach.pdf'),'-dpdf');
fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));

function mPlotContinuousEID(dat)
global PLOT_EER_BAND
if ~PLOT_EER_BAND
    return;
end
%% Plot Parameters
tScale = 5;   % Interval of EID plot update, set to 1 will plot all of the EID map
nBins = 256;  % Color resolution in the y axis
alpha = 0.5;  % Transparency of the EID color
cmap = [0.1 0.4 0.9];
eidList = dat.eidList;
tRes = length(dat.oTrajList) / (size(eidList,2)-1);
sRes = size(eidList,1);
s = 1 / sRes;
faces = 1:4;

idxList = tScale:tScale:floor(length(dat.oTrajList) / tRes);
for idx = 1:length(idxList)
    i = idxList(idx);
    [~,~,bin] = histcounts(eidList(:,i), nBins);
    for k = 1:sRes
        if bin(k) <= 2
            continue;
        end
        verts = [(i-tScale)*tRes, (k-1)*s;...
            (i-0)*tRes, (k-1)*s;...
            (i-0)*tRes, (k-0)*s;...
            (i-tScale)*tRes, (k-0)*s];
        patch('Faces',faces,'Vertices',verts,...
            'FaceColor', cmap,...
            'FaceAlpha', alpha*bin(k)/nBins,...
            'EdgeColor', 'none');
    end
end
for idx = sRes:sRes:size(eidList,2)
    line('XData', [idx, idx], ...
        'YData', [0, 1] , ...
        'LineStyle', '--', 'LineWidth', 0.5);
end



function outList = flattenResultList(list)
outList = zeros(size(list,2)*size(list,3), size(list,1));
for i = 1:size(list,3)
    for j = 1:size(list,2)
        outList((i-1)*size(list,2) + j,:) = list(:,j,i)';
    end
end