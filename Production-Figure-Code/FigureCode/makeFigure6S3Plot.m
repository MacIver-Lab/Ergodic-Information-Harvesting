function makeFigure6S3Plot(dataPath, savePath)
%% Plot figure 6---figure supplement 3

close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
GEN_DATA_PATH = @(fname) fullfile(dataPath, 'search_sim', fname);
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile('./PublishedData/', 'animal_behavior_data', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
sTrajHighCutFreq = 2.10; % Hz
global PLOT_EER_BAND
PLOT_EER_BAND = 1;

%% Re-engage searching
% Load and Process Data
% Ergodic Harvesting data
dat = load(GEN_DATA_PATH('EIH-Rat-SearchSim.mat'), ...
    'dt', 'oTrajList', 'sTrajList', 'phi');
sTraj = dat.sTrajList(1:end);
oTraj = dat.oTrajList(1:end);
dat.eidList = flattenResultList(dat.phi(:,:,1:end-1))';
dt = dat.dt;
blindIdx = [810, 1500];
% Filter Sensor Trajectory
sTraj = LPF(sTraj, 1/dt, sTrajHighCutFreq);

% Rat behavior data
ratDat = load(GEN_BEHAVIOR_DATA_PATH('Rat/Khan12a-fig1c.mat'));
ratDat.rat = ratDat.rat(17:end-7, :);
ratDat.odor = ratDat.odor(5:end-12, :);

% Plot
figure(1), clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
% Part #1 - Ergodic Harvesting Search
hold on;
hS = plot(sTraj, 'LineWidth', 2, ...
    'Color', [238, 46, 47]/255.0);
hO = plot(oTraj, 'LineWidth', 2, ...
    'Color', [57, 83, 164]/255.0);
rectangle('Position',[blindIdx(1) 0 diff(blindIdx) 1], ...
    'LineWidth', 4, 'EdgeColor', [15, 128, 64]/255.0);
set(gca,'YLim',[0, 1]);
% Layered EER Patch
mPlotContinuousEID(dat);
% Prettify Figure
opt = [];
opt.BoxDim = [8,4] * 0.5;
opt.YLim = [0, 1];
opt.XLim = [130, length(sTraj)];
opt.XTick = [];
opt.YTick = [];
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.ShowBox = 'off';
opt.FontName = 'Helvetica';
setAxesProp(opt, gca);
legend(gca, 'off');
hS.Color = [238, 46, 47]/255.0;
hO.Color = [57, 83, 164]/255.0;
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.4, 0.2];
set(gca, 'Position', axesPosition);

% Part #2 - Rat Search
h = axes; hold on;
hS = plot(ratDat.rat(:, 1), ratDat.rat(:, 2), 'LineWidth', 2, ...
    'Color', [238, 46, 47]/255.0);
hO = plot(ratDat.odor(:, 1), ratDat.odor(:, 2), 'LineWidth', 2, ...
    'Color', [57, 83, 164]/255.0);
rectangle('Position',[179, mean(ratDat.rat(:, 2))-120, 75.8, 240], ...
    'LineWidth', 4, 'EdgeColor', [15, 136, 64]/255.0);
% Prettify Figure
opt = [];
opt.BoxDim = [8,4] * 0.5;
opt.YLim = [mean(ratDat.rat(:, 2))-120, mean(ratDat.rat(:, 2))+120];
opt.XLim = [min(ratDat.rat(:, 1)), max(ratDat.rat(:, 1))-10];
opt.XTick = [];
opt.YTick = [];
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.ShowBox = 'off';
opt.FontName = 'Helvetica';
setAxesProp(opt, h);
hS.Color = [238, 46, 47]/255.0;
hO.Color = [57, 83, 164]/255.0;
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.4, 0.5];
set(gca, 'Position', axesPosition);

print(GEN_SAVE_PATH('fig6s3A.pdf'), '-dpdf');

%% Initial searching
% Load and Process Data
% Ergodic Harvesting data
dat = load(GEN_DATA_PATH('EIH-Mole-SearchSim.mat'), ...
    'dt', 'oTrajList', 'sTrajList', 'phi');
sTraj = dat.sTrajList(1:1860);
oTraj = dat.oTrajList(1:1860);
dat.eidList = flattenResultList(dat.phi(:,:,1:end-1))';
dt = dat.dt;
% Filter Sensor Trajectory
sTraj = LPF(sTraj, 1/dt, sTrajHighCutFreq);

% Rat behavior data
moleDat = load(GEN_BEHAVIOR_DATA_PATH('Mole/fig2-rBlock-green.mat'), 'angleData');
moleDat.angleData = moleDat.angleData;

% Plot
figure(2), clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
% Part #1 - Ergodic Harvesting Search
hold on;
blindIdx = [0, 1120];
hS = plot(sTraj, 'LineWidth', 2, ...
    'Color', [238, 46, 47]/255.0);
hO = plot(oTraj, 'LineWidth', 2, ...
    'Color', [57, 83, 164]/255.0);
rectangle('Position',[blindIdx(1) 0 diff(blindIdx) 1], ...
    'LineWidth', 4, 'EdgeColor', [15, 128, 64]/255.0);
set(gca,'YLim',[0, 1]);
% Layered EID Patch
mPlotContinuousEID(dat);
% Prettify Figure
opt = [];
opt.BoxDim = [8,4] * 0.5;
opt.YLim = [0, 1];
opt.XLim = [0, length(sTraj)];
opt.XTick = [];
opt.YTick = [];
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.ShowBox = 'off';
opt.FontName = 'Helvetica';
setAxesProp(opt, gca);
legend(gca, 'off');
hS.Color = [238, 46, 47]/255.0;
hO.Color = [57, 83, 164]/255.0;
hO.LineStyle = '--';
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.4, 0.2];
set(gca, 'Position', axesPosition);

% Part #2 - Mole Search
h = axes; hold on;
hS = plot(moleDat.angleData, 'LineWidth', 2, ...
    'Color', [238, 46, 47]/255.0);
hO = line([0, length(moleDat.angleData)], [0, 0], 'LineWidth', 2, ...
    'Color', [57, 83, 164]/255.0);
rectangle('Position',[0, -80, 50, 180], ...
    'LineWidth', 4, 'EdgeColor', [15, 136, 64]/255.0);
% Prettify Figure
opt = [];
opt.BoxDim = [8,4] * 0.5;
opt.YLim = [-80, 100];
opt.XLim = [0, 83];
opt.XTick = [];
opt.YTick = [];
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.ShowBox = 'off';
opt.FontName = 'Helvetica';
setAxesProp(opt, h);
hS.Color = [238, 46, 47]/255.0;
hO.Color = [57, 83, 164]/255.0;
hO.LineStyle = '--';
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.4, 0.5];
set(gca, 'Position', axesPosition);

print(GEN_SAVE_PATH('fig6s3B.pdf'), '-dpdf');

function mPlotContinuousEID(dat)
global PLOT_EER_BAND
if ~PLOT_EER_BAND
    return;
end
%% Plot Parameters
tScale = 5;   % Interval of EID plot update, set to 1 will plot all of the EID map
nBins = 256;   % Color resolution in the y axis
alpha = 0.5;  % Transparency of the EID color
cmap = [0.7 0 0.4];
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


function outList = flattenResultList(list)
outList = zeros(size(list,2)*size(list,3), size(list,1));
for i = 1:size(list,3)
    for j = 1:size(list,2)
        outList((i-1)*size(list,2) + j,:) = list(:,j,i)';
    end
end