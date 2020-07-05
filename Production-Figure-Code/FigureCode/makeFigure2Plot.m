function makeFigure2Plot(dataPath, savePath)
%% Plot individual panels for figure 2
% 

warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, 'behavior_sim', fname);
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile('./PublishedData/', 'animal_behavior_data', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
barColor = [72, 110, 181;...
    50, 180, 74; ...
    236, 29, 36] / 255;
% set to 1 to enable EER plot overlay
% note that due to the complexity of EER bands, it's can be fairly slow 
% to plot the EER bands
global PLOT_EER_BAND
PLOT_EER_BAND = 1;

%% Load data
EH_lSNR = load(GEN_DATA_PATH('EIH-ElectricFish-WeakSignal-RandSeed-1.mat'), ...
    'oTrajList', 'sTrajList', 'dt', 'phi', 'pB');
EH_hSNR = load(GEN_DATA_PATH('EIH-ElectricFish-StrongSignal-RandSeed-1.mat'), ...
    'oTrajList', 'sTrajList', 'dt', 'phi', 'pB');
EH_lSNR.eidList = flattenResultList(EH_lSNR.phi(:,:,1:end-1))';
EH_hSNR.eidList = flattenResultList(EH_hSNR.phi(:,:,1:end-1))';


% Load fish behavioral data
ss = load(GEN_BEHAVIOR_DATA_PATH('/ElectricFish/Eigenmannia-sp-StrongSignal.mat'), 'strongSigData');
fish.hSNR.fishTraj = ss.strongSigData{8, 1};
fish.hSNR.refugeTraj = ss.strongSigData{8, 2};
dt = 1 / 60;
trajSegStr = 31.5;
trajSegLen = 30;
trajSegStrIdx = trajSegStr * 60;
trajSegIdx = (trajSegStr+trajSegLen) * 60;
timeIdx = 0:dt:trajSegLen*60*dt;
fish.hSNR.fishTraj = fish.hSNR.fishTraj(trajSegStrIdx:trajSegIdx);
fish.hSNR.fishTraj = fish.hSNR.fishTraj - mean(fish.hSNR.fishTraj);
fish.hSNR.refugeTraj = fish.hSNR.refugeTraj(trajSegStrIdx:trajSegIdx);
fish.hSNR.refugeTraj = fish.hSNR.refugeTraj - mean(fish.hSNR.refugeTraj);

ws = load(GEN_BEHAVIOR_DATA_PATH('/ElectricFish/Eigenmannia-sp-WeakSignal.mat'), 'weakSigData');
fish.lSNR.fishTraj = ws.weakSigData{1, 1};
fish.lSNR.fishTraj = fish.lSNR.fishTraj - mean(fish.lSNR.fishTraj);
fish.lSNR.refugeTraj = ws.weakSigData{1, 2};
fish.lSNR.refugeTraj = fish.lSNR.refugeTraj - mean(fish.lSNR.refugeTraj);
% Filter trajectory
simTrajHighCutFreq = 2.10;
EH_lSNR.sTrajList = LPF(EH_lSNR.sTrajList, 1/EH_lSNR.dt, simTrajHighCutFreq);
EH_hSNR.sTrajList = LPF(EH_hSNR.sTrajList, 1/EH_hSNR.dt, simTrajHighCutFreq);

% calibrate distance unit
pix2cm = 1.7 / (max(fish.hSNR.refugeTraj) - min(fish.hSNR.refugeTraj));

%--------- Fish Sinusoidal Tracking ---------%
% Strong Signal Trajectory plot
figure(1);clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
plot(timeIdx, fish.hSNR.refugeTraj*pix2cm, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(timeIdx, fish.hSNR.fishTraj*pix2cm, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time (s)');
ylabel('Position (cm)');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [0, 30];
opt.YTick = -2:2:2;
opt.YLim = [-2, 2];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor(1:2,:);
setAxesProp(opt, gca);
legend(gca, 'off');
set(gca,  'Position', [1    4    2.8320    1.7700]);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.75];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
plot(timeIdx, fish.lSNR.refugeTraj*pix2cm, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(timeIdx, fish.lSNR.fishTraj*pix2cm, ...
    'LineWidth', 2, ...
    'Color', barColor(3,:));
xlabel('Time (s)');
ylabel('Position (cm)');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XTick = [0, 30];
opt.YTick = -2:2:2;
opt.YLim = [-2, 2];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3],:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.75];
set(gca, 'Position', axesPosition);


%--------- Ergodic Harvesting Simulation Trajectory ---------%
%--------------------- EID Plot --------------------------%
% Strong Signal Trajectory plot
trajAxes = axes; hold on;
plot(EH_hSNR.oTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(EH_hSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(EH_hSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor(1:2,:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(EH_hSNR, 0);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.2, 0.2];
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
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(EH_lSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3],:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(EH_lSNR, 0);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.2];
set(gca, 'Position', axesPosition);

%--------------------- Belief Plot --------------------------%
EH_lSNR.eidList = flattenResultList(EH_lSNR.pB(:,:,1:end-1))';
EH_hSNR.eidList = flattenResultList(EH_hSNR.pB(:,:,1:end-1))';
% Strong Signal Trajectory plot
trajAxes = axes; hold on;
plot(EH_hSNR.oTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(EH_hSNR.sTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(2,:));
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(EH_hSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor(1:2,:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(EH_hSNR, 1);
set(gca, 'units', 'normalized');
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
opt.BoxDim = [8,4] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [];
opt.YTick = [];
opt.XLim = [0, length(EH_lSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3],:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(EH_lSNR, 1);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.5];
set(gca, 'Position', axesPosition);
drawnow;
print(GEN_SAVE_PATH('fig2.pdf'), '-dpdf');
fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));


function mPlotContinuousEID(dat, isBelief)
global PLOT_EER_BAND
if ~PLOT_EER_BAND
    return;
end
%% Plot Parameters
tScale = 5;   % Interval of EID plot update, set to 1 will plot all of the EID map
nBins = 256;   % Color resolution in the y axis
alpha = 0.5;  % Transparency of the EID color
if isBelief
    cmap = [0.1 0.4 0.9];
else
    cmap = [0.7 0 0.4];
end
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