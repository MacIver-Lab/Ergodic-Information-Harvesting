function makeFigure6S2Plot(dataPath, savePath)
%% Plot figure 6---figure supplement 2

close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, 'behavior_sim', fname);
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile('./PublishedData/', 'animal_behavior_data', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
barColor = [72, 110, 181;...
    50, 180, 74; ...
    236, 29, 36] / 255;
% Lambda function handle for computing cumulative 1D distance travelled
cumDist = @(x) sum(abs(diff(x)));
% Whether nor not to plot EER band
% set to 1 to enable EER plot overlay
% note that due to the complexity of EER bands, it's can be fairly slow 
% to plot the EER bands
global PLOT_EER_BAND
PLOT_EER_BAND = 1;

%% Rat Odor Tracking
% Load Data
lSNR = load(GEN_DATA_PATH('EIH-Rat-WeakSignal.mat'));
hSNR = load(GEN_DATA_PATH('EIH-Rat-StrongSignal.mat'));
lSNR.eidList = flattenResultList(lSNR.phi(:,:,1:end-1))';
hSNR.eidList = flattenResultList(hSNR.phi(:,:,1:end-1))';

khan = load(GEN_BEHAVIOR_DATA_PATH('Rat/Khan12a_fig2.mat'));
khan.lSNR.sTraj = khan.fig2b_nose;
khan.lSNR.oTraj = khan.fig2b_trail;
khan.hSNR.sTraj = khan.fig2a_nose;
khan.hSNR.oTraj = khan.fig2a_trail;

% Adjust time horizon to fit into the actual data length (provided in Khan12a)
lSNR.dt = 7.8947 / length(lSNR.sTrajList);
hSNR.dt = 6.8421 / length(hSNR.sTrajList);
khan.lSNR.dt = 7.8947 / length(khan.lSNR.sTraj);
khan.hSNR.dt = 6.8421 / length(khan.hSNR.sTraj);

% Compute cumulative 1D distance travelled
% Exclude initial global search before converge to ensure data consistency
% The criteria is to crop the initial searching trajectory until the sensor
% crosses the target for the first time
dist_hSNR_Sensor = cumDist(hSNR.sTrajList(1:end));
dist_hSNR_Trail = cumDist(hSNR.oTrajList(1:end));
dist_lSNR_Sensor = cumDist(lSNR.sTrajList(1:end));
dist_lSNR_Trail = cumDist(lSNR.oTrajList(1:end));
dist_rat_hSNR_Sensor = cumDist(khan.hSNR.sTraj);
dist_rat_hSNR_Trail = cumDist(khan.hSNR.oTraj);
dist_rat_lSNR_Sensor = cumDist(khan.lSNR.sTraj);
dist_rat_lSNR_Trail = cumDist(khan.lSNR.oTraj);

%--------- Rat Odor Tracking data from Khan12a ---------%
% Strong Signal Trajectory plot
figure(2);clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [30, 30]);
plot(khan.hSNR.oTraj, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(khan.hSNR.sTraj, ...
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
opt.XLim = [0, length(khan.hSNR.oTraj)];
opt.YLim = [mean(khan.hSNR.oTraj)-50, mean(khan.hSNR.oTraj)+50];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor(1:2,:);
setAxesProp(opt, gca);
legend(gca, 'off');
set(gca,  'Position', [1    4    2.8320    1.7700]);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.6];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
plot(khan.lSNR.oTraj, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(khan.lSNR.sTraj, ...
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
opt.XLim = [0, length(khan.lSNR.oTraj)];
opt.YLim = [mean(khan.lSNR.oTraj)-50, mean(khan.lSNR.oTraj)+50];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3],:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.6];
set(gca, 'Position', axesPosition);
% Relative Exploration bar plot
barAxes = axes; hold on;
barData = [1, dist_rat_hSNR_Sensor/dist_rat_hSNR_Trail, dist_rat_lSNR_Sensor/dist_rat_lSNR_Trail];
for i = 1:3
    bar(i, barData(i), 0.4, 'BaseValue', -1, ...
        'FaceColor', barColor(i,:));
end
opt = [];
opt.BoxDim = [8,5] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [1, 2, 3];
opt.YTick = [1, 5];
opt.YLim = [-1, 5];
opt.FontSize = 8;
opt.FontName = 'Helvetica';
setAxesProp(opt, barAxes);
set(gca,'YTickLabel', {'1x', '5x'});
set(gca,'XTickLabel', {'Target', 'Strong Signal', 'Weak Signal'});
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.6];
set(gca, 'Position', axesPosition);

%--------- Ergodic Harvesting Simulation Trajectory ---------%
% Strong Signal Trajectory plot
trajAxes = axes; hold on;
plot(hSNR.oTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(hSNR.sTrajList, ...
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
opt.XLim = [0, length(hSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor(1:2,:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(hSNR);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.2, 0.3];
set(gca, 'Position', axesPosition);
% Weak Signal Trajectory
trajAxes = axes; hold on;
plot(lSNR.oTrajList, ...
    'LineWidth', 2, ...
    'Color', barColor(1,:));
plot(lSNR.sTrajList, ...
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
opt.XLim = [0, length(lSNR.oTrajList)];
opt.YLim = [0, 1];
opt.FontSize = 10;
opt.FontName = 'Helvetica';
opt.Colors = barColor([1,3],:);
setAxesProp(opt, trajAxes);
legend(gca, 'off');
mPlotContinuousEID(lSNR);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.4, 0.3];
set(gca, 'Position', axesPosition);
% Bar plot
barAxes = axes; hold on;
barData = [1, dist_hSNR_Sensor/dist_hSNR_Trail, dist_lSNR_Sensor/dist_lSNR_Trail];
for i = 1:3
    bar(i, barData(i), 0.4, 'BaseValue', -1, ...
        'FaceColor', barColor(i,:));
end
opt = [];
opt.BoxDim = [8,5] * 0.354;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [1, 2, 3];
opt.YTick = [1, 5];
opt.YLim = [-1, 5];
opt.FontSize = 8;
opt.FontName = 'Helvetica';
setAxesProp(opt, barAxes);
set(gca,'YTickLabel', {'1x', '5x'});
set(gca,'XTickLabel', {'Target', 'Strong Signal', 'Weak Signal'});
legend(gca, 'off');
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.3];
set(gca, 'Position', axesPosition);

% All set, now print the first section into PDF
drawnow;
print(gcf, GEN_SAVE_PATH('fig6s2.pdf'),'-dpdf');


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