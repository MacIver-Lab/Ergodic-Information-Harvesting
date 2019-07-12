function makeFig2Plot(dataPath, savePath)
%% Plot simulated trajectory in figure 2
% Chen Chen

warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_DATA_PATH = @(fname) fullfile(dataPath, 'behavior_sim', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
% Whether nor not to plot EER band
% set to 1 to enable EER plot overlay
% note that due to the complexity of EER bands, it's can be fairly slow 
% to plot the EER bands
global PLOT_EER_BAND
PLOT_EER_BAND = 1;

%% Load data
EH_lSNR = load(GEN_DATA_PATH('EIH-ElectricFish-WeakSignal-RandSeed-1.mat'), 'oTrajList', 'sTrajList', 'phi');
EH_hSNR = load(GEN_DATA_PATH('EIH-ElectricFish-StrongSignal-RandSeed-1.mat'), 'oTrajList', 'sTrajList', 'phi');

%% Make plot
close all;
% weak signal condition
figure(1); clf; hold on;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
plot(EH_lSNR.oTrajList(1:1414), ...
    'LineWidth', 2);
plot(EH_lSNR.sTrajList(1:1414), ...
    'LineWidth', 2);
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [7,5] * 0.6;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [0];
opt.YTick = [0, 1];
opt.XLim = [0, 1414];
opt.YLim = [0, 1];
opt.FontSize = 16;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; 216, 58, 58; 73, 109, 181] / 255.0;
setAxesProp(opt, gca);
legend(gca, 'off');
mPlotContinuousEID(EH_lSNR);
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.4];
set(gca, 'Position', axesPosition);


% strong signal condition
axes; hold on;
plot(EH_hSNR.oTrajList(1:1414), ...
    'LineWidth', 2);
plot(EH_hSNR.sTrajList(1:1414), ...
    'LineWidth', 2);
xlabel('Time');
ylabel('Position');
opt = [];
opt.BoxDim = [7,5] * 0.6;
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off'; 
opt.XTick = [0];
opt.YTick = [0, 1];
opt.XLim = [0, 1414];
opt.YLim = [0, 1];
opt.FontSize = 16;
opt.FontName = 'Helvetica';
opt.Colors = [0, 0, 0; 216, 58, 58; 73, 109, 181] / 255.0;
setAxesProp(opt, gca);
legend(gca, 'off');
mPlotContinuousEID(EH_hSNR);
set(gca, 'units', 'normalized');
axesPosition(1:2) = [0.6, 0.4];
set(gca, 'Position', axesPosition);

print(GEN_SAVE_PATH('fig2-EIH-intro.pdf'), '-dpdf');
fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));


function mPlotContinuousEID(dat)
global PLOT_EER_BAND
if ~PLOT_EER_BAND
    return;
end
%% Plot Parameters
tScale = 10;   % Interval of EID plot update, set to 1 will plot all of the EID map
nBins = 40;   % Color resolution in the y axis
alpha = 0.5;  % Transparency of the EID color
cmap = [0.7 0 0.4];

%% Prepare data
eidList = flattenResultList(dat.phi(:,:,1:end-1))';
tRes = length(dat.oTrajList) / (size(eidList,2)-1);
sRes = size(eidList,1);
s = 1 / sRes;
faces = 1:4;
idxList = tScale:tScale:floor(length(dat.oTrajList) / tRes);

%% Plot
for idx = 1:length(idxList)
    i = idxList(idx);
    if i > 1501
        drawnow;
        return;
    end
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
    drawnow;
end

function outList = flattenResultList(list)
    outList = zeros(size(list,2)*size(list,3), size(list,1));
    for i = 1:size(list,3)
        for j = 1:size(list,2)
            outList((i-1)*size(list,2) + j,:) = list(:,j,i)';
        end
    end