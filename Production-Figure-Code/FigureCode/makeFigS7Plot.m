function makeFigS7Plot(dataPath, savePath)

close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage');
GEN_BEHAVIOR_DATA_PATH = @(fname) fullfile(dataPath, 'animal_behavior_data', 'ElectricFish', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
%% Load data
load(GEN_BEHAVIOR_DATA_PATH('JammingData.mat'));

%% Make plot
figure(1); clf;
set(gcf, ...
    'units','normalized','outerposition',[0 0 1 1], ...
    'PaperPositionMode','auto', ...
    'PaperOrientation','landscape', ...
    'PaperSize', [13 8]);
eodVar = [jarCorrData.mean_MFD];
plot(eodVar,'o-','LineWidth',2,'MarkerSize',12);

% change settings
opt = [];
opt.BoxDim = [8, 5] * 0.5;
opt.LineWidth = 2;
opt.Markers = {'o'};
opt.XLabel = 'Jamming Amplitude (mA)'; % xlabel
opt.YLabel = 'Maximum Frequency Shift'; %ylabel
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XTick = 1:5;
opt.YLimit = [0,6];
opt.FontName = 'Helvetica';
opt.FontSize = 14;
% apply the settings
setAxesProp(opt, gca);
legend('off');
set(gca,'XTickLabel',num2str(sortedJamAmp'));
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.6, 0.5];
set(gca, 'Position', axesPosition);


%% Make EOD plot
load(GEN_BEHAVIOR_DATA_PATH('FishEOD_Jamming.mat'),'Data');

% Process Data
fFishRaw = Data.ElectricFieldData.FishFreq;
fJamRaw = Data.ElectricFieldData.JammingFreq;
FishTimeStamp = Data.ElectricFieldData.FishFreqTimeStamps;
JamTimeStamp = Data.ElectricFieldData.JammingFreqTimeStamps;
JamAmp = Data.ElectricFieldData.JammingAmps;

% Filter Data
tsFish = mean(diff(FishTimeStamp));
tsJam = mean(diff(JamTimeStamp));
fFishLPF = LPF(fFishRaw, 1/tsFish, 0.1);
fJamLPF = LPF(fJamRaw, 1/tsJam, 0.1);
JamSwitchTime = floor(JamTimeStamp(abs(diff(JamAmp))~=0) - 1);
if isempty(JamSwitchTime)
    return;
end
pJTS = [find(JamTimeStamp > JamSwitchTime(1),1),...
    find(JamTimeStamp >= JamSwitchTime(2),1)]-1;

axes;
plot(JamTimeStamp(pJTS(1):pJTS(2)),fJamLPF(pJTS(1):pJTS(2)), 'LineWidth', 3);
hold on;
plot(FishTimeStamp,fFishLPF, 'LineWidth', 3);

% Shading
pVerts = [JamTimeStamp(pJTS(1)),min([fJamLPF,fFishLPF])-2;...
    JamTimeStamp(pJTS(1)),max([fJamLPF,fFishLPF])+2;
    JamTimeStamp(pJTS(2)),max([fJamLPF,fFishLPF])+2;...
    JamTimeStamp(pJTS(2)),min([fJamLPF,fFishLPF])-2;];
pFaces = [1,2,3,4];
patchinfo.Vertices = pVerts;
patchinfo.Faces = pFaces;
patchinfo.FaceColor = 'c';
patchinfo.FaceAlpha = 0.1;
patch(patchinfo);

% Adjust Plot
opt = [];
opt.BoxDim = [8, 5] * 0.5;
opt.LineWidth = [3, 3];
opt.XLabel = 'Time (s)';
opt.YLabel = 'Frequency (Hz)';
opt.ShowBox = 'off';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.XLim = [0,max([FishTimeStamp,JamTimeStamp])+1];
opt.YLim = [min([fJamLPF,fFishLPF])-0.5,...
    max([fJamLPF,fFishLPF])+1];
opt.YTick = [462, 465, 468, 471];
opt.YLim = [461.5, 471];
opt.FontName = 'Helvetica';
opt.FontSize = 14;
opt.Colors = [64/255, 0, 202/255;...
    179/255, 48/255, 0];
% apply the settings
setAxesProp(opt, gca);
legend('off');
set(gca, 'units', 'normalized');
axesPosition = get(gca, 'Position');
axesPosition(1:2) = [0.2, 0.5];
set(gca, 'Position', axesPosition);
print(GEN_SAVE_PATH('figS7.pdf'), '-dpdf');