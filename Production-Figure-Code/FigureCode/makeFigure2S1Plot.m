function makeFigure2S1Plot(dataPath, savePath)
%% Plot individual panels for figure 2---figure supplement 1
% 
warning('off', 'MATLAB:print:FigureTooLargeForPage');
global GEN_SAVE_PATH
GEN_DATA_PATH = @(fname) fullfile(dataPath, 'distractor_sim', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
% Whether nor not to plot EER band
% set to 1 to enable EER plot overlay
% note that due to the complexity of EER bands, it's can be fairly slow 
% to plot the EER bands

global PLOT_EER_BAND BAR_COLOR
PLOT_EER_BAND = 1;
BAR_COLOR = [72, 110, 181;...
    50, 180, 74; ...
    236, 29, 36] / 255;

%% Cockroach odor source localization
mMakePanel(...
    GEN_DATA_PATH('EIH-DistractorSim.mat'), ...
    'EIH', 1, ...
    'fig2_s1_a1.pdf');
mMakePanel(...
    GEN_DATA_PATH('EIH-DistractorSim.mat'), ...
    'EIH', 0, ...
    'fig2_s1_a2.pdf');
mMakePanel(...
    GEN_DATA_PATH('IF-DistractorSim.mat'), ...
    'Infotaxis', 1, ...
    'fig2_s1_b1.pdf');
mMakePanel(...
    GEN_DATA_PATH('IF-DistractorSim.mat'), ...
    'Infotaxis', 0, ...
    'fig2_s1_b2.pdf');
fprintf('Figure panels created at %s\n', GEN_SAVE_PATH(''));


function mMakePanel(file, alg, showBelief, outFilename)
global GEN_SAVE_PATH
global BAR_COLOR
figure; clf; hold on;
if showBelief
    title(alg + " (Belief)");
    isBelief = 1;
else
    title(alg + " (EID)");
    isBelief = 0;
end
%--------- Ergodic Harvesting Simulation Trajectory ---------%
lsnr = load(file, 'dt', 'sTrajList', 'oTrajList', 'phi', 'pB', 'objCenter');
trajLen = 1600;
lsnr.oTrajList = lsnr.oTrajList(1:trajLen);
lsnr.sTrajList = LPF(lsnr.sTrajList(1:trajLen), 1/lsnr.dt, 2);
if showBelief
    lsnr.eidList = flattenResultList(lsnr.pB(:,:,1:end-1))';
else
    lsnr.eidList = flattenResultList(lsnr.phi(:,:,1:end-1))';
end
lsnr.eidList = lsnr.eidList(:, 1:trajLen);
% Plot
% Reference path
line([0, length(lsnr.oTrajList)], [lsnr.objCenter, lsnr.objCenter], ...
    'LineStyle', '-', 'LineWidth', 4, 'Color', BAR_COLOR(1, :));
if ~contains(file, 'DualTarget')
    line([0, length(lsnr.oTrajList)], [1-lsnr.objCenter, 1-lsnr.objCenter], ...
        'LineStyle', '--', 'LineWidth', 4, 'Color', BAR_COLOR(1, :));
else
    line([0, length(lsnr.oTrajList)], [1-lsnr.objCenter, 1-lsnr.objCenter], ...
        'LineStyle', '-', 'LineWidth', 4, 'Color', BAR_COLOR(1, :));
end
% Weak signal
plot(lsnr.sTrajList, ...
    'LineWidth', 2, ...
    'Color', BAR_COLOR(3, :));
% Configure figure
ylim([0, 1]);
xticks([]);
yticks([]);
xlabel('Time');
ylabel('Lateral Position');
mPlotContinuousEID(lsnr, isBelief);
legend(gca, 'off');
drawnow;
print(gcf, GEN_SAVE_PATH(outFilename), '-dpdf');
   
function mPlotContinuousEID(dat, isBelief)
global PLOT_EER_BAND
if ~PLOT_EER_BAND
    return;
end
%% Plot Parameters
tScale = 10;   % Interval of EID plot update, set to 1 will plot all of the EID map
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
