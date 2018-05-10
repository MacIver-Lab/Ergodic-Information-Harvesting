function [snrErg, snrInf, RE_Erg, RE_Inf] = SMFig1ProcessData(dataPath, savePath)

warning('off', 'MATLAB:print:FigureTooLargeForPage');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
GEN_DATA_PATH = @(fname) fullfile(dataPath, fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);
%% Load data
eihTrials = dir(GEN_DATA_PATH('*ErgodicHarvest*.mat'));
ifTrials = dir(GEN_DATA_PATH('*Infotaxis*.mat'));
if isempty(eihTrials) || isempty(ifTrials)
    error(['Cannot found simulation data, please make sure sm-fig1 ', ...
           'simulation has finished or use USE_PREV_DATASET = 1']);
end

%% Process data
[ifData, ~] = mProcessData(ifTrials);
[eihData, ~] = mProcessData(eihTrials);

%% Plot result
snrErg = double(eihData.SNR);
% eihData.SNR(mod(eihData.SNR, 5) ~= 0) = eihData.SNR(mod(eihData.SNR, 5) ~= 0) - 0.5;
snrInf = double(ifData.SNR);
RE_Erg = eihData.RE;
RE_Inf = ifData.RE;

save(GEN_SAVE_PATH('sm-fig1-EH_IF_Data.mat'), ...
    'snrErg', 'snrInf', 'gainRatioErg', 'gainRatioInf');

function Perf = mEvalPerf(dat)
% Estimate target position using mean belief
if size(dat.pB, 2) == 1
    % InfoMax
    dat.pB(:,1,1) = dat.pB(:,1,1) / sum(dat.pB(:,1,1));
    posEstimate = squeeze([sum(dat.pB .* linspace(0,1,101)')])';  
else
    % Ergodic Harvesting
    posEstimate = [];
    for i = 1:size(dat.pB, 3)
        posEstimate = [posEstimate, squeeze([sum(dat.pB(:,:,i) .* linspace(0,1,101)')])];
    end
end

errEstimate = abs(posEstimate(1:length(dat.oTrajList)) - dat.oTrajList);
absPosErr = dat.sTrajList(1:length(dat.oTrajList)) - dat.oTrajList;
Perf.rmsEstErr = rms(errEstimate);
Perf.meanEstErr = mean(errEstimate);
Perf.varEstErr = var(errEstimate);
binEdgesErr = linspace(-0.4, 0.4, 21);
binEdges = linspace(0.2, 0.8, 21);
[Perf.posErrHist.N, Perf.posErrHist.edges] = histcounts(absPosErr, binEdgesErr);
[Perf.posHistSensor.N, Perf.posHistSensor.edges] = histcounts(dat.sTrajList, binEdges);
[Perf.posHistTarget.N, Perf.posHistTarget.edges] = histcounts(dat.oTrajList, binEdges);

% FFT
dat.sTrajList = movmean(dat.sTrajList, 5);
dat.sTrajList = LPF(dat.sTrajList, 1/dat.dt, 2);
Perf.RE = sum(abs(diff(dat.sTrajList))) / sum(abs(diff(dat.oTrajList)));


function [snrData, snrStruct] = mProcessData(trials)
snrData.SNR = [];
snrData.rmsEstErr = [];
snrData.meanEstErr = [];
snrData.varEstErr = [];
snrData.meanEntropy = [];
snrData.rmsEntropy = [];
snrData.varEntropy = [];
snrData.posErrHist = [];
snrData.posErrHistData = [];
snrData.posHistSensor = [];
snrData.posHistSensorData = [];
snrData.posHistTarget = [];
snrData.posHistTargetData = [];
snrData.RE = [];
for i = 1:length(trials)
    dat = load([trials(i).folder, '/', trials(i).name]);
    
    % Evaluate tracking performance
    Perf = mEvalPerf(dat);
    
    % Append data
    snrData.SNR = [snrData.SNR, dat.SNR];
    snrData.rmsEstErr = [snrData.rmsEstErr, Perf.rmsEstErr];
    snrData.meanEstErr = [snrData.meanEstErr, Perf.meanEstErr];
    snrData.varEstErr = [snrData.varEstErr, Perf.varEstErr];
    snrData.posErrHist = [snrData.posErrHist, Perf.posErrHist];
    snrData.posErrHistData = [snrData.posErrHistData; Perf.posErrHist.N];
    snrData.posHistSensor = [snrData.posHistSensor, Perf.posHistSensor];
    snrData.posHistSensorData = [snrData.posHistSensorData; Perf.posHistSensor.N];
    snrData.posHistTarget = [snrData.posHistTarget, Perf.posHistTarget];
    snrData.posHistTargetData = [snrData.posHistTargetData; Perf.posHistTarget.N];
    snrData.RE = [snrData.RE, Perf.RE];
    
    if size(dat.pB, 2) == 1
        % InfoMax trial
        snrData.meanEntropy = [snrData.meanEntropy, mean(dat.enpList(1:end-1))];
        snrData.rmsEntropy = [snrData.rmsEntropy, rms(dat.enpList(1:end-1))];
        snrData.varEntropy = [snrData.varEntropy, var(dat.enpList(1:end-1))];
    else
        % ErgInfo trial
        snrData.meanEntropy = [snrData.meanEntropy, mean(dat.enpList(2:end))];
        snrData.rmsEntropy = [snrData.rmsEntropy, rms(dat.enpList(2:end))];
        snrData.varEntropy = [snrData.varEntropy, var(dat.enpList(2:end))];    
    end
end

% Summarize SNR conditions if there are more than 1 samples per SNR
% condition
SNRsamps = unique(snrData.SNR);
nSamps = length(SNRsamps);
snrStruct(nSamps).SNR = [];
snrStruct(nSamps).rmsEstErr = [];
snrStruct(nSamps).meanEstErr = [];
snrStruct(nSamps).varEstErr = [];
snrStruct(nSamps).meanEntropy = [];
snrStruct(nSamps).rmsEntropy = [];
snrStruct(nSamps).varEntropy = [];
snrStruct(nSamps).posErrHist = [];
snrStruct(nSamps).posErrHistData = [];
snrStruct(nSamps).posErrHistDataMean = [];
snrStruct(nSamps).posHistSensorData = [];
snrStruct(nSamps).posHistSensorDataMean = [];
snrStruct(nSamps).posHistTargetData = [];
snrStruct(nSamps).posHistTargetDataMean = [];
snrStruct(nSamps).RE = [];

for i = 1:nSamps
    snrStruct(i).SNR = SNRsamps(i);
    snrStruct(i).rmsEstErr = snrData.rmsEstErr(snrData.SNR == SNRsamps(i));
    snrStruct(i).meanEstErr = snrData.meanEstErr(snrData.SNR == SNRsamps(i));
    snrStruct(i).varEstErr = snrData.varEstErr(snrData.SNR == SNRsamps(i));
    snrStruct(i).meanEntropy = snrData.meanEntropy(snrData.SNR == SNRsamps(i));
    snrStruct(i).rmsEntropy = snrData.rmsEntropy(snrData.SNR == SNRsamps(i));
    snrStruct(i).varEntropy = snrData.varEntropy(snrData.SNR == SNRsamps(i));
    snrStruct(i).posErrHist = snrData.posErrHist(snrData.SNR == SNRsamps(i));
    snrStruct(i).posErrHistData = snrData.posErrHistData(snrData.SNR == SNRsamps(i),:);
    snrStruct(i).RE = snrData.RE(snrData.SNR == SNRsamps(i));
    if size(snrData.posErrHistData(snrData.SNR == SNRsamps(i),:),1) > 1
        snrStruct(i).posErrHistDataMean = mean(snrData.posErrHistData(snrData.SNR == SNRsamps(i),:));
    else
        snrStruct(i).posErrHistDataMean = snrData.posErrHistData(snrData.SNR == SNRsamps(i),:);
    end
    if size(snrData.posHistSensorData(snrData.SNR == SNRsamps(i),:),1) > 1
        snrStruct(i).posHistSensorDataMean = mean(snrData.posHistSensorData(snrData.SNR == SNRsamps(i),:));
    else
        snrStruct(i).posHistSensorDataMean = snrData.posHistSensorData(snrData.SNR == SNRsamps(i),:);
    end
    if size(snrData.posHistTargetData(snrData.SNR == SNRsamps(i),:),1) > 1
        snrStruct(i).posHistTargetDataMean = mean(snrData.posHistTargetData(snrData.SNR == SNRsamps(i),:));
    else
        snrStruct(i).posHistTargetDataMean = snrData.posHistTargetData(snrData.SNR == SNRsamps(i),:);
    end

end
