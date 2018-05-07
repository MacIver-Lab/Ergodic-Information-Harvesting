function snrData = WiggleAttenuationAnalysis(varargin)
%% Analyze the necessity of wiggling observed in animals and EEDI simulations
% 
% Chen Chen
% 11/22/2017

%% Load dataset
if nargin > 0
    dataFiles = varargin{1};
    %     sourcePath = varargin{1};
else
    sourcePath = './SimDataSource/SimulatedAttenuateTraj/*.mat';
    dataFiles = dir(sourcePath);
end
nFiles = length(dataFiles);

%% Analyze wiggle attenuation
snrData.SNR = [];
snrData.attenGain = [];
snrData.estObjPos = [];
snrData.varBelief = [];
snrData.estError = [];

snrData.meanVarBelief = [];
snrData.rmsVarBelief = [];
snrData.varEstError = [];
snrData.rmsEstError = [];
snrData.meanEstError = [];
for idx = 1:nFiles
    % Load simulated data
    simData = load([dataFiles(idx).folder, '/', dataFiles(idx).name]);
    
    % Evaluate performance
    Perf = mCalcPerf(simData.pB, simData.phi, simData.oTrajList);
     
    fprintf('File %s, gAttn = %d, meanVarBelief = %.4e, rmsVarBelief = %.4e\n', ...
        dataFiles(idx).name, simData.AttenuateMetrics.attenGain, ...
        Perf.meanVarBelief, Perf.rmsVarBelief);
    
    if isempty(find(simData.SNR == [snrData.SNR], 1))
        if isempty([snrData.SNR])
            snrData(1).SNR = simData.SNR;
        else
            snrData(length(snrData)+1).SNR = simData.SNR;
        end
    end
    
    snrData(simData.SNR == [snrData.SNR]).attenGain = [...
        snrData(simData.SNR == [snrData.SNR]).attenGain, simData.AttenuateMetrics.attenGain];
    snrData(simData.SNR == [snrData.SNR]).estObjPos = [...
        snrData(simData.SNR == [snrData.SNR]).estObjPos, Perf.estObjPos];
    snrData(simData.SNR == [snrData.SNR]).varBelief = [...
        snrData(simData.SNR == [snrData.SNR]).varBelief, Perf.varBelief];
    snrData(simData.SNR == [snrData.SNR]).estError = [...
        snrData(simData.SNR == [snrData.SNR]).estError, Perf.estError];

    snrData(simData.SNR == [snrData.SNR]).meanVarBelief = [...
        snrData(simData.SNR == [snrData.SNR]).meanVarBelief, Perf.meanVarBelief];
    snrData(simData.SNR == [snrData.SNR]).rmsVarBelief = [...
        snrData(simData.SNR == [snrData.SNR]).rmsVarBelief, Perf.rmsVarBelief];
    snrData(simData.SNR == [snrData.SNR]).varEstError = [...
        snrData(simData.SNR == [snrData.SNR]).varEstError, Perf.varEstError];
    snrData(simData.SNR == [snrData.SNR]).rmsEstError = [...
        snrData(simData.SNR == [snrData.SNR]).rmsEstError, Perf.rmsEstError];
    snrData(simData.SNR == [snrData.SNR]).meanEstError = [...
        snrData(simData.SNR == [snrData.SNR]).meanEstError, Perf.meanEstError];
end

% Sort and organize data fileds
for i = 1:length([snrData.SNR])
    [snrData(i).attenGain, idx] = sort(snrData(i).attenGain);
    snrData(i).estObjPos = snrData(i).estObjPos(idx);
    snrData(i).varBelief = snrData(i).varBelief(idx);
    snrData(i).estError = snrData(i).estError(idx);
    
    snrData(i).meanVarBelief = snrData(i).meanVarBelief(idx);
    snrData(i).rmsVarBelief = snrData(i).rmsVarBelief(idx);
    snrData(i).varEstError = snrData(i).varEstError(idx);
    snrData(i).rmsEstError = snrData(i).rmsEstError(idx);
    snrData(i).meanEstError = snrData(i).meanEstError(idx);
end


function Perf = mCalcPerf(pB, phi, objPos)
SampleGrid = linspace(0,1,size(phi,1))';
idx = 1;
for eidIter = 2:size(phi,3)
    for simIter = 1:size(phi,2)
        % Max belief estimate
        [~, maxIdx] = max(pB(:, simIter, eidIter));
        Perf.maxBelief(idx) = SampleGrid(maxIdx);
        
        % Mean belief estimate
        Perf.estObjPos(idx) = mean(sum(pB(:, simIter, eidIter) .* SampleGrid));
        
        % Variance of belief
        Perf.varBelief(idx) = var(pB(:, simIter, eidIter));
        
        idx = idx + 1;
    end
end

% Estimation error
Perf.estError = abs(Perf.estObjPos - objPos);
Perf.varEstError = var(Perf.estError);
Perf.rmsEstError = rms(Perf.estError);
Perf.meanEstError = mean(Perf.estError);

Perf.meanVarBelief = mean(Perf.varBelief);
Perf.rmsVarBelief = rms(Perf.varBelief);