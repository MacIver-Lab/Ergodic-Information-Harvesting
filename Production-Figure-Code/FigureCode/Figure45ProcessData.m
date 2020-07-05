function Figure45ProcessData(dataPath, savePath)
%% Process simulated data for figure 4-5

GEN_DATA_PATH = @(fname) fullfile(dataPath, 'wiggle_attenuation_sim', 'rawdata', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);

%% STEP #1 - Performance Evaluation
fprintf('********* STEP #1 - Performance Evaluation *********\n');
load('./FigureCode/figure45DataAssociation.mat');
trackingSimFileNamePrefixList = {attenTrialData.trackingSimFileNamePrefix};
nRawSimFiles = length(attenTrialData);
% Tracking simulation data
for idx = 1:nRawSimFiles
    trackingSimFileName = dir(GEN_DATA_PATH([trackingSimFileNamePrefixList{idx}, '*']));
    
    if ~isempty(trackingSimFileName)
        attenTrialData(idx).trackingSimFileName = trackingSimFileName;
        attenTrialData(idx).simData = load(GEN_DATA_PATH(trackingSimFileName.name), ...
            'pB', 'phi', 'ergList', 'eidList', 'enpList', 'sTrajList', 'oTrajList', ...
            'dt');
        attenTrialData(idx).Perf = mCalcPerf(...
            attenTrialData(idx).simData.pB, ...
            attenTrialData(idx).simData.phi, ...
            attenTrialData(idx).simData.ergList, ...
            attenTrialData(idx).simData.enpList, ...
            attenTrialData(idx).srcObjPos);
        attenTrialData(idx).Perf.reEffort = mCalcEffortCalibrated(...
            attenTrialData(idx).simData.oTrajList, ...
            attenTrialData(idx).simData.sTrajList, ...
            attenTrialData(idx).simData.dt);
        fprintf('File %s processed\n...', trackingSimFileName.name);
    else
        warning('Could not locate any file with prefix %s', trackingSimFileNamePrefixList{idx});
    end
end

nTrackingSimTrials = length(attenTrialData);

% Remove unnecessary fields
attenTrialData = rmfield(attenTrialData, 'srcObjPos');
attenTrialData = rmfield(attenTrialData, 'trackingSimFileNamePrefix');

fprintf('\n--------- STEP #1 Success, %d trials analyzed ---------\n', nTrackingSimTrials);

%% STEP #2 - Group and Analyze Data
fprintf('********* STEP #2 - Group and Analyze Data *********\n');

% Preallocate data structure
PerfData.gAtten = [];
PerfData.SNR = [];
PerfData.Rc = [];
PerfData.varEstErrorMaxBelief = [];
PerfData.rmsEstErrorMaxBelief = [];
PerfData.meanEstErrorMaxBelief = [];
PerfData.varEstErrorMeanBelief = [];
PerfData.rmsEstErrorMeanBelief = [];
PerfData.meanEstErrorMeanBelief = [];
PerfData.meanVarBelief = [];
PerfData.rmsVarBelief = [];
PerfData.meanErgodicity = [];
PerfData.rmsErgodicity = [];
PerfData.varErgodicity = [];

gAtten = [attenTrialData.gAtten];
SNR = attenTrialData(end).SNR;
Rc = attenTrialData(end).Rc;
perfDataField = [attenTrialData.Perf];
varEstErrorMaxBelief = [perfDataField.varEstErrorMaxBelief];
rmsEstErrorMaxBelief = [perfDataField.rmsEstErrorMaxBelief];
meanEstErrorMaxBelief = [perfDataField.meanEstErrorMaxBelief];
varEstErrorMeanBelief = [perfDataField.varEstErrorMeanBelief];
rmsEstErrorMeanBelief = [perfDataField.rmsEstErrorMeanBelief];
meanEstErrorMeanBelief = [perfDataField.meanEstErrorMeanBelief];
meanVarBelief = [perfDataField.meanVarBelief];
rmsVarBelief = [perfDataField.rmsVarBelief];
reEffort = [perfDataField.reEffort];
% Sort gain
[gAtten, sortIdx] = sort(gAtten);
PerfData.gAtten = gAtten;
PerfData.SNR = SNR;
PerfData.Rc = Rc;
PerfData.varEstErrorMaxBelief = varEstErrorMaxBelief(sortIdx);
PerfData.rmsEstErrorMaxBelief = rmsEstErrorMaxBelief(sortIdx);
PerfData.meanEstErrorMaxBelief = meanEstErrorMaxBelief(sortIdx);
PerfData.varEstErrorMeanBelief = varEstErrorMeanBelief(sortIdx);
PerfData.rmsEstErrorMeanBelief = rmsEstErrorMeanBelief(sortIdx);
PerfData.meanEstErrorMeanBelief = meanEstErrorMeanBelief(sortIdx);
PerfData.meanVarBelief = meanVarBelief(sortIdx);
PerfData.rmsVarBelief = rmsVarBelief(sortIdx);
PerfData.reEffort = reEffort(sortIdx);
% Ergodicity
meanErgodicity = [perfDataField.meanErgodicity];
rmsErgodicity = [perfDataField.rmsErgodicity];
varErgodicity = [perfDataField.varErgodicity];
PerfData.meanErgodicity = meanErgodicity(sortIdx);
PerfData.rmsErgodicity = rmsErgodicity(sortIdx);
PerfData.varErgodicity = varErgodicity(sortIdx);
% Entropy
meanEntropy = [perfDataField.meanEntropy];
rmsEntropy = [perfDataField.rmsEntropy];
varEntropy = [perfDataField.varEntropy];
PerfData.meanEntropy = meanEntropy(sortIdx);
PerfData.rmsEntropy = rmsEntropy(sortIdx);
PerfData.varEntropy = varEntropy(sortIdx);

fprintf('--------- STEP #2 Success, %d trials analyzed ---------\n', nTrackingSimTrials);

%% Clean up and save data
% Save data
fprintf('All done, saving data... ');
attenTrialData = rmfield(attenTrialData, 'simData');
save(GEN_DATA_PATH('../TrajData.mat'), 'attenTrialData', '-v7.3');
save(GEN_DATA_PATH('../PerfData.mat'), 'PerfData', '-v7.3');
fprintf('Success\n');


function relativeEffort = mCalcEffortCalibrated(oTraj, sTraj, dt)
%% assuming the same spatial and temporal unit of the actual fish
fs = 1.0 / dt;
% refuge's peak to peak length in [cm].
refugeSineLengthP2P = 3.312;
[~, locs] = findpeaks(oTraj, 'NPeaks', 3);
targetLengthP2P = mean(oTraj(locs));
calibratedLengthRatio = refugeSineLengthP2P / targetLengthP2P;

% Fish body weight
% This is estimated by using black ghost knifefish as a reference
% According to MacI10a Energy-Information Trade-Offs between Movement and Sensing
% black ghost knifefish weights 23 grams with a body length of 19 cm.
% Added mass is taken from Clai08a (cited in the paper).
addedMass = 6.04;
fishWeightGram = 23 + addedMass;

%% Approximate velocity through time delta
% Unit [cm/sec]
v_oTraj = diff(oTraj) * fs * calibratedLengthRatio;
v_sTraj = diff(sTraj) * fs * calibratedLengthRatio;
% fprintf('Maximum fish velocity       = %.4f cm/sec\n', max(abs(fVel)));

%% Approximate acceleration through time delta
% Unit [cm/sec^2]
a_oTraj = [0, diff(v_oTraj)] * fs;
a_sTraj = [0, diff(v_sTraj)] * fs;
% fprintf('Maximum fish acceleration   = %.4f cm/sec^2\n', max(abs(fAcc)));

%% Estimate mechanical force for locomotion
% Force unit [g*cm/sec^2] = [1e-5 Newton]
F_oTraj = fishWeightGram * a_oTraj;
F_sTraj = fishWeightGram * a_sTraj;

%% Estimate Mechanical Effort
% Compute instant power, unit [g*cm^2/sec^3] = [1e-7 W]
P_oTraj = abs(F_oTraj .* v_oTraj);
P_sTraj = abs(F_sTraj .* v_sTraj);

%% Integrate to get final mechanical effort estimation
% unit [g*cm^2/sec^2] = [1e-7 J]
time_oTraj = dt:dt:(length(F_oTraj) * dt);
time_sTraj = dt:dt:(length(F_sTraj) * dt);
W_oTraj = trapz(time_oTraj, P_oTraj);
W_sTraj = trapz(time_sTraj, P_sTraj);
% Convert unit to micro-J
% unit 1e-1 * [g*cm^2/sec^2] = 1e-1 * [uJ]
W_oTraj = 1e-1 * W_oTraj;
W_sTraj = 1e-1 * W_sTraj;
relativeEffort = W_sTraj / W_oTraj;


function Perf = mCalcPerf(pB, phi, erg, enp, objPos)
SampleGrid = linspace(0,1,size(phi,1))';
erg = erg(1:end-1); % Remove the last ergodicity value, as it's not valid
enp = enp(2:end);   % Remove the first entropy value, as it's a constant (entropy of a uniform prior belief)

idx = 1;
for eidIter = 2:size(phi,3)
    for simIter = 1:size(phi,2)
        % Max belief estimate
        [~, maxIdx] = max(pB(:, simIter, eidIter));
        Perf.estObjPosMaxBelief(idx) = SampleGrid(maxIdx);

        % Mean belief estimate
        Perf.estObjPosMeanBelief(idx) = mean(sum(pB(:, simIter, eidIter) .* SampleGrid));

        % Variance of belief
        Perf.varBelief(idx) = var(pB(:, simIter, eidIter));
        
        idx = idx + 1;
    end
end

% Estimation error
% Max belief
Perf.estErrorMaxBelief = abs(Perf.estObjPosMaxBelief - objPos);
Perf.varEstErrorMaxBelief = var(Perf.estErrorMaxBelief);
Perf.rmsEstErrorMaxBelief = rms(Perf.estErrorMaxBelief);
Perf.meanEstErrorMaxBelief = mean(Perf.estErrorMaxBelief);
% Mean belief
Perf.estErrorMeanBelief = abs(Perf.estObjPosMeanBelief - objPos);
Perf.varEstErrorMeanBelief = var(Perf.estErrorMeanBelief);
Perf.rmsEstErrorMeanBelief = rms(Perf.estErrorMeanBelief);
Perf.meanEstErrorMeanBelief = mean(Perf.estErrorMeanBelief);

% Variance of belief
Perf.meanVarBelief = mean(Perf.varBelief);
Perf.rmsVarBelief = rms(Perf.varBelief);

% Ergodicicity
Perf.meanErgodicity = mean(erg);
Perf.rmsErgodicity = rms(erg);
Perf.varErgodicity = var(erg);

% Entropy
Perf.meanEntropy = mean(enp);
Perf.rmsEntropy = rms(enp);
Perf.varEntropy = var(enp);