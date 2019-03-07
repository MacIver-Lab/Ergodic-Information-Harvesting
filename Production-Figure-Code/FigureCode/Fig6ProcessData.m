function Fig6ProcessData(dataPath, savePath)
%% Process simulated data for figure 6
% Chen Chen

GEN_DATA_PATH = @(fname) fullfile(dataPath, 'wiggle_attenuation_sim', 'rawdata', fname);
GEN_SAVE_PATH = @(fname) fullfile(savePath, fname);

%% Global Constants
% Directory
% DIR_ROOT = fullfile(srcPath, 'SimData/');
% DIR_TRIAL_PREFIX = 'SimTrial-';
% DIR_TRIAL_ATTEN_TRAJ_PREFIX = 'AttenTraj/';
% DIR_TRIAL_SIM_TRAJ_PREFIX = 'SimTraj/';
% DIR_FIND_MATFILE = @(d)dir([d,'*.mat']);
% DIR_GET_FULL_PATH = @(l,idx)[l(idx).folder, '/', l(idx).name];
% DIR_TRIAL_ATTEN_TRAJ = REPDIR([DIR_TRIAL_ROOT, DIR_TRIAL_ATTEN_TRAJ_PREFIX]);
% DIR_TRIAL_SIM_TRAJ = REPDIR([DIR_TRIAL_ROOT, DIR_TRIAL_SIM_TRAJ_PREFIX]);


%% STEP #1 - Performance Evaluation
fprintf('********* STEP #1 - Performance Evaluation *********\n');
load('./FigureCode/fig6DataAssociation.mat');
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
        oEffort = mCalcEffort(attenTrialData(idx).simData.oTrajList, attenTrialData(idx).simData.dt);
        sEffort = mCalcEffort(attenTrialData(idx).simData.sTrajList, attenTrialData(idx).simData.dt);
        attenTrialData(idx).Perf.reEffort = sEffort / oEffort;
        fprintf('File %s processed\n...', trackingSimFileName.name);
    else
        warning('Could not locate any file with prefix %s', trackingSimFileNamePrefixList{idx});
    end
end

nTrackingSimTrials = length(attenTrialData);

% Remove unnecessary fields
attenTrialData = rmfield(attenTrialData, 'srcObjPos');
attenTrialData = rmfield(attenTrialData, 'trackingSimFileNamePrefix');

fprintf('\n--------- STEP #3 Success, %d trials analyzed ---------\n', nTrackingSimTrials);

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

fprintf('--------- STEP #4 Success, %d trials analyzed ---------\n', nTrackingSimTrials);


%% Clean up and save data
% Save data
fprintf('All done, saving data... ');
attenTrialData = rmfield(attenTrialData, 'simData');
save(GEN_DATA_PATH('../TrajData.mat'), ...
    'attenTrialData', '-v7.3');
save(GEN_DATA_PATH('../PerfData.mat'), 'PerfData', '-v7.3');
fprintf('Success\n');


function w = mCalcEffort(x, dt)
% ignoring the units since we only need the ratio eventually
v = diff(x);
a = [0, diff(v)]; % acceleration
f = abs(a .* v);  % net force
time = dt:dt:dt*(length(x)-1);
w = trapz(time, f);


function idx = strfindCell(cellIn, str)
idx = find(contains(cellIn, str));

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