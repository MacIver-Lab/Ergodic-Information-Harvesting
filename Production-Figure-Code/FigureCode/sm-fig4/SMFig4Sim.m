function SMFig4Sim(srcPath, nThreads)
%% Batch simulation for supplement figure 4
% Chen Chen

%% Global Constants
% Control Parameters
GEN_PYTHON_JOB_LIST = 1;

% Directory
DIR_ROOT = fullfile(srcPath, 'SimData/');
DIR_TRIAL_PREFIX = 'SimTrial-';
DIR_TRIAL_ATTEN_TRAJ_PREFIX = 'AttenTraj/';
DIR_TRIAL_SIM_TRAJ_PREFIX = 'SimTraj/';
DIR_FIND_MATFILE = @(d)dir([d,'*.mat']);
DIR_GET_FULL_PATH = @(l,idx)[l(idx).folder, '/', l(idx).name];

% Filter Parameters
Filter.AttenuateGainList = [5,10:20:150];
Filter.AttenuateGainSamps = 8;
Filter.AttenFreqWin = [0.15, 1.5];
Filter.RawTrajLPFHighCut = 1.75;

% EIH Parameters
EIH.MaxTime = 50;

%% Initialization
% Utility lambda function to replace directory separator to backslash
REPDIR = @(d) strrep(d, '\', '/');

% Determine Actual Folder Name
dataFolders = dir([DIR_ROOT,DIR_TRIAL_PREFIX,'*']);
dataFolders = dataFolders([dataFolders.isdir]);
nTrials = 1;%length(dataFolders);
DIR_TRIAL_ROOT = REPDIR([DIR_ROOT, DIR_TRIAL_PREFIX, sprintf('%03d/',nTrials)]);
DIR_TRIAL_ATTEN_TRAJ = REPDIR([DIR_TRIAL_ROOT, DIR_TRIAL_ATTEN_TRAJ_PREFIX]);
DIR_TRIAL_SIM_TRAJ = REPDIR([DIR_TRIAL_ROOT, DIR_TRIAL_SIM_TRAJ_PREFIX]);

% Python Simulation API
PYTHON_CODE_PATH = './FigureCode/sm-fig4/EIH_Sim_Main.py';
PYTHON_JOB_LIST = REPDIR(fullfile(DIR_TRIAL_ROOT, 'SimJobList.txt'));
EIH_API_SIM_TRACKING = @()system(sprintf('python %s %s %d', ...
    PYTHON_CODE_PATH,PYTHON_JOB_LIST,nThreads));
EIH_API_SIM_TRACKING_JOB_STRING = @(rs, src, ofile)sprintf('trackingOnly %d %s %s %s', ...
    rs,src,ofile,DIR_TRIAL_SIM_TRAJ);

% Make new data folders
mkdir(DIR_TRIAL_ATTEN_TRAJ);
mkdir(DIR_TRIAL_SIM_TRAJ);


%% STEP #1 - Generate Attenuated Trajectories
fprintf('********* STEP #1 - Generate Attenuated Trajectories *********\n');
rawSimFileNameList = dir(fullfile(srcPath, 'fig2-ErgodicHarvest-ElectricFish-SNR-30*'));
nRawSimFiles = length(rawSimFileNameList);

nAttenSimFiles = nRawSimFiles * Filter.AttenuateGainSamps * length(Filter.AttenuateGainList);
attenTrialData(nAttenSimFiles).srcName = [];
attenTrialData(nAttenSimFiles).fileName = [];
attenTrialData(nAttenSimFiles).srcObjPos = [];
attenTrialData(nAttenSimFiles).trackingSimFileName = [];
attenTrialData(nAttenSimFiles).gAtten = [];
attenTrialData(nAttenSimFiles).SNR = [];
attenTrialData(nAttenSimFiles).Rc = [];
attenTrialData(nAttenSimFiles).randSeed = [];

for idx = 1:nRawSimFiles
    simData = load(DIR_GET_FULL_PATH(rawSimFileNameList, idx));
    rawTraj = LPF(simData.sTrajList, 1/simData.dt, Filter.RawTrajLPFHighCut);
    
    % Remove unused fileds
    simData = rmfield(simData, 'eidList');
    simData = rmfield(simData, 'timeDone');
    
    for gainIdx = 1:length(Filter.AttenuateGainList)
        % Generate new trajectory
        [simData.sTrajList, simData.AttenuateMetrics] = AdaptiveWiggleAttenuator(...
            rawTraj, simData.dt, Filter.AttenFreqWin, Filter.AttenuateGainList(gainIdx), 0, 0);
        simData.gAtten = Filter.AttenuateGainList(gainIdx);
        
        % Save new trajectory
        fileName = sprintf('EIH-SNR-%2d-gAttn-%d-RawTraj.mat',...
            simData.SNR, Filter.AttenuateGainList(gainIdx));
        save([DIR_TRIAL_ATTEN_TRAJ, fileName], '-struct', 'simData');
        
        for sampIdx = 1:Filter.AttenuateGainSamps
            attenTrialIdx = (idx-1)*nRawSimFiles + (gainIdx-1)*Filter.AttenuateGainSamps + sampIdx;
            
            % Save some additional fields
            attenTrialData(attenTrialIdx).srcName = rawSimFileNameList(idx).name;
            attenTrialData(attenTrialIdx).fileName = fileName;
            attenTrialData(attenTrialIdx).trackingSimFileNamePrefix = sprintf('EIH-SNR-%2d-gAttn-%d-RandSeed-%d',...
                simData.SNR, Filter.AttenuateGainList(gainIdx), sampIdx-1);
            attenTrialData(attenTrialIdx).srcObjPos = simData.oTrajList;
            attenTrialData(attenTrialIdx).gAtten = simData.gAtten;
            attenTrialData(attenTrialIdx).SNR = simData.SNR;
            attenTrialData(attenTrialIdx).Rc = simData.wControl;
            attenTrialData(attenTrialIdx).randSeed = sampIdx-1;
        end
        
        % Move on
        fprintf('File %s created, attenuation gain %d\n', fileName, Filter.AttenuateGainList(gainIdx));
    end
end
nAttenSimFiles = length(attenTrialData);

% Add baseline tracking case
for idx = 1:nRawSimFiles
    simData = load(DIR_GET_FULL_PATH(rawSimFileNameList, idx));
    simData.sTrajList = LPF(simData.sTrajList, 1/simData.dt, Filter.RawTrajLPFHighCut);
    
    % Remove unused fileds
    simData = rmfield(simData, 'eidList');
    simData = rmfield(simData, 'timeDone');
    
    % Add AttenuationMetrics
    simData.AttenuateMetrics.attenGain = 0;

    % Save new trajectory
    fileName = sprintf('EIH-SNR-%2d-gAttn-%d-RawTraj.mat', simData.SNR, 0);
    save([DIR_TRIAL_ATTEN_TRAJ, fileName], '-struct', 'simData');

    for sampIdx = 1:Filter.AttenuateGainSamps
        attenTrialIdx = nAttenSimFiles + (idx-1)*nRawSimFiles*Filter.AttenuateGainSamps + sampIdx;
        
        % Save some additional fields
        attenTrialData(attenTrialIdx).srcName = rawSimFileNameList(idx).name;
        attenTrialData(attenTrialIdx).fileName = fileName;
        attenTrialData(attenTrialIdx).trackingSimFileNamePrefix = sprintf('EIH-SNR-%2d-gAttn-%d-RandSeed-%d',...
            simData.SNR, 0, sampIdx-1);
        attenTrialData(attenTrialIdx).srcObjPos = simData.oTrajList;
        attenTrialData(attenTrialIdx).gAtten = 0;
        attenTrialData(attenTrialIdx).SNR = simData.SNR;
        attenTrialData(attenTrialIdx).Rc = simData.wControl;
        attenTrialData(attenTrialIdx).randSeed = sampIdx-1;
    end
    
    % Move on
    fprintf('File %s created, attenuation gain %d\n', fileName, 0);
end
nAttenSimFiles = length(attenTrialData);

fprintf('--------- STEP #1 Success, %d attenuated trials generated ---------\n', nAttenSimFiles);

%% STEP #2 - Tracking Simulation
fprintf('********* STEP #2 - Tracking Simulation *********\n');
if GEN_PYTHON_JOB_LIST
    fileID = fopen([DIR_TRIAL_ROOT, 'SimJobList.txt'], 'w');
end
parFileList = {attenTrialData.fileName};
parRandSeed = {attenTrialData.randSeed};
parSimFileList = {attenTrialData.trackingSimFileNamePrefix};

for idx = 1:nAttenSimFiles
    if GEN_PYTHON_JOB_LIST
        fileCmd = EIH_API_SIM_TRACKING_JOB_STRING(...
            parRandSeed{idx}, ...
            [DIR_TRIAL_ATTEN_TRAJ, parFileList{idx}], ...
            [parSimFileList{idx}, '.mat']);
        fprintf(fileID, [REPDIR(fileCmd), '\n']);
    end
end
fclose(fileID);
% Submit simulation
fprintf('Submitting simulation with %d attenuated trials using %d CPU threads\n', nAttenSimFiles, nThreads);
disp(['Note that this process could take a LONG time to finish depends ', ...
    'on the number of threads and CPU']);
fprintf('%s\n%s\n%s\n',["And it's normal see the simulation running for a while without", ...
    "any output in the MATLAB command window as MATLAB surpresses them.", ...
    "They will show up eventually"]);
EIH_API_SIM_TRACKING();
fprintf('--------- STEP #2 Success, %d attenuated trials simulated ---------\n', nAttenSimFiles);

%% STEP #3 - Performance Evaluation
fprintf('********* STEP #3 - Performance Evaluation *********\n');
trackingSimFileNamePrefixList = {attenTrialData.trackingSimFileNamePrefix};

% Tracking simulation data
for idx = 1:nAttenSimFiles
    trackingSimFileName = DIR_FIND_MATFILE([DIR_TRIAL_SIM_TRAJ,trackingSimFileNamePrefixList{idx}]);
    trackingSimFileName = trackingSimFileName.name;
    
    if ~isempty(trackingSimFileName)
        attenTrialData(idx).trackingSimFileName = trackingSimFileName;
        attenTrialData(idx).simData = load([DIR_TRIAL_SIM_TRAJ,trackingSimFileName], ...
            'pB', 'phi', 'ergList', 'eidList', 'enpList');
        attenTrialData(idx).Perf = mCalcPerf(...
            attenTrialData(idx).simData.pB, ...
            attenTrialData(idx).simData.phi, ...
            attenTrialData(idx).simData.ergList, ...
            attenTrialData(idx).simData.enpList, ...
            attenTrialData(idx).srcObjPos);
        fprintf('File %s processed\n...', trackingSimFileName);
    else
        warning('Could not locate any file with prefix %s', trackingSimFileNamePrefixList{idx});
    end
end

nTrackingSimTrials = length(attenTrialData);

% Remove unnecessary fields
attenTrialData = rmfield(attenTrialData, 'srcObjPos');
attenTrialData = rmfield(attenTrialData, 'trackingSimFileNamePrefix');

fprintf('\n--------- STEP #3 Success, %d trials analyzed ---------\n', nTrackingSimTrials);

%% STEP #4 - Group and Analyze Data
fprintf('********* STEP #4 - Group and Analyze Data *********\n');

% Preallocate data structure
PerfData(nRawSimFiles).gAtten = [];
PerfData(nRawSimFiles).SNR = [];
PerfData(nRawSimFiles).Rc = [];
PerfData(nRawSimFiles).varEstErrorMaxBelief = [];
PerfData(nRawSimFiles).rmsEstErrorMaxBelief = [];
PerfData(nRawSimFiles).meanEstErrorMaxBelief = [];
PerfData(nRawSimFiles).varEstErrorMeanBelief = [];
PerfData(nRawSimFiles).rmsEstErrorMeanBelief = [];
PerfData(nRawSimFiles).meanEstErrorMeanBelief = [];
PerfData(nRawSimFiles).meanVarBelief = [];
PerfData(nRawSimFiles).rmsVarBelief = [];
PerfData(nRawSimFiles).meanErgodicity = [];
PerfData(nRawSimFiles).rmsErgodicity = [];
PerfData(nRawSimFiles).varErgodicity = [];

srcFileNameList = {attenTrialData.srcName};
for idx = 1:nRawSimFiles
    fileIdx = strfindCell(srcFileNameList, rawSimFileNameList(idx).name);
    gAtten = [attenTrialData(fileIdx).gAtten];
    
    SNR = attenTrialData(fileIdx(end)).SNR;
    Rc = attenTrialData(fileIdx(end)).Rc;
    
    perfDataField = [attenTrialData(fileIdx).Perf];
    
    varEstErrorMaxBelief = [perfDataField.varEstErrorMaxBelief];
    rmsEstErrorMaxBelief = [perfDataField.rmsEstErrorMaxBelief];
    meanEstErrorMaxBelief = [perfDataField.meanEstErrorMaxBelief];
    
    varEstErrorMeanBelief = [perfDataField.varEstErrorMeanBelief];
    rmsEstErrorMeanBelief = [perfDataField.rmsEstErrorMeanBelief];
    meanEstErrorMeanBelief = [perfDataField.meanEstErrorMeanBelief];
    
    meanVarBelief = [perfDataField.meanVarBelief];
    rmsVarBelief = [perfDataField.rmsVarBelief];
        
    % Sort gain
    [gAtten, sortIdx] = sort(gAtten);
    PerfData(idx).gAtten = gAtten;
    
    PerfData(idx).SNR = SNR;
    PerfData(idx).Rc = Rc;
    PerfData(idx).FreqWin = Filter.AttenFreqWin;
    
    PerfData(idx).varEstErrorMaxBelief = varEstErrorMaxBelief(sortIdx);
    PerfData(idx).rmsEstErrorMaxBelief = rmsEstErrorMaxBelief(sortIdx);
    PerfData(idx).meanEstErrorMaxBelief = meanEstErrorMaxBelief(sortIdx);
    
    PerfData(idx).varEstErrorMeanBelief = varEstErrorMeanBelief(sortIdx);
    PerfData(idx).rmsEstErrorMeanBelief = rmsEstErrorMeanBelief(sortIdx);
    PerfData(idx).meanEstErrorMeanBelief = meanEstErrorMeanBelief(sortIdx);
    
    PerfData(idx).meanVarBelief = meanVarBelief(sortIdx);
    PerfData(idx).rmsVarBelief = rmsVarBelief(sortIdx);
    
    % Ergodicity
    meanErgodicity = [perfDataField.meanErgodicity];
    rmsErgodicity = [perfDataField.rmsErgodicity];
    varErgodicity = [perfDataField.varErgodicity];
    PerfData(idx).meanErgodicity = meanErgodicity(sortIdx);
    PerfData(idx).rmsErgodicity = rmsErgodicity(sortIdx);
    PerfData(idx).varErgodicity = varErgodicity(sortIdx);
    
    % Entropy
    meanEntropy = [perfDataField.meanEntropy];
    rmsEntropy = [perfDataField.rmsEntropy];
    varEntropy = [perfDataField.varEntropy];
    PerfData(idx).meanEntropy = meanEntropy(sortIdx);
    PerfData(idx).rmsEntropy = rmsEntropy(sortIdx);
    PerfData(idx).varEntropy = varEntropy(sortIdx);
end

fprintf('--------- STEP #4 Success, %d trials analyzed ---------\n', nTrackingSimTrials);


%% Clean up and save data
FreqWin = Filter.AttenFreqWin;

% Save data
fprintf('All done, saving data... ');
% save(sprintf([DIR_TRIAL_ROOT,'BatchEEDISim-Trial%03d-AllData.mat'],nTrials), ...
%     'attenTrialData', 'FreqWin', '-v7.3');
% Slim file, only save performance data
attenTrialData = rmfield(attenTrialData, 'simData');
% save(sprintf([DIR_TRIAL_ROOT,'BatchEEDISim-Trial%03d-Data.mat'],nTrials), ...
%     'attenTrialData', 'FreqWin', '-v7.3');
% 
% save(sprintf([DIR_TRIAL_ROOT,'BatchEEDISim-Trial%03d-Result.mat'],nTrials), 'PerfData', '-v7.3');
save(fullfile(srcPath, 'sm-fig4-Data.mat'), 'PerfData', '-v7.3');
fprintf('Success\n');

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