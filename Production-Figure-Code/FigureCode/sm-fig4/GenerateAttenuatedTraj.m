function attenTrialFileList = GenerateAttenuatedTraj(varargin)
%% Analyze the necessity of wiggling observed in animals and EEDI simulations
% 
% Chen Chen
% 11/22/2017

%% Load dataset
if nargin > 0
    sourcePath = varargin{1};
    savePath = varargin{2};
    AttenuateGainList = varargin{3};
    freqWin = varargin{4}; 
    nAttenTrials = length(AttenuateGainList);
else
    sourcePath = './SimDataSource/*.mat';
    savePath = './SimDataSource/AttenuatedTraj/';
    AttenuateGainList = [5, 10, 20:20:200];
    freqWin = [0.2, 2];
    nAttenTrials = length(AttenuateGainList)+1;
end

dataFiles = dir(sourcePath);

%% Generate Attenuated Trajectories
attenTrialFileList(nAttenTrials).source = [];
attenTrialFileList(nAttenTrials).name = [];
attenTrialFileList(nAttenTrials).path = [];
attenTrialIdx = 1;

for idx = 1:length(dataFiles)
    simData = load([dataFiles(idx).folder, '/', dataFiles(idx).name]);
    rawTraj = simData.sTrajList;
    
    % Add field
    simData.maxTime = ceil(length(simData.oTrajList) * simData.dt);
    
    % Remove unused fileds
    simData = rmfield(simData, 'eidList');
    simData = rmfield(simData, 'timeDone');
    
    for gainIdx = 1:length(AttenuateGainList)
        % Generate new trajectory
        [simData.sTrajList, simData.AttenuateMetrics] = AdaptiveWiggleAttenuator(...
            rawTraj, simData.dt, freqWin, AttenuateGainList(gainIdx), 0);
        
        % Save new trajectory
        fileName = sprintf('EID-SNR-%2d-Sigma-%.3f-oAmp-%.2f-wCtrl-%2d-dt-%.3f-gAttn-%d.mat',...
            simData.SNR, simData.Sigma, simData.objAmp, simData.wControl, simData.dt, AttenuateGainList(gainIdx));
        save([savePath, fileName], '-struct', 'simData');
        attenTrialFileList(attenTrialIdx).source = dataFiles(idx).name;
        attenTrialFileList(attenTrialIdx).path = savePath;
        attenTrialFileList(attenTrialIdx).name = fileName;
        attenTrialIdx = attenTrialIdx + 1;
        fprintf('File %s created, attenuation gain %d\n', fileName, AttenuateGainList(gainIdx));
    end
end