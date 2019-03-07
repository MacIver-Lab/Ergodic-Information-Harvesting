function out = LPF(data, Fs, varargin)
% Low Pass Filter
% Chen Chen
% Aug 13, 2014

%% Basic Parameters
PlotMode = 0;
if nargin == 2
    CutOffFreq = 1; % Hz
elseif nargin == 3
    CutOffFreq = varargin{1};
elseif nargin == 4
    CutOffFreq = varargin{1};
    switch varargin{2}
        case 'disp'
            PlotMode = 1;
    end
end
% Fs = 60; % Sample Freq in Hz

%% Calc Offset
data = double(data);
% MedData = medfilt1(data,20);
% MedData(1:2) = MedData(3:4);
Offset = mean(data);
DC_data = data - Offset;

%% Build LPF
% Low-Pass Filter
if length(Fs) > 384
    order = 128;
elseif length(Fs) > 192
    order = 64;
elseif length(Fs) > 64
    order = 32;
elseif length(Fs) > 48
    order = 16;
else
    order = 8; % fallback
end
[B,A] = fir1(order,CutOffFreq/(Fs*.5));

%% Filter Data
out = filtfilt(B,A,DC_data) + Offset;
% out = WaveletSmooth(data');

%% Plot
if PlotMode
    figure(1), clf, plot(data,'b'), hold on; grid on;
    plot(out,'r','LineWidth',2);
end
