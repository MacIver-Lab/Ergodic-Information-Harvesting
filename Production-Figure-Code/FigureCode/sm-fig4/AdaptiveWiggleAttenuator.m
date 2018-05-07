function [sPosOut, AttenuateMetrics, distGain] = AdaptiveWiggleAttenuator(sPos, dt, freqWin, stopBandAtt, plotSW, savePDF)
%% Filter Sensor Trajectory within the Input Frequency Window
% 
% Chen Chen
% 11/22/2017

Fs = 1.0 / dt;  % Sampling Frequency
if isempty(savePDF)
    savePDF = 0;
end

%% Configure Filter - Method #1 - FIR LPF
% Fpass = freqWin(1);        % Passband Frequency
% Fstop = freqWin(2);        % Stopband Frequency
% rp    = 0.5;              % Passband Ripple
% rs    = 200;              % Stopband Ripple
% fir_a = [1, 10^(-attFactor/20)];   % Desired amplitude on the pass and stop band
% 
% dev = [(10^(rp/20)-1)/(10^(rp/20)+1)  10^(-rs/20)]; 
% 
% % Calculate the order from the parameters using FIRPMORD
% % Parks-McClellan optimal FIR filter order estimation
% [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), fir_a, dev);
% 
% % Calculate the coefficients using the FIRPM function.
% b  = firpm(N, Fo, Ao, W);
% Hd = dsp.FIRFilter('Numerator', b);
% 
% % Convert to Transfer Function
% [b, a] = tf(Hd);

%% Configure Filter - Methods #2 - IIR LPF
lpFilt = designfilt('lowpassiir',...
         'PassbandFrequency', freqWin(1),...
         'StopbandFrequency', freqWin(2), ...
         'StopbandAttenuation', stopBandAtt, ...
         'PassbandRipple', 5e-1, ...
         'DesignMethod', 'butter', ...
         'SampleRate', Fs);

%% FIR Filter and FFT
sPosOut = filtfilt(lpFilt, sPos);

% Parameters
FreqResolution = 0.001; % [Hz]
% MaxFrequency = 2.0; % [Hz]
PlotMaxFrequency = 1.5; % [Hz]
sPosFFT = fft(sPos-mean(sPos), round(Fs/FreqResolution));
sPosFFT_filt = fft(sPosOut-mean(sPosOut), round(Fs/FreqResolution));
distGain = sum(abs(diff(sPosOut))) / sum(abs(diff(sPos)));

% Calculate RMS Gain for the Filter
% rmsGain = rms( gainFactor );
% fprintf('Dist Gain = %f\n',distGain);

%% Evaluate Filter's Magnitude Response
nPts = 2^12;
% freqIdx = 0:Fs/nPts:2;
% [cutoffIdx, ~] = findNearestIdx(freqIdx, freqWin(1));

% h = freqz(b, a, nPts, 'whole', Fs);
% h = freqz(lpFilt, nPts, 'whole', Fs);
% magResponsedB = 20.0 * log10(abs(h(1:length(freqIdx))));

freqWin(2) = 1.5;
[h, f] = freqz(lpFilt, nPts, Fs);
[cutoffIdx, ~] = findNearestIdx(f, freqWin(1));
freqIdx = f(f < freqWin(2));
magResponsedB = 20.0 * log10(abs(h(f < freqWin(2))));


%% Assemble Output Data Structure
AttenuateMetrics.Fs = Fs;
AttenuateMetrics.attenGain = stopBandAtt;
AttenuateMetrics.freqWin = freqWin;
AttenuateMetrics.magResponsedB = magResponsedB;
AttenuateMetrics.magResponseFreqIdx = freqIdx;

%% Plot
if plotSW
    figure(1); clf;
    % First segment - pass band
    plot(freqIdx(1:cutoffIdx), magResponsedB(1:cutoffIdx), 'LineWidth', 2); hold on;
    hLine = line([freqWin(1), freqWin(1)], [-200, 0], ...
        'LineStyle', '--', 'LineWidth', 2);
    
    % Second segment - attenuating band
    plot(freqIdx(cutoffIdx:end), magResponsedB(cutoffIdx:end), 'LineWidth', 2)
    
    opt = [];
    opt.BoxDim = [8,5];
    opt.ShowBox = 'off';
    opt.XMinorTick = 'off';
    opt.YMinorTick = 'off';
    opt.XLabel = 'Frequency (Hz)';
    opt.YLabel = 'Magnitude (dB)';
    opt.XLim = [0, freqWin(2)];
    opt.XTick = [0, 0.2, 0.5:0.5:2];
    opt.YLim = [-200, 2];
    opt.YTick = -200:40:0;
    opt.FontName = 'Helvetica';
    setPlotProp(opt);
    hLine.LineStyle = '--';
    if savePDF
        print(sprintf('Filter-%ddB-MagResponse.pdf', stopBandAtt),'-dpdf');
    end
    
    figure(2);clf;
    plot(FreqResolution:FreqResolution:PlotMaxFrequency, abs(sPosFFT(2:PlotMaxFrequency/FreqResolution+1)),...
        'LineWidth',2); hold on;
    plot(FreqResolution:FreqResolution:PlotMaxFrequency, abs(sPosFFT_filt(2:PlotMaxFrequency/FreqResolution+1)),...
        'LineWidth',2);
    hLine = line([freqWin(1), freqWin(1)], [0, 104], ...
        'LineStyle', '--', 'LineWidth', 2);
    legend('Original FFT', sprintf('Filtered (StopBandAttenuation = %ddB)',stopBandAtt));
    xlabel('Frequency (Hz)');
    set(gca,'YTick',[]);
    set(gca,'XTick',[0, 0.2, 0.5:0.5:1.5]);
    ylabel('Normalized Gain');
    opt = [];
    opt.BoxDim = [8,5];
    opt.ShowBox = 'off';
    opt.XMinorTick = 'off';
    opt.YMinorTick = 'off';
    opt.XLabel = 'Frequency (Hz)';
    opt.YLabel = 'Magnitude';
    opt.YLim = [0, 105];
    opt.XLim = [0, freqWin(2)];
    opt.FontName = 'Helvetica';
    setPlotProp(opt);
    hLine.LineStyle = '--';
    if savePDF
        print(sprintf('Filter-%ddB-FFT.pdf', stopBandAtt),'-dpdf');
    end
    
    clf;
    simTimestamp = (1:length(sPos)) * dt;
    plot(simTimestamp, sPos, 'LineWidth',2); hold on;
    plot(simTimestamp, sPosOut, 'LineWidth',2);
    legend('Original Trajectory', sprintf('Filtered (StopBandAttenuation = %ddB)',stopBandAtt));
    xlabel('Time');
    ylabel('Position');
    set(gca,'XLim',[0,dt*length(sPos)+0.01]);
    set(gca,'YLim',[0.1,0.9]);
    opt = [];
    opt.BoxDim = [8,5];
    opt.ShowBox = 'off';
    opt.XMinorTick = 'off';
    opt.YMinorTick = 'off';
    opt.XLabel = 'Time (Seconds)';
    opt.YLabel = 'Position';
    opt.YTick = [];
    opt.FontName = 'Helvetica';
    setPlotProp(opt);
    if savePDF
        print(sprintf('Filter-%ddB-Traj.pdf', stopBandAtt),'-dpdf');
    end
end

function [idx, val] = findNearestIdx(x, val)
[val, idx] = min(abs(x-val));
