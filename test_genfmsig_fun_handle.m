%% Plot the frequency modulated sinusoid
%Signal parameters
f0 = 30;
f1 = 5;
b = 3;

% Maximum frequency of FM signal and sample interval
maxFreq = f0 + b*f1;
samplFreq = maxFreq*5;
samplIntrval = 1/samplFreq;

% Time samples
timeVec = 0:samplIntrval:1;

% Number of samples
nsamples = length(timeVec);

% Function input parameters in structure 
inParams = struct('freq0', f0, 'freq1', f1, 'b0', b);

% Create the function handle
H = @(snr) genfmsig_struct(timeVec, snr, inParams);

% SNR values
snr = [10, 12, 15];

% Plot the signal for different SNRs
for i = 1:length(snr)
    figure;
    plot(timeVec, H(snr(i)), 'Marker', '.', 'Markersize', 24);
    title(['Signal for SNR = ', num2str(snr(i))]);
    xlabel('Time (sec)');
    ylabel('Signal Amplitude');
end
