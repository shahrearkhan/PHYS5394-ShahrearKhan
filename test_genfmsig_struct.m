%% Plot the frequency modulated sinusoid
%Signal parameters
f0 = 30;
f1 = 5;
b = 3;
snr = 10;

% Maximum frequency of FM signal and sample interval
maxFreq = f0 + b*f1;
samplFreq = maxFreq*5;
samplIntrval = 1/samplFreq;

% Time samples
timeVec = 0:samplIntrval:1;

% Number of samples
nsamples = length(timeVec);

inParams = struct('freq0', f0, 'freq1', f1, 'b0', b);
sigVec = genfmsig_struct(timeVec, snr, inParams);
% Plot the Signal time series
plot(timeVec, sigVec, 'Marker', '.', 'Markersize', 24);
title('Signal time series');
xlabel('Time (sec)');
ylabel('Signal Amplitude');

% Plot the periodogram
% Length of data
dataLen = timeVec(end) - timeVec(1);

% DFT sample corresponding to Nyquist frequency
kNyq = floor(nsamples/2) + 1;

% Positive Fourier frequencies
posFreq = (0:(kNyq - 1))*(1/dataLen);

% FFT of signal
fftSig = fft(sigVec);

% Discard negative frequencies
fftSig = fftSig(1:kNyq);

% Plot periodogram
figure;
plot(posFreq, abs(fftSig)); 
title('Periodogram of the Signal');
xlabel('Positive Frequencies (Hz)');
ylabel('Magnitude');

