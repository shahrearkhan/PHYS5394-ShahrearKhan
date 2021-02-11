%% Filter Design
% Signal parameters
samplFreq = 1024;
nSamples = 2048;

% Time samples
timeVec = (0:(nSamples-1))/samplFreq;

% Generate the signals
sigVec1 = 10*sin(2*pi*100*timeVec);
sigVec2 = 5*sin(2*pi*200*timeVec + pi/6);
sigVec3 = 2.5*sin(2*pi*300*timeVec + pi/4);

% Add the signals
sigVec = sigVec1 + sigVec2 +  sigVec3;

% Length of data
dataLen = timeVec(end) - timeVec(1);

% DFT sample corresponding to Nyquist frequency
kNyq = floor(nSamples/2) + 1;

% Positive Fourier frequencies
posFreq = (0:(kNyq - 1))*(1/dataLen);

% FFT of the input signal
fftSig = fft(sigVec);
% Discard negative frequencies
fftSig = fftSig(1:kNyq);

% Design low-pass filter for sigVec1
filtOrdr = 30;
b1 = fir1(filtOrdr, 150/(samplFreq/2));
% Apply filter
filtSig1 = fftfilt(b1, sigVec);
% FFT of signal
fftSig1 = fft(filtSig1);
% Discard negative frequencies
fftSig1 = fftSig1(1:kNyq);

% Plot the signal from low pass filter
figure;
hold on;
subplot(2,3,1);
plot(timeVec, sigVec);
xlabel('time(sec)');
ylabel('signal amplitude');
title('Input Signal');
subplot(2,3,4);
plot(timeVec, sigVec);
xlabel('time(sec)');
ylabel('signal amplitude');
hold on;
plot(timeVec,filtSig1);
title('Signal from low-pass filter');

% Design band-pass filter for sigVec2
filtOrdr = 30;
b2 = fir1(filtOrdr, [190/(samplFreq/2) 250/(samplFreq/2)], 'bandpass');
% Apply filter
filtSig2 = fftfilt(b2, sigVec);
fftSig2 = fft(filtSig2);
% Discard negative frequencies
fftSig2 = fftSig2(1:kNyq);

% Plot the signal from band pass filter
subplot(2,3,2);
plot(timeVec, sigVec);
xlabel('time(sec)');
ylabel('signal amplitude');
title('Input Signal');
subplot(2,3,5);
plot(timeVec, sigVec);
hold on;
plot(timeVec, filtSig2);
xlabel('time(sec)');
ylabel('signal amplitude');
title('Signal from band-pass filter');

% Design high-pass filter for sigVec3
filtOrdr = 30;
b3 = fir1(filtOrdr, 0.59, 'high');
% Apply filter
filtSig3 = fftfilt(b3, sigVec);
fftSig3 = fft(filtSig3);
% Discard negative frequencies
fftSig3 = fftSig3(1:kNyq);

% Plot the signal from high pass filter
subplot(2,3,3);
plot(timeVec, sigVec);
xlabel('time(sec)');
ylabel('signal amplitude');
title('Input signal');
subplot(2,3,6);
plot(timeVec, sigVec);
hold on;
plot(timeVec, filtSig3);
xlabel('time(sec)');
ylabel('signal amplitude');
title('Signal from high-pass filter');
hold off;

figure;
hold on;
% Plot the periodograms
subplot(2,3,1);
plot(posFreq, abs(fftSig));
xlabel('Positive Frequency');
ylabel('Absolute value of fft signal');
title('Periodogram of input signal');
subplot(2,3,4);
plot(posFreq, abs(fftSig1)); 
xlabel('Positive Frequency');
ylabel('Absolute value of fft signal');
title('Periodogram of output signal from low-pass filter');

subplot(2,3,2);
plot(posFreq, abs(fftSig));
xlabel('Positive Frequency');
ylabel('Absolute value of fft signal');
title('Periodogram of input signal');
subplot(2,3,5);
plot(posFreq, abs(fftSig2));
xlabel('Positive Frequency');
ylabel('Absolute value of fft signal');
title('Periodogram of output signal from band-pass filter');

subplot(2,3,3);
plot(postFreq, abs(fftSig));
xlabel('Positive Frequency');
ylabel('Absolute value of fft signal');
title('Periodogram of input signal');
subplot(2,3,6);
plot(postFreq, abs(fftSig3));
xlabel('Positive Frequency');
ylabel('Absolute value of fft signal');
title('Periodogram of output signal from high-pass filter');
