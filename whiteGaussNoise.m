%% White Gaussian noise generation

% Load the data file
data = load('testData.txt');
% Signal parameters
% Sampling frequency for noise realization
sampFreq = 1024; %Hz
% Number of samples in the data file
nSamples = 16384;
% Time samples
timeVec = data(:,1);

%% Estimate PSD with signal free part of the noise
inSigFree = data(data(:,1) <= 5, 2);
[psdVals, freqVec] = pwelch(inSigFree, 256, 128, [], sampFreq);
% Plot the estimated PSD
figure;
plot(freqVec,psdVals);
title('Estimated PSD');
xlabel('Frequency (Hz)');
ylabel('PSD');

%% Generate White Gaussian noise realization
inNoise = data(:,2);
% Design FIR filter with T(f)= 1/square root of estimated PSD
% Estimate PSD with signal free part of the noise
sqrtPSD = 1./sqrt(psdVals);
fltrOrdr = 500;
b = fir2(fltrOrdr,freqVec/(sampFreq/2),sqrtPSD);
outNoise = sqrt(sampFreq)*fftfilt(b,inNoise);

% Plot the Colored Gaussian Noise Realization
figure;
plot(timeVec, inNoise);
title('Colored Noise');
xlabel('sampling time(sec)');
ylabel('data value');

% Plot the White Gaussian noise realization
figure;
plot(timeVec, outNoise);
title('White Noise');
xlabel('sampling time(sec)');
ylabel('data value');

%% Plot the spectrograms
% Spectrogram of colored noise
% Specify window length and overlap
winLen = 0.2;
ovrlp = 0.1;
% Convert to integer number of samples
winLenSampls = floor(winLen*sampFreq);
ovrlpSampls = floor(ovrlp*sampFreq);
% Generate the spectrogram
[S,F,T] = spectrogram(inNoise, winLenSampls, ovrlpSampls, [], sampFreq);
% Plot the spectrogram of Colored Noise
figure;
imagesc(T, F, abs(S)); 
axis xy;
title('Spectrogram of Colored Noise');
xlabel('Time(sec)');
ylabel('Frequency(Hz)');

%% Spectrogram of white noise
% Specify window length and overlap
winLen = 0.2;
ovrlp = 0.1;
% Convert to integer number of samples
winLenSampls = floor(winLen*sampFreq);
ovrlpSampls = floor(ovrlp*sampFreq);
% Generate the spectrogram
[S,F,T] = spectrogram(outNoise, winLenSampls, ovrlpSampls, [], sampFreq);
% Plot the spectrogram of White Gaussian Noise
figure;
imagesc(T, F, abs(S)); 
axis xy;
title('Spectrogram of White Noise');
xlabel('Time(sec)');
ylabel('Frequency(Hz)');
