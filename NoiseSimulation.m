%% Simulation of LIGO Noise
% Signal parameters
sampFreq = 1024;
nSamples = 16384;

% Load data from the file
data = load('iLIGOSensitivity.txt');

%% Insert missing frequencies and corresponding PSD values
y = zeros(length(data)+2, 2);
y(1, 1) = 0;
y(1, 2) = data(1, 1);
ind_1 = find(data(:,1) < 512);
y(2:length(ind_1)+1, 1) = data(1:length(ind_1), 1);
y(2:length(ind_1)+1, 2) = data(1:length(ind_1), 2);
y(length(ind_1)+2, 1) = 512;
y(length(ind_1)+2, 2) = interp1(data(:, 1), data(:, 2), 512);
y(length(ind_1)+3:end, 1) = data(length(ind_1)+1:end, 1);
y(length(ind_1)+3:end, 2) = data(length(ind_1)+1:end, 2);

% Target PSD
psdVals = y(:,2);
freqVec = y(:,1);
% Plot the target PSD
figure;
loglog(freqVec, psdVals);
title('Target PSD');
xlabel('Frequencies(Hz)');
ylabel('Target PSD');

%% Generate noise realization
inNoise = randn(1, nSamples);
% Design filter with Transfer function T(f) = sqrt(targetPSD)
% Filter parameters
fltrOrdr = 500;
freq = y(:,1)*512/max(y(:,1));
b = fir2(fltrOrdr, freq/(sampFreq/2), y(:,2));

% Filter the input noise
outNoise = sqrt(sampFreq)*fftfilt(b, inNoise);

% Estimate the PSD of output noise
[pxx, f] = pwelch(outNoise, 256, 128, [], sampFreq);
v = [f pxx];

% Specify psd values for certain frequencies
indx_1 = find(v(:,1) < 50);
indx_2 = find(v(:,1) > 700);
v(indx_1, 2) = v(v(:,1) == 48, 2);
v(indx_2, 2) = v(v(:,1) == 700, 2);
sqrt_v = sqrt(v(:,2));

% Plot the estimated psd
figure;
loglog(v(:, 1), sqrt_v);
title('Estimated PSD');
xlabel('Frequencies(Hz)');
ylabel('Estimated PSD values');
