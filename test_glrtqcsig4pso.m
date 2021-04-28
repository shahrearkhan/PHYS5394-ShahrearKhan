%% Script to test the fitness function glrtqcsig4pso and compute the fitness
%% values

% Time samples
dataX = (0:(nSamples-1))/sampFreq;
% Squares and Cubes of the time values
dataXSq = dataX.^2;
dataXCb = dataX.^3;
% Quadratic Chirp signal parameters
% True Parameters
a1 = 10;
a2 = 5;
a3 = 2;
% Construct standardized coordinate values
a1_min = 9;
a1_max = 11;
a2_min = 4;
a2_max = 6;
a3_min = 1;
a3_max = 3;
% Create an array of values for a1
A = a1_min:0.001:a1_max;
% SNR
snr = 10;

X = zeros(length(A), 3);
for i = 1:length(A)
    X(i, 1) = (A(i) - a1_min)/(a1_max - a1_min);
    X(i, 2) = (a2 - a2_min)/(a2_max - a2_min);
    X(i, 3) = (a3 - a3_min)/(a3_max - a3_min);
end

%% Colored Gaussian Noise Generation
% Parameters
% Sampling frequency for noise realization
sampFreq = 2048; %Hz
% Number of samples to generate
nSamples = 4096;

% Produce PSD Vector
freqVec = 0:0.5:1024;
rng(10);
psdVec = rand(1, 2049);

%%
% Design FIR filter with T(f)= square root of target PSD
fltrOrdr = 500;
% Generate Output Noise
outNoise = statgaussnoisegen(nSamples,[freqVec(:), psdVec(:)],fltrOrdr,sampFreq);
% Define the coefficients
qcCoefs = [a1, a2, a3];
% Generate Signal Vector
sigVec = crcbgenqcsig(dataX,snr,qcCoefs);
% Generate Data Realization
dataVec = outNoise + sigVec;

%% Compute fitness value of GLRT
fitVal = zeros(1, length(A));
% Create a struct to pass necessary data
params = struct('dataX', dataX, 'dataXSq', dataXSq, 'dataXCb', dataXCb, 'dataVec', dataVec, 'psdVec', psdVec, 'snr', snr);
for i = 1:length(A)
    fitVal(i) = glrtqcsig4pso(X(i,:), params);
end
% Display Global minimum of the fitness function
disp('Global minimum of the fitness function: ');
disp(min(fitVal));
% Plot the Fitness values
plot(A, fitVal);
title('Fitness Value vs Signal Parameter');
xlabel('Signal Parameter values');
ylabel('Fitness Values');
