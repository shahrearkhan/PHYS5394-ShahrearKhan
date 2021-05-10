% addpath 'D:\Spring 2021\Statistical Methods\Lab 7'
% addpath 'D:\Spring 2021\Statistical Methods\lab 11'
% addpath 'D:\Spring 2021\Statistical Methods\Final 4'

%% Script to apply PSO to the data realization
% Load the noise and data realizations
noise1 = load('TrainingData.mat');
data1 = noise1.trainData;
noise2 = load('analysisData.mat');
data2 = noise2.dataVec;
% Signal parameters
nSamples = 2048;
sampFreq = 1024; %Hz

% Estimate noise PSD
[pxx, f] = pwelch(data1, 2048, 1024, [], sampFreq);
psdVec = pxx';
freqVec = f';
% Plot PSD
plot(freqVec, psdVec);
title('Power Spectral Density vs Frequency');
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density');

% Signal to Noise ratio
snr = 10;
% Time samples
dataX = (0:(nSamples-1))/sampFreq;
% Squares and Cubes of the time values
dataXSq = dataX.^2;
dataXCb = dataX.^3;
% Search range of phase coefficients
rmin = [40, 1, 1];
rmax = [100, 50, 15];

% Generate Data Realization
dataVec = data2;
% Number of parallel PSO runs
nRuns = 8;
% Create the input parameter struct
inParams = struct('dataX', dataX,...
                  'dataXSq', dataXSq,...
                  'dataXCb', dataXCb,...
                  'dataVec', dataVec,...
                  'psdVec', psdVec,...
                  'snr', snr,...
                  'sampFreq', sampFreq,...
                  'rmin', rmin,...
                  'rmax', rmax);
              
% Outputs from glrtqcpso
outStruct = glrtqcpso(inParams,struct('maxSteps',2000),nRuns);

%% Plots
figure;
hold on;
plot(dataX, dataVec, '.');
for lpruns = 1:nRuns
    plot(dataX,outStruct.allRunsOutput(lpruns).estSig,'Color',[51,255,153]/255,'LineWidth',4.0);
end
plot(dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
title('Data and Signals');
xlabel('Time(sec)');
ylabel('Amplitude');
legend('Data','Signal',...
       ['Estimated signal: ',num2str(nRuns),' runs'],...
       'Estimated signal: Best run');
disp(['Estimated parameters: a1=',num2str(outStruct.bestQcCoefs(1)),...
                             '; a2=',num2str(outStruct.bestQcCoefs(2)),...
                             '; a3=',num2str(outStruct.bestQcCoefs(3))]);

% Estimated quadratic chirp coefficients
qcCoefs = [outStruct.bestQcCoefs(1), outStruct.bestQcCoefs(2), outStruct.bestQcCoefs(3)];
                         
% GLRT value of the data realization
%SDM**************************
%glrt = outStruct.allRunsOutput.fitVal;
%Use the fitness value from the best PSO run
glrt = outStruct.bestFitness;
%Use the negative of this fitness value (PSO used the negative for
%minimization)
glrt = -glrt;
%******************************
% Likelihood Ratio
llr = sqrt(glrt);

%% Generate m colored noise data realizations
m = 1000;
gamma = zeros(1, m);
count = 0;
% Filter Order
fltrOrdr = 100;
for i = 1:m
    noise = statgaussnoisegen(nSamples,[freqVec(:),psdVec(:)],fltrOrdr,sampFreq);
    gamma(i) = glrtqcsig(noise, psdVec, qcCoefs);
    if gamma(i) >= glrt
        count = count + 1;
    end
end

%% Calculate significance
signf = count / m; % significance of the GLRT
disp(['Significance = ', num2str(signf)]);
disp('No detection');
