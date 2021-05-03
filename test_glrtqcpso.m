%% Script to apply PSO to the data realization with Colored Gaussian Noise
% Signal parameters
nSamples = 512;
sampFreq = 512; % Hz

% Function handle to noise PSD
noisePSD = @(f) (f>=50 & f<=100).*(f-50).*(100-f)/625 + 1;
% Generate PSD vector for all positive DFT frequencies
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdVec = noisePSD(posFreq);

% Signal to Noise ratio of the true signal
snr = 10;
% True Signal Parameters
a1 = 10;
a2 = 3;
a3 = 3;
% Time samples
dataX = (0:(nSamples-1))/sampFreq;
% Squares and Cubes of the time values
dataXSq = dataX.^2;
dataXCb = dataX.^3;
% Search range of phase coefficients
rmin = [1, 1, 1];
rmax = [180, 10, 10];

%% Data realization
% Design FIR filter with T(f)= target PSD
fltrOrdr = 300;
% Generate Output Noise
outNoise = statgaussnoisegen(nSamples,[posFreq(:), psdVec(:)],fltrOrdr,sampFreq);
% Signal phase coefficients
qcCoefs = [a1, a2, a3];
% Generate Signal Vector
sigVec = crcbgenqcsig(dataX,snr,qcCoefs);
% Generate Data Realization
dataVec = outNoise + sigVec;
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
outStruct = glrtqcpso(inParams,struct('maxSteps',1000),nRuns);

%% Plots
figure;
hold on;
plot(dataX, dataVec, '.');
plot(dataX, sigVec);
for lpruns = 1:nRuns
    plot(dataX,outStruct.allRunsOutput(lpruns).estSig,'Color',[51,255,153]/255,'LineWidth',4.0);
end
plot(dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
legend('Data','Signal',...
       ['Estimated signal: ',num2str(nRuns),' runs'],...
       'Estimated signal: Best run');
disp(['Estimated parameters: a1=',num2str(outStruct.bestQcCoefs(1)),...
                             '; a2=',num2str(outStruct.bestQcCoefs(2)),...
                             '; a3=',num2str(outStruct.bestQcCoefs(3))]);
                         