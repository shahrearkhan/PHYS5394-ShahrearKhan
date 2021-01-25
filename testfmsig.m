%% Plot the frequency modulated sinusoid
%Signal parameters
A = 10;
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

% Generate the signal
sigVec = genfmsig(timeVec,A,b,f0,f1);

% Plot the signal
figure;
plot(timeVec, sigVec, 'Marker', '.', 'Markersize', 24);
