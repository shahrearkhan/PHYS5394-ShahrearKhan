%% Plotting the signals using two different sampling rates
%Signal parameters
A = 10;
f0 = 30;
f1 = 5;
b = 3;
% Instantaneous frequency after 3/20 sec is
maxFreq = f0 + b*f1;
% Sampling frequency using 5 times the maximum frequency
samplFreq_5 = 5*maxFreq;
samplIntrval_5 = 1/samplFreq_5;

% Time samples
timeVec_5 = 0:samplIntrval_5:1;

% Number of samples
nsamples_5 = length(timeVec_5);

% Generate the signal
sigVec1 = genfmsig(timeVec_5,A,b,f0,f1);

%Sampling frequency using half the maximum frequency
samplFreq_half = maxFreq/2;
samplIntrval_half = 1/samplFreq_half;

% Time samples
timeVec_half = 0:samplIntrval_half:1;

% Number of samples
nsamples = length(timeVec_half);

% Generate the signal
sigVec2 = genfmsig(timeVec_half,A,b,f0,f1);

% Plot both signals
figure;
subplot(2,1,1);
plot(timeVec_5, sigVec1, 'Marker', '.', 'Markersize', 24);
xlabel('time(sec)');
ylabel('signal value');
title(['Signal with sampling frequency ',num2str(samplFreq_5),' Hz']);
subplot(2,1,2);
plot(timeVec_half,sigVec2, 'Marker', '.', 'Markersize', 24);
xlabel('time(sec)');
ylabel('signal value');
title('Signal with sampling frequency half the maxmimum frequency');
