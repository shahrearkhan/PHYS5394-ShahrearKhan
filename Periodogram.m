%% Plot the frequency modulated sinusoid and it's periodogram
%Signal parameters
A = 10;
f0 = 100;
f1 = 2;
b = 1;

% Instantaneous frequcncy after 3/20 sec is
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
plot(timeVec,sigVec, 'Marker', '.', 'Markersize', 24);

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
