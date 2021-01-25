function sigVec = genfmsig(dataX,snr,b,f0,f1)
% Generate a Freuqency Modulated (FM) sinusoid
% S = GENFMSIG(X, SNR, B, F0, F1)
% Generate a frequency modulated signal S. X is the vector
% of time stamps at which the samples of the signal are to be computed.
% SNR is the matched filtering signa-to-noise ratio of S and B is the
% modulator index of the FM signal. F0 is the carrier frequency and F1 is
% the modulator frequency.

% Shahrear Khan Faisal, January 2021

    phaseVec = f0*dataX + b*cos(2*pi*f1*dataX)/(2*pi);
    sigVec = sin(2*pi*phaseVec);
    sigVec = snr*sigVec/norm(sigVec);