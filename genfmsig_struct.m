function sigVec = genfmsig_struct(dataX, snr, inParams)
% Generate a Freuqency Modulated (FM) sinusoid
% S = GENFMSIG(X, SNR, inParams)
% Generate a frequency modulated signal S. X is the vector
% of time stamps at which the samples of the signal are to be computed.
% SNR is the matched filtering signa-to-noise ratio.
% inParams is the structure that holds the values of signal parameters f0,
% f1 and b.

% Shahrear Khan Faisal, April 2021

% Set signal parameters in different fields of a structure
f0 = inParams.freq0;
f1 = inParams.freq1;
b = inParams.b0;
    
% Generate the Signal
phaseVec = f0*dataX + b*cos(2*pi*f1*dataX)/(2*pi);
sigVec = sin(2*pi*phaseVec);
sigVec = snr.*sigVec/norm(sigVec);
    