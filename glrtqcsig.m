%% Function to compute GLRT of Quadratic Chirp(QC) signal
function glrt = glrtqcsig(dataVec, psdVals, a)

% Number of samples and sampling frequency.
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

    %% Compute GLRT
    % Quadratic Chirp signal
    sigVec = crcbgenqcsig(timeVec,10,a);
    % Template vector using the normsig4psd function
    [templateVec,~] = normsig4psd(sigVec,sampFreq,psdVals,1);
    % Calculate inner product of data with template
    llr = innerprodpsd(dataVec,templateVec,sampFreq,psdVals);
    %GLRT is its square
    glrt = llr^2;
    