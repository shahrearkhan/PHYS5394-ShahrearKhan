%% Function to compute GLRT of Quadratic Chirp(QC) signal
function [fitVal, varargout] = glrtqcsig4pso(xVec, params)

%% Compute GLRT
%rows: points
%columns: coordinates of a point
[nVecs,~]=size(xVec);

%storage for fitness values
fitVal = zeros(nVecs,1);

%Check for out of bound coordinates and flag them
validPts = crcbchkstdsrchrng(xVec);
%Set fitness for invalid points to infty
fitVal(~validPts)=inf;

for lpc = 1:nVecs
    if validPts(lpc)
        x = xVec(lpc,:);
        % Compute fitness value using ssrqc function
        fitVal(lpc) = ssrqc(x, params);
    end
end

%Return real coordinates if requested
if nargout > 1
    varargout{1}=xVec;
end

%% Sum of squared residuals after maximizing over amplitude parameter
function ssrVal = ssrqc(x,params)
%Generate normalized quadratic chirp
phaseVec = x(1)*params.dataX + x(2)*params.dataXSq + x(3)*params.dataXCb;
sigVec = sin(2*pi*phaseVec);
sampFreq = 2048;
% Template vector using the normsig4psd function
[templateVec,~] = normsig4psd(sigVec,sampFreq,params.psdVec,params.snr);
% Calculate inner product of data with template
llr = innerprodpsd(params.dataVec,templateVec,sampFreq,params.psdVec);
%GLRT is its square
ssrVal = -llr^2;
