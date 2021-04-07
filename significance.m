%%
addpath 'D:\Spring 2021\Statistical Methods\lab 11'

data1 = load('data1.txt');
data2 = load('data2.txt');
data3 = load('data3.txt');

a1 = 10;
a2 = 3;
a3 = 3;
a = [a1, a2, a3];

%% Parameters for data realization
% Number of samples and sampling frequency.
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

%% Supply PSD values
noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);

%% Determine GLRT values of 3 data realizations
glrt1 = glrtqcsig(data1', psdPosFreq, a);
glrt2 = glrtqcsig(data2', psdPosFreq, a);
glrt3 = glrtqcsig(data3', psdPosFreq, a);

%% Generate m data realizations
m = 90000;
gamma = zeros(1, m);
count_1 = 0;
count_2 = 0;
count_3 = 0;
for i = 1:m
    noise = randn(1, nSamples)*80;
    gamma(i) = glrtqcsig(noise, psdPosFreq, a);
end
for i = 1:m
    if gamma(i) >= glrt1
        count_1 = count_1 + 1;
    end
    if gamma(i) >= glrt2
        count_2 = count_2 + 1;
    end
    if gamma(i) >= glrt3
        count_3 = count_3 + 1;
    end
end

%% Calculate significance of the data realizations
signf_1 = count_1 / m; % significance of the first GLRT
disp(signf_1); 
signf_2 = count_2 / m; % significance of the first GLRT
disp(signf_2); 
signf_3 = count_3 / m; % significance of the first GLRT
disp(signf_3);
    