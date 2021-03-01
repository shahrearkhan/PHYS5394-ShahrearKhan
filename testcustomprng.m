%% Histograms of uniform and normal probability distribution functions

% Generate uniform pdf
u = customrand(-2, 1, 10000);
% Generate theoretical uniform pdf
a = 1/(1-(-2));
u_th = ones(1,10000)*a;

% Generate normal pdf
n = customrandn(1.5,2.0,10000);
% Generate theoretical normal pdf
x = -4:0.2:8; % create bins
n_th = normpdf(x, 2.0, 1.5);

% Plot the coresponding histograms
% Histogram of uniform pdf
figure;
histogram(u_th, 'normalization', 'pdf', 'BinLimits', [-2, 1]);
hold on;
histogram(u, 'normalization', 'pdf');
title('Histogram of trial values from a uniform pdf');
xlabel('Trial values');
ylabel('Number of trial values in each bin');
legend('Theoretical pdf','Real pdf');
hold off;

% Histogram of normal pdf
figure;
plot(x, n_th);                      
hold on;
histogram(n, 'normalization', 'pdf');
title('Histogram of trial values from a normal pdf');
xlabel('Trial values');
ylabel('Number of trial values in each bin');
legend('Theoretical pdf','Real pdf');
hold off;
