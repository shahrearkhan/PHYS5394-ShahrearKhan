%% Function to generate 10000 trial values with uniform probability distributions
function y = customrandn(a,b)
    x = randn(1,10000);
    y = a*x + b;
end