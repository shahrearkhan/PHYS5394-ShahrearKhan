%% Function to generate trial values with uniform probability distributions
function y = customrandn(a,b,c)
    x = randn(1, c);
    y = a*x + b;
end
