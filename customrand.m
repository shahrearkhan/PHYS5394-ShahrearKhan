%% Function to generate 10000 trial values with uniform probability distributions
function y = customrand(a,b)
    x = rand(1,10000);
    y = x*(b-a) + a;
end