%% Reconstructing the signal by using Hermite functions.
function [aprx]=hermite_generate(x,co)
n=length(co);
H=hermite_system(x,n);
aprx=H*co;