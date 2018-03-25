%% Calculating the related Hermite coefficients 'co' of the signal 'signal' defined at 'x'.
function co = hermite_coeff(signal,x,n)
N=length(signal);
s=reshape(signal,N,1);
H=hermite_system(x,n);
h=x(2)-x(1); %stepsize of the uniform discretization
co=H'*s*h;

