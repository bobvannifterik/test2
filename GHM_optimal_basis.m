%
%  Bob van Nifterik -
%
% - Normalized Gaussian–Hermite basis function  -
% implemented as described by  https://doi.org/10.1016/j.sigpro.2011.04.012
%
% This implementation optemizes the varince parameter sigma over the orders


function [h_ha] = GHM_optimal_basis(Nimg,order)
N_sample = (1/2)*Nimg;
sig = 0.9*order^-0.52;
h(1,:) = ones(1,length(-1:(1/N_sample):1));
for x = -1:(1/N_sample):1
    xi = round((x+1)*N_sample+1);

    h_ha(1,xi) = (h(1,xi))*exp((-x^2)/(2*sig^2))/(sig*sqrt(pi)*factorial(0)*2^0)^0.5;
end

for x = -1:(1/N_sample):1
    xi = round((x+1)*N_sample+1);

    h(2,xi) = 2*(x/sig);
    h_ha(2,xi) = (h(2,xi))*exp((-x^2)/(2*sig^2))/(sig*sqrt(pi)*factorial(1)*2^1)^0.5;
    
end

for n = 2:order
    ni = n+1;

    for x = -1:(1/N_sample):1
        xi = round((x+1)*N_sample+1);
        
        h(ni,xi) =   2*(x/sig)*h(ni-1,xi) - 2*(n-1)* h(ni-2,xi);
        h_ha(ni,xi) = (h(ni,xi))*exp((-x^2)/(2*sig^2))/(sig*sqrt(pi)*factorial(n)*2^n)^0.5;
        
    end
    
end


end

