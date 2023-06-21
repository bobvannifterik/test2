function [JFMB] = JFM_basis(order,rho,p,q)
% Implementation based on Fast computation of Jacobi-Fourier moments for 
% invariantimage recognition - Rahul Upneja, Chandan Singh 

n = order;
alpha = 1;

%---------------- recursive calculation of pn --------------

%- initialize for n=0, n=1
PN = (n >= 2) * ((gamma(p) / gamma(q)) * ones(size(rho))) + ...
     (n == 1) * ((gamma(p+1) / gamma(q)) * (1 - ((p+1)/q) * (rho))) + ...
     (n == 0) * ((gamma(p) / gamma(q)) * ones(size(rho)));

if n > 1
    PN1 = (gamma(p+1) / gamma(q)) * (1 - ((p+1)/q) * (rho));
    PN2 = (gamma(p) / gamma(q)) * ones(size(rho));
    for k = 2:n
        L1 = -((2*k+p-1)*(2*k+p-2))/(k*(q+k-1));
        L2 = (p+2*k-2) + (((k-1)*(q+k-2)*L1) / (p+2*k-3));
        L3 = ((p+2*k-4)*(p+2*k-3)/2) + ((q+k-3)*(k-2)*L1/2) - ((p+2*k-4)*L2);
        PN = (L1*(rho) + L2) .* PN1 + L3 .* PN2;
        PN2 = PN1;
        PN1 = PN;
    end
end

%---------------- recursive calculation of An --------------

%- initialize for n=0, n=1
AN = (n >= 1) * sqrt(gamma(q) / (gamma(p) * gamma(p-q+1))) + ...
     (n == 0) * sqrt(gamma(q) / (gamma(p) * gamma(p-q+1)));

if n > 0
    AN1 = sqrt(gamma(q) / (gamma(p) * gamma(p-q+1)));
    for k = 1:n
        AN = sqrt(k*(q+k-1)/((p+k-1)*(p-q+k))) * AN1;
        AN1 = AN;
    end
end



JFMB = sqrt((p + 2*n) * alpha * ((1 - rho.^alpha).^(p - q)) .* (rho.^(alpha*q - 1)) ./ rho) .* AN .* PN;

JFMB(rho == 0) = 0; % Set NaN values to 0

end
