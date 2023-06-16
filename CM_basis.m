%
% Bob van Nifterik - MS3 - TU Delft
%
% img: image NxN gray scale
% order: max order
% a1: free parameter 
% moments: matrix containing [0,max order] moments 
% Method used: 
% https://doi.org/10.1007/s11042-018-6757-z

function c = CM_basis(a1,N,order)


for x = 0:N-1
    xi = x+1;
    %w = (exp(-a1)*a1^x)/factorial(x) not stable 
    
    wlog = log( exp(-a1)) + x*log(a1)-gammaln(x+1);
    w = exp(wlog);
    
    rho_0 = 1/1;
    c(1,xi) = sqrt(w/rho_0);
end

for x = 0:N-1
    xi = x+1;
    %w = (exp(-a1)*a1^x)/factorial(x); not stable 
        
    wlog = log( exp(-a1)) + x*log(a1)-gammaln(x+1);
    w = exp(wlog);
    
    rho_1 = 1/a1;
    c(2,xi) = ((a1-x)/a1)*sqrt(w/rho_1);
end

for x = 0:N-1 
    xi = 1+x;
    for n = 2:order
        ni = n+1;
            
        c(ni,xi) = ((a1-x+n-1)/a1)*sqrt(a1/(n))* c(ni-1,xi)- sqrt((n-1)/n)* c(ni-2,xi);
        
    end
end

end




