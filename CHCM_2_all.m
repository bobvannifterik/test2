
%
% --- Bob van Nifterik --- TU delft - ms3
% 
%  https://doi.org/10.7153/jca-03-02
%
% F: image NxN gray scale
% n: max order
% moments: matrix containing [0,max order] moments 



function [moments] = CHCM_2_all(F,n)

a = 1; % lazy method 
%
% F: image NxN gray scale
% n: max order
% moments: matrix containing [0,max order] moments 
%
 
[N,N] = size(F);
 nmax = n; 
 


for i = 1:N
    x= (i*2/N) -1;
    G(0+1,i) = 1;
    G(1+1,i) = 2*a*x;
    for n = 2:nmax
        G(n+1,i) = (2*(n+a-1)*x*G(n-1+1,i)-(n+2*a-2)*G(n-2+1,i))/n;
        %G(n+1,i) = gegenbauerC(n, ((i/N)*2-1), a);
    end
end

x= [-1+(2/N): (2/N): 1];
for n = 0:nmax
    ni = n+1;
    
    for m = 0:nmax
        mi = m+1;
        
        sum = 0 ;
        for i = 1:N
            for j = 1:N
                sum = sum + G(ni,i)*G(mi,j)*F(i,j)*((1-x(j)^2)^(a-0.5))*((1-x(i)^2)^(a-0.5));
            end
        end
        moments(ni,mi) = sum * (4/(N-1)^2)*((factorial(n)*(n+a)*abs(gamma(a))^2)/(pi*gamma(n+2*a)*(2^(1-2*a))))* (factorial(m)*(m+a)*abs(gamma(a))^2)/(pi*gamma(m+2*a)*(2^(1-2*a)))    ;
    end
    
end








