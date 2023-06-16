
%
% --- Bob van Nifterik --- TU delft - ms3
% 
%
% F: image NxN gray scale
% n: max order
% moments: matrix containing [0,max order] moments 

function [img] = CHCM_1_R(F,n, moments)


[N,N] = size(F);
 nmax = n; 
 
for i = 1:N
    x= (i*2/N) -1;
    T(0+1,i) = 1;
    T(1+1,i) = x;
    for n = 2:nmax
        T(n+1,i) = (2*x*T(n-1+1,i)-T(n-2+1,i));
        
        %T2(n+1,i) = chebyshevT(n,x);
    end
end


x= [-1+(2/N): (2/N): 1];
for n = 0:nmax
    ni = n+1;
    
    for m = 0:nmax
        mi = m+1;
        
        sum = 0 ;
        for i = 1:N-1
            for j = 1:N-1
                
                sum = sum + T(ni,i)*T(mi,j)*F(i,j)*(1/sqrt(1-x(i)^2))*(1/sqrt(1-x(j)^2));
                
            end
        end
        
        if m == 0 && n ==0
            moments(ni,mi) = sum *(1/pi)*(1/pi)* (4/(N-1)^2);
        elseif m ~= n
            moments(ni,mi) = sum* (4/(N-1)^2);
        else
            
            moments(ni,mi) = sum* (4/(N-1)^2) *(2/pi)*(2/pi) ;
            
        end
    end
end


for i = 1:N
    for j = 1:N
        
        sum = 0 ;
        for n = 0:nmax
             ni = n+1;
            for m = 0:nmax
                mi = m+1;
                
                sum = sum + moments(ni,mi)*T(ni,i)*T(mi,j);
                Fr(i,j) = sum;
            end
        end
        %moments(ni,mi) = sum * (2*n+1)*(2*m+1)/(N-1)^2;
    end
    
end

img = abs(Fr);
