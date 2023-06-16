% --- Bob van Nifterik --- TU delft - ms3
% 
%
% F: image NxN gray scale
% n: max order
% moments: matrix containing [0,max order] moments 

function [img] = LM_R(F,n,moments)
 
 [N,N] = size(F);
 nmax = n; 
 


for i = 1:N
    x= (i*2/N) -1;
    P(0+1,i) = 1;
    P(1+1,i) = x;
    for n = 2:nmax
        P(n+1,i) = ((2*n-1)*x*P(n-1+1,i)-(n-1)*P(n-2+1,i))/n;
        %PP(n+1,i) = legendreP(n, ((i/N)*2-1));
    end
end


for i = 1:N
   
    
    for j = 1:N
        
        
        sum = 0 ;
        for n = 0:nmax
             ni = n+1;
            for m = 0:nmax
                mi = m+1;
                
                sum = sum + moments(ni,mi)*P(ni,i)*P(mi,j);
                Fr(i,j) = sum;
            end
        end
        %moments(ni,mi) = sum * (2*n+1)*(2*m+1)/(N-1)^2;
    end
    
end


img = (abs(Fr));


%imshow(uint8((abs(Fr))));