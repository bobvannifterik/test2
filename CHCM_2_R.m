
%
% --- Bob van Nifterik --- TU delft - ms3
% 
%  https://doi.org/10.7153/jca-03-02
%
% F: image NxN gray scale
% n: max order
% moments: matrix containing [0,max order] moments 


function [img] = CHCM_2_R(F,n, moments)
a = 1; % set a to 1 for CHCM first kind 

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
for i = 1:N
       
    for j = 1:N
        
        
        sum = 0 ;
        for n = 0:nmax
             ni = n+1;
            for m = 0:nmax
                mi = m+1;
                
                sum = sum + moments(ni,mi)*G(ni,i)*G(mi,j);
                Fr(i,j) = sum;
            end
        end
        %moments(ni,mi) = sum * (2*n+1)*(2*m+1)/(N-1)^2;
    end
    
end

img = (abs(Fr));

imshow(uint8((abs(Fr))));
