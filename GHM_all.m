%
% Bob van Nifterik - MS3 - TU Delft
% img: image NxN gray scale
% order: max order
% Moments: matrix containing [0,max order] moments 


function [moments] = GHM_all(img,order)
 
[N,N] = size(img);

       
Hn =  GHM_optimal_basis(N,order);
Hm = GHM_optimal_basis(N,order);

for m = 0:order 
   for n= 0:order % change for other dimensions  
        sum1 = 0;
        
         
        for s = 0:N-1
            for t = 0:(N-1) 
                
                sum1 = sum1 + Hn(n+1,s+1)*Hm(m+1,t+1)*img(s+1,t+1);
            end
        end
        GHM(n+1,m+1) = sum1;
       
   end
end

moments = GHM;