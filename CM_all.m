%
% Bob van Nifterik - MS3 - TU Delft
%
% img: image NxN gray scale
% order: max order
% a1: free parameter 
% moments: matrix containing [0,max order] moments 

function CM = CM_all(a1,img, order)

[N, N] = size(img);


for m = 0:order 
   for n= 0:order % change for other dimensions  
        sum1 = 0;
    
        Hn =  CM_basis(a1,N,order);
        Hm = CM_basis(a1,N,order);
         
        for s = 0:(N-1) 
            for t = 0:(N-1) 
                
                sum1 = sum1 + Hn(n+1,s+1)*Hm(m+1,t+1)*img(s+1,t+1);
            end
        end
        CM(n+1,m+1) = sum1;
       
   end
end
end