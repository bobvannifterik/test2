
%
% Bob van Nifterik - MS3 - TU Delft
%
% img: image NxN gray scale
% order: max order
% a1: free parameter 
% M: matrix containing [0,max order] moments 
% 

function image = CM_R(a1,M,img, order)

[N, N] = size(img);

Hn =  CM_basis(a1,N,order);
Hm = CM_basis(a1,N,order);

for s = 0:N-1
    for t = 0:N-1
        sum2 = 0;
        for m = 0:order-1
            for n= 0:order-1
                sum2 = sum2 + M(n+1,m+1)*Hn(n+1,s+1)*Hm(m+1,t+1);
                
            end
        end
        image(s+1,t+1) = sum2;
    end
end