%
%
% Bob van Nifterik - MS3 - TU Delft
% img: image NxN gray scale
% order: max order
% DHMF: matrix containing [0,max order] moments 
%

function image = GHM_R(img,GHM,order)
[N,N] = size(img);

Hn =    GHM_optimal_basis(N,order);
Hm =  GHM_optimal_basis(N,order);

for s = 0:N-1
    for t = 0:N-1
        sum2 = 0;
        for m = 0:order
            for n= 0:order
                sum2 = sum2 + GHM(n+1,m+1)*Hn(n+1,s+1)*Hm(m+1,t+1);

            end
        end
        image(s+1,t+1) = sum2;
    end
end
end