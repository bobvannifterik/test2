% --- Bob van Nifterik --- TU delft - ms3
%
% F: image NxN gray scale
% n: max order
% b,mu: free parameters
% m: matrix containing [0,max order] moments

function image = MM_R(b,mu,img, order,M)
[N, N] = size(img);
Hn =  MM_basis(b,mu,N,order);
Hm = MM_basis(b,mu,N,order);

for s = 0:N-1
    for t = 0:N-1
        sum2 = 0;
        for m = 0:order
            for n= 0:order
                sum2 = sum2 + M(n+1,m+1)*Hn(n+1,s+1)*Hm(m+1,t+1);

            end
        end
        image(s+1,t+1) = sum2;
    end
end
