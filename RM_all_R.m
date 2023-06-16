%
% --- Bob van Nifterik --- TU delft - ms3
%
%
% img: image NxN gray scale
% n: max order
% a,b,alph,bet: free parameters
% RM: matrix containing [0,max order] moments
%


function image = RM_all_R(RM,img, a, b, alph, bet, order)

[N, N] = size(img);
Rn =  RM_basis_R(N,order, alph,bet,a,b);
Rm = RM_basis_R(N,order, alph,bet,a,b);

for s = a:(b-1)
    for t = a:(b-1)
        sum2 = 0;
        for m = 0:order
            for n= 0:order
                sum2 = sum2 + RM(n+1,m+1)*Rn(n+1,s+1)*Rm(m+1,t+1);

            end
        end
        image(s+1,t+1) = sum2;
    end
end


end