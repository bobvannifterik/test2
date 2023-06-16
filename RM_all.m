
%
%  Script calcualting Racah moments
% --- Bob van Nifterik --- TU delft - ms3
% 
%
% img: image NxN gray scale
% n: max order
% a,b,alph,bet: free parameters
% RM: matrix containing [0,max order] moments 
%  


function RM= RM_all(img, a, b, alph, bet, order)

[N, N] = size(img);


for m = 0:order 
   for n= 0:order
        sum1 = 0;
        tic
        
        Rn =  RM_basis_R(N,order, alph,bet,a,b);
        Rm = RM_basis_R(N,order, alph,bet,a,b);
         t1 = toc;
        for s = a:(b-1) 
            for t = a:(b-1) 
                
                sum1 = sum1 + Rn(n+1,s+1)*Rm(m+1,t+1)*img(s+1-a,t+1-a);
            end
        end
        RM(n+1,m+1) = sum1;
       
   end
end


end