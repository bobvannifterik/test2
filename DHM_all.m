%
% Dual Hahn implementation -- s(1+s)
% 
% Bob van Nifterik - MS3 - TU Delft
%
% img: image NxN gray scale
% order: max order
% a,b,c: free parameter 
% DHM: matrix containing [0,max order] moments 
% 


function DHM = DHM_all(img, order,a,b,c)



for m = 0:order 
   for n= 0:order 
        sum1 = 0;
        tic
        Hn =  DHM_basis_alt(img,n,a,b,c);
        Hm = DHM_basis_alt(img,m,a,b,c);
         t1 = toc;
        for s = a:(b-1)
            for t = a:(b-1) 
                
                sum1 = sum1 + Hn(n+1,s+1)*Hm(m+1,t+1)*img(s+1-a,t+1-a);
            end
        end
        DHM(n+1,m+1) = sum1;
       
   end
end


end