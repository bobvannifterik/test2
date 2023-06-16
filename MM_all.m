% --- Bob van Nifterik --- TU delft - ms3
% 
%
% F: image NxN gray scale
% n: max order
% b,mu: free parameters
% moments: matrix containing [0,max order] moments 



function MM = MM_all(b,mu,img, order)

[N, N] = size(img);


for m = 0:order
   for n= 0:order % change for other dimensions  
        sum1 = 0;
    
        Hn =  MM_basis(b,mu,N,order);
        Hm = MM_basis(b,mu,N,order);
         
        for s = 0:(N-1) 
            for t = 0:(N-1) 
                
                sum1 = sum1 + Hn(n+1,s+1)*Hm(m+1,t+1)*img(s+1,t+1);
            end
        end
        MM(n+1,m+1) = sum1;
       
   end
end

end