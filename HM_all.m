
%
% Bob van Nifterik - MS3 - TU Delft
% img: image NxN gray scale
% order: max order
% a,b: free parameter
% HM: matrix containing [0,max order] moments 

function HM = HM_all(a,b,img, order)

[N, N] = size(img);


for m = 0:order 
   for n= 0:order 
        sum1 = 0;
        tic
        Hn =  HM_basis(a,b,N, order);
        Hm = HM_basis(a,b,N, order);
        
        for s = 0:(N-1) 
            for t = 0:(N-1) 
                
                sum1 = sum1 + Hn(n+1,s+1)*Hm(m+1,t+1)*img(s+1,t+1);
            end
        end
        HM(n+1,m+1) = sum1;
       
   end
end

end