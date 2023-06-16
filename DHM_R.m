%
% Dual Hahn implementation -- s(1+s)
% 
% Bob van Nifterik - MS3 - TU Delft
% img: image NxN gray scale
% order: max order
% a,b,c: free parameter 
% DHMF: matrix containing [0,max order] moments 


function image = DHM_R(img, DHM, order,a,b,c)

                Hn =  DHM_basis_alt(img,order,a,b,c);
                Hm = DHM_basis_alt(img,order,a,b,c);
                
for s = a:(b-1)
    for t = a:(b-1)
        sum2 = 0;
        for m = 0:order
            for n= 0:order
                sum2 = sum2 + DHM(n+1,m+1)*Hn(n+1,s+1)*Hm(m+1,t+1);
                
            end
        end
        image(s+1-a,t+1-a)= sum2;%image(s+1,t+1) = sum2;
    end
end



end