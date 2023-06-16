
%
% Bob van Nifterik - MS3 - TU Delft
%
% img: image NxN gray scale
% order: max order
% CHM: matrix containing [0,max order] moments 

function image = CHDM_all_R(img, CHM, order)
[N,N] = size(img);
tm = tcheb_rtn(order,N);
tn = tcheb_rtn(order,N);


for s = 0:N-1
    for t = 0:N-1
        sum2 = 0;
        for m = 0:order
            for n= 0:order
                sum2 = sum2 + CHM(n+1,m+1)*tn(n+1,s+1)*tm(m+1,t+1);
                
            end
        end
        image(s+1,t+1) = sum2;
    end
end

end 
function t = tcheb_rtn(order,N)

for x = 0:N-1
    t(1,x+1) = 1;
    t(2,x+1) = (2*x+1-N)/N;
    for n = 2:order
        t(n+1,x+1) =((2*n-1)*t(1+1,x+1)*(t(n,x+1))- (n-1)*(1-((n-1)^2/(N^2)))*t(n-1,x+1))/n;
    end
end

end