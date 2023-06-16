
% Chebychef Moment calcualtion 
%
% Bob van Nifterik - MS3 - TU Delft
%
% img: image NxN gray scale
% order: max order
% moments: matrix containing [0,max order] moments 

function CHM = CHDM_all(img, order)

[N,N] = size(img);

tm = tcheb_rtn(order,N);
tn = tcheb_rtn(order,N);


for m = 0:order
    for n= 0:order % change for other dimensions
        sum1 = 0;
        for s = 0:N-1
            for t = 0:N-1
                sum1 = sum1 + tn(n+1,s+1)*tm(m+1,t+1)*img(s+1,t+1);
            end
        end
        CHM(n+1,m+1) = sum1/(rho_r(N,n)*rho_r(N,m));
        
    end
end
end


function rho = rho_r(N,n)

rho=1;
for i = 0:n
    %rho = (1/(N))*(N+i)*(i*2)*(N-i)/((2*i+1));
    rho = rho*(1/(N^2))*(N^2-i^2);
end
rho = (rho)/(2*n+1);
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