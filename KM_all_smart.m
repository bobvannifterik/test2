%
%  Script calcualting Krauwchouk moments 
%  Implemented as discribed by https://link.springer.com/article/10.1007/s10851-017-0758-9
%
% script saves radial polynominals of same size to speed up the calcualtion
% of moments 
%
%  Bob van Nifterik -- TU Delft -- MS3 

function KM = KM_all_smart(img,p,order)


[N, N] = size(img);


%%

b = N;
% Kn =  KM_basis(p,order,N);
% Km = KM_basis(p,order,N);


% ----------------- save basis functions ---------------------
    if(~exist('Kn','var') && ~exist('Km','var'))
        flag =1;
        filecheckn =  strjoin(["Rn",int2str(order),"_",int2str(p),"_",int2str(N),".mat"],'');
        filecheckm =  strjoin(["Rm",int2str(order),"_",int2str(p),"_",int2str(N),".mat"],'');
        if (isfile(filecheckn) && isfile(filecheckm))
            load(filecheckn)
            load(filecheckm)
        else
            Kn =  KM_basis(p,order,N);
            Km = KM_basis(p,order,N);
            filenamen = strjoin(["Rn",int2str(order),"_",int2str(p),"_",int2str(N),".mat"],'');
            save(filenamen,'Kn')
            filenamem = strjoin(["Rm",int2str(order),"_",int2str(p),"_",int2str(N),".mat"],'');
            save(filenamem,'Km')
        end
    end
% -----------------------------------------------------------

for m = 0:order
    for n= 0:order % change for other dimensions
        sum1 = 0;
       
       
        for s = 0:(b-1)
            for t = 0:(b-1)
                
                sum1 = sum1 + Kn(n+1,s+1)*Km(m+1,t+1)*img(s+1, t+1);
            end
        end
        KM(n+1,m+1) = sum1;
        
    end
end


end





%%

function wbin = w(x,p,N)

wbin= binopdf(x,N,p);

end

function rhon = rho(n, p, N)
rhon = ((-1)^n)*((1-p)/p)^n * factorial(n)/pochhammer((-N),n);
end

function KmF = KM_basis(p,order,N)
for  x = 0:N-1
    K(0+1,x+1) = sqrt(w(x,p,N)/rho(0, p, N));
    K(0+2,x+1) = (1-(x/(N*p)))*sqrt(w(x,p,N)/rho(1, p, N));
    
    
    for n = 2:order

        ni = n;
        n = n-1;
        A = sqrt( (p*(N-n))  /((1-p)*(n+1)));
        B = sqrt(((p)^2 *(N-n)*(N-n+1))/((1-p)^2 *(n+1)*n ));
        K(ni+1,x+1) = (1/(p*(N-n)))*A*(N*p -2*n*p+n-x)*K(ni,x+1)- (1/(p*(N-n)))*B*n*(1-p)*K(ni-1,x+1);
    end
end
KmF = K;
end

%%


