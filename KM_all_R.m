%
%  Script calcualting image reconsturction using Krauwchouk moments
% Bob van Nifterik - MS3 - TU Delft
% img: image NxN gray scale
% order: max order
% p: free parameter
% KM: matrix containing [0,max order] moments


function imgRec = KM_all_R(KM,img,p,order)
[N, N] = size(img);
Rn =  KM_basis(p,order,N);
Rm = KM_basis(p,order,N);
b =N;
for s = 0:(b-1)
    for t = 0:(b-1)
        sum2 = 0;
        for m = 0:order-1
            for n= 0:order-1
                sum2 = sum2 + KM(n+1,m+1)*Rn(n+1,s+1)*Rm(m+1,t+1);

            end
        end
        image(s+1,t+1) = sum2;
    end
end
imgRec = image;
end

function wbin = w(x,p,N)

wbin= binopdf(x,N,p);

end

function rhon = rho(n, p, N)
rhon = ((-1)^n)*((1-p)/p)^n * factorial(n)/pochhammer((-N),n);
end

function KMF = KM_basis(p,order,N)
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
KMF = K;
end