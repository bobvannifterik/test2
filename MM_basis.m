%    meixner basis function 
%    Implemented as  10.1049/iet-ipr.2009.0195
%    recusrive over N With modificaitions for calculation of w recursilavly  
%
%    Bob van Nifterik
%    
   
function h = MM_basis(b,mu,N,order)

for x = 0:N-1
    xi = x+1;
    
    w = wX(x,b,mu);
    
    ro_0 = 1*(1)/((mu^0)*(1-mu)^b);
    
    h(1,xi) = sqrt(w/ro_0); 
end

for x = 0:N-1
   xi = x+1;
    
    w = wX(x,b,mu);
    ro_1 = 1*(b)/((mu^1)*(1-mu)^b);
    
    h(2,xi) = (b+x-(x/mu))*sqrt(w/ro_1); 
end

for x = 0:N-1 
    xi = 1+x;
    for n = 2:order
        ni = n+1;
        
        A = mu/(mu-1);
        B = (x*(1-mu)-mu*(n+b-1)-n+1)/(1-mu);
        C = ((n-1)*(n-2+b))/(mu-1);
        D = sqrt(mu/(n*(b+n-1)));
        E = sqrt((mu^2)/(n*(n-1)*(b+n-2)*(b+n-1)));
       
        
        h(ni,xi) = B*D*h(ni-1,xi)/A - C*E*h(ni-2,xi)/A;
        
        
    end
end
end

function w = wX(xmax,b,mu)
    w=factorial(b-1);
    for x = 1:xmax
        w =  ((x +b-1)/(x))* w;
    end
    w = ((mu^xmax)/gamma(b))*w;
end