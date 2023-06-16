%%
% Recursive computation of racah poly
%
% Bob van Nifterik


function RMF = RM_basis_R(N,order, alph,bet,a,b)

rs = rhoRe(bet,alph,a,b);
dnsq = dnRe(bet,alph,a,b,N);

d0sq = dnsq(1);
d1sq = dnsq(2);
dndn1 = 1./dn0dn1(bet,alph,a,b,N);
dndn2 = 1./dn0dn2(bet,alph,a,b,N);

for i = 0:(N-1)
    
    s = a+i ;
    
    R(0+1,s+1) =  sqrt(rs(s+1)*(2*s+1)/d0sq); 
    R(1+1, s+1) =  (((-(a+s+1)*(s-a+bet+1)*(b+alph+s+1)*(b-s-1))/(2*s+1))+(((b+alph-s)*(a-bet+s)*(s-a)*(b+s))/(s*2+1)))*sqrt(rs(s+1)*(2*s+1)/d1sq);
    
end


for x = 0:N-1
    for n = 2:order
        s = x+a;
        
        An= n*(alph+bet+n)/ ((alph+bet+2*n-1)*(alph+bet+2*n));
        
        Bn = s*(s+1) - (a^2 +b^2 +(a-bet)^2+(b+alph)^2 -2)/4 + (alph+bet+2*n-2)*(alph+bet+2*n)/8 ...
            - ((bet^2-alph^2)*((b+alph/2)^2 - (a-(bet/2))^2))/(2*(alph+bet+n*2-2)*(alph+bet+2*n));
        
        Cn = -((alph +n -1)*(bet+n-1)/((alph +bet+2*n-2)*(alph+bet+2*n-1)))* ...
            ((a+b+((alph-bet)/2))^2 - (n-1+((alph+bet)/2))^2)*...
            ((b-a+((alph+bet)/2))^2 - (n-1+((alph+bet)/2))^2);
         
        
          R(n+1,s+1) = (Bn/An)*sqrt(dndn1(n+1))*R(n,s+1)+(Cn/An)*sqrt(dndn2(n+1))*R(n-1,s+1);
         
        %R(n+1,s+1) = (Bn/An)*sqrt(dnsq(n)/dnsq(n+1))*R(n,s+1)+(Cn/An)*sqrt(dnsq(n-1)/dnsq(n+1))*R(n-1,s+1);
        
    end
end

RMF = R;%.*(1./An);
end

function rs = rhoRe(bet,alph,a,b)
rho_a = gamma(alph+1)*gamma(alph+2*a+2)*gamma(bet+1)/((2*a+1)*gamma(2*a-bet+1));

for bi = (a+2):b
    rho_a = rho_a*(bi+alph-a-1)*(bi+alph+a)/((bi-a-1)*(bi+a));
end

rho_rec(a+1) = rho_a;
for s = (a+1):(b-1)
    rho_rec(s+1) =  rho_rec(s)*(a + s)*(s - a + bet)*(b + alph + s)*(b - s)/...
        ((a - bet + s)*(s - a)*(b + s)*(b + alph - s));
end
rs =rho_rec;
end


function dn = dnRe(bet,alph,a,b,n)

d02 = gamma(alph + 1)*gamma(bet + 1)*gamma(alph + bet + 2)*gamma(2*a + alph + 2)/...
    ((alph + bet + 1)*gamma(alph + bet + 1)*gamma(2*a + 1 - bet));

for bi = (a+2):(b)
    d02 = d02*(bi - a + alph + bet)*(a + bi + alph) / ((bi - a- 1)*(a + bi - bet - 1) );
end

dn2_rec(1) = d02;
for  ni = 1:n
    dn2_rec(ni+1) = dn2_rec(ni) *(alph + ni)*(bet + ni)*(b-a+alph+bet+ni)...
        *(a + b + alph + ni)*(b - a - ni)*(a + b - bet - ni) / (ni*(alph + bet + ni));
end
scale = (alph + bet + 1)./ (alph + bet + 1+ 2*[0:1:n]);

dn2_rec = scale.*dn2_rec;

dn = dn2_rec;
end




function d1d22 = dn0dn2(bet,alph,a,b,n)

d1d22 = zeros(1,n);
for ni = 2:n
    scale = (((alph + bet + 1)/ (alph + bet + 1+ 2*ni))/((alph + bet + 1)/ (alph + bet + 1+ 2*(ni-2))));
    nii = ni-1;
    d1d22(ni+1) = scale*((alph + ni)*(bet + ni)*(b-a+alph+bet+ni)...
        *(a + b + alph + ni)*(b - a - ni)*(a + b - bet - ni) / (ni*(alph + bet + ni)))*...
        (alph + nii)*(bet + nii)*(b-a+alph+bet+nii)...
        *(a + b + alph + nii)*(b - a - nii)*(a + b - bet - nii) / (nii*(alph + bet + nii));

end


end


function d1d2 = dn0dn1(bet,alph,a,b,n)

for ni = 1:n
    scale = ((alph + bet + 1)/ (alph + bet + 1+ 2*ni))/((alph + bet + 1)/ (alph + bet + 1+ 2*(ni-1)));
    d1d2(ni+1) = scale*(alph + ni)*(bet + ni)*(b-a+alph+bet+ni)...
        *(a + b + alph + ni)*(b - a - ni)*(a + b - bet - ni) /(ni*(alph + bet + ni)) ;

end

end

