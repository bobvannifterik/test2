function HMF = HM_basis_alt(N,order, alpha,beta,c)
% clear
% N =50;
% order = 4;
% a = 10;
% c = 10;
% b = 50;
a = (alpha+beta)/2;
b = a+N;
c = (beta-alpha)/2;

E = gammaln(2*a+1)+gammaln(b-c);
F = gammaln(b+a+1)+gammaln(a-c+1);
Hhat0a = sqrt(2*a+1) * exp((E-F)/2);
Hhat0s = zeros(1,(N+a));
Hhat0s(a+1) = Hhat0a;
for i = 1:(N-1) 
    s = i+a;
    Hhat0s(s+1) = sqrt(((a+s)*(c+s)*(b-s)*(2*s+1))/((s-a)*(b+s)*(s-c)*(2*s-1)))*Hhat0s(s-1+1);
end

A = sqrt(1/((a+c+1)*(b-c-1)*(b-a-1)));
Hhat1s = 0;
for i = 0:N-1
     s = i+a;
     Hhat1s(s+1) = Hhat0s(s+1)*A*((s-a)*(s+b)*(s-c)+(s+a+1)*(s+c+1)*(s-b+1))/(2*s+1) ;

end

Hhatns = zeros(order,(N+a));
Hhatns(1,:) = Hhat0s;
Hhatns(2,:) = Hhat1s;


for i = 0:N-1
     for n = 2:order
         s = i+a;
         A = (1/ n)*( s*(s + 1) - a*b + a*c - b*c - (b - a - c - 1)*(2*n - 1) + 2*(n - 1)^2 );
         B = ( n / (( a + c + n )*( b-a-n )*( b-c-n )));
         C = (-1/n)*(a + c +n - 1)*( b-a-n +1)*( b-c-n+ 1);
         D = ( n*( n-1 )/(( a + c + n )*( a + c + n- 1)*(b-a-n + 1)*(b-a -n)*( b-c -n + 1)*( b-c-n ) ) );
         
         Hhatns(n+1, s+1) = A*sqrt(B)* Hhatns(n-1+1,s+1) + C*sqrt(D)* Hhatns(n-2+1,s+1);
         
     end
     
end
 
DHMF = Hhatns;
end