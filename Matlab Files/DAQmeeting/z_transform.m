clc
clear all
load('roosload.mat');

p1=-150;
p2=-150;
p3=-150;
p4=-150;
a=(- p1 - p2 - p3 - p4);
b=(p1*p2 + p1*p3 + p1*p4 + p2*p3 + p2*p4 + p3*p4);
c=(- p1*p2*p3 - p1*p2*p4 - p1*p3*p4 - p2*p3*p4);
d= p1*p2*p3*p4;
Ts = 0.008;
dotnum = (length(t)-1)/t(length(t));
syms   s z ka
ztf1 = subs(d*(1+ka*s)/(s^4+a*s^3+b*s^2+c*s+d), s, 2*(z-1)/(Ts*(z+1)) );
ztf2 = simplifyFraction(ztf1,'Expand', true);
[N,D] = numden(ztf2);%vpa(partfrac(sumoutpos),6)
ztfnum = collect(N,z);
ztfden = collect(D,z);
cnum = coeffs(ztfnum, z);  % cnum(i+1)=z^(i)
cden = coeffs(ztfden, z);  % cden(i+1)=z^(i)
y1=0;
y2=0;
y3=0;
y4=0;
ISE_TS = 0;
for n=5:Ts*dotnum:length(t)*(length(t)/length(t))
    y = 1/cden(5)*(-cden(4)*y1-cden(3)*y2-cden(2)*y3-cden(1)*y4+cnum(5)*pos(n)+cnum(4)*pos(n-1)+cnum(3)*pos(n-2)+cnum(2)*pos(n-3)+cnum(1)*pos(n-4));
    %y = 1/d5*(-d4*y1-d3*y2-d2*y3-d1*y4+n5*pos(n)+n4*pos(n-1)+n3*pos(n-2)+n2*pos(n-3)+n1*pos(n-4));
    y4 = y3;
    y3 = y2;
    y2 = y1;
    y1 = y;
    ISE_TS = ISE_TS+(pos(n)-y)^2; 
end

ISE = ISE_TS*Ts;
ISE_sim = vpa(partfrac(ISE),6)
k = 1.8e4;
subs(ISE_sim, ka, 20/k)





