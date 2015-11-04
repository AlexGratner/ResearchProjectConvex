s=tf('s');
p1=-150;
p2=-150;
p3=-150;
p4=-150;
k=5000;%
jl=1.1;
jm=10/9*Cmj*l*r^4;
jg=0;
l=0.05;
rho=7800;
js=k*rho*l^2/75e9;
eta=1;
d=20;
kt=0.182;

Gur = (p1*p2*p3*p4*s^2*(4*k*(jl+n^2*eta*(jm+jg+js))+4*d*(jl+n^2*eta*(jm+jg+js))*s+(2*jm+js)*(2*jl+n^2*eta*(2*jg+js))*s^2))/(kt*n*eta*k*(p1-s)*(p2-s)*(p3-s)*(p4-s));
Gurs = (p1*p2*p3*p4*s^2*4*(jl+n^2*eta*(jm+jg+js)))/(kt*n*eta*(p1-s)*(p2-s)*(p3-s)*(p4-s));




