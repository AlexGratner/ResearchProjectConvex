function [amp, ang] = dynamic_tf (rm, lm, ll, kl, n, p1, p2, p3, p4, fre)

Cmj = 1056;
inertia_ratio = 10/9; 
steel_den = 8000;
eta = 1;
d = 20;

jl = 1.1;
jg = 0;
jm = inertia_ratio * Cmj * lm * rm^4;
js = kl*steel_den*ll^2/75e9;
kt=0.182;
omega = 2*pi*fre;

u_dynamic = (p1*p2*p3*p4*(omega*j)^2*(4*kl*(jl+n^2*eta*(jm+jg+js))+4*d*(jl+n^2*eta*(jm+jg+js))*omega*j+(2*jm+js)*(2*jl+n^2*eta*(2*jg+js))*(omega*j)^2))/(kt*n*eta*kl*(p1+omega*j)*(p2+omega*j)*(p3+omega*j)*(p4+omega*j)*4);;
amp = abs (u_dynamic);
ang = angle(u_dynamic);

end


