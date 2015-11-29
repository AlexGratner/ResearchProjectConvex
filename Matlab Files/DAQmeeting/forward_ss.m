clc
clear all
load('roosload.mat');

% sys = ss(a,b,c,d) 
% syms k;
% a = [0 1 0;-1 0 1;0 0 1] ;
% b = [0;0;1];
% c = [1 0 0];
% d = [0];
% sys = ss(a,b,c,d)
% Ts = 0.1;
% sysd = c2d(sys,Ts,'ForwardEuler')
% 506250000*(1+Ka*s)/(s^4 + 600 s^3 + 135000 s^2 + 1.35e07 s + 5.063e08);


syms Ka;
% Ka = 20/1.8e4;
Ts = 0.01;
dotnum = (length(t)-1)/t(length(t));
x1=0;
x2=0;
x3=0;
x4=0;
x1_next=0;
x2_next=0;
x3_next=0;
x4_next=0;
ISE_ts=0;
tm=0;
for n=1:Ts*dotnum:length(t)
    y = tm;
    ISE_ts = ISE_ts+(pos(n)-y)^2; 
    x1_next = x1+Ts*x2;
    x2_next = x2+Ts*x3;
    x3_next = x3+Ts*x4;
    x4_next = -506250000*Ts*x1-1.35e07*Ts*x2-135000*Ts*x3+(-600*Ts+1)*x4+Ts*pos(n);
    tm = 506250000*x1_next+506250000*Ka*x2_next;
    x1 = x1_next;
    x2 = x2_next;
    x3 = x3_next;
    x4 = x4_next; 
end
ISE = ISE_ts*Ts;
ISE_sim = vpa(partfrac(ISE),6)
k = 1.8e4;
subs(ISE_sim, Ka, 20/k)