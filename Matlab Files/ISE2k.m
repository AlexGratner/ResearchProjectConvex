function[k_lower,k_upper] = ISE2k(ISE, pos, t) 
% load('roosload.mat');
poles = [150 150 150 150];
numHarm = 5;                                % Number of harmonies
duration = t(end);                          % Get end from input
fs = (length(t)-1)/duration;            
noSamples = length(t);                  
f = 0 : fs/noSamples : fs - fs/(noSamples);
x_fft = fft(pos);
pw = x_fft.*conj(x_fft)/noSamples;
m = 2.*abs(x_fft)/noSamples;
p = unwrap(angle(x_fft)); 
pw(noSamples/2:noSamples)=zeros(1,floor(noSamples/2+1));
[Xs,IX]=sort(pw,'descend');
f_sorted = f(IX);
m_sorted = m(IX);
p_sorted = p(IX);
m_sorted(f_sorted==0) = 0.5*m_sorted(f_sorted==0);
rms_array = (1/sqrt(2)*m_sorted(1:numHarm)).^2;
rms_value = sqrt(sum(rms_array));
for i=1:numHarm
   x_app(i,:) = m_sorted(i).*cos(2*pi.*f_sorted(i).*t + p_sorted(i));
end
x_approx = (sum (x_app,1))';
rms_approx = sqrt(integrator(pos.^2,t)/t(end));
             harmonics = f_sorted (1:numHarm);
             phase = p_sorted (1:numHarm);
             amplitude = m_sorted (1:numHarm);
l = length(harmonics);
for i=1:l
        valuescl(i,:) = amplitude(i)*cos(2*pi*harmonics(i)*t + (phase(i)))'; % Y_harmonic
end



refsum = sum(valuescl,1);                       % Sum of reference harmonics
%plot(t, refsum);



%  c = real((-150-2*pi*harmonics(i)*1i)^4)  
%  d = imag((-150-2*pi*harmonics(i)*1i)^4)

syms  k t0
p1=-poles(1);
p2=-poles(2);
p3=-poles(3);
p4=-poles(4);

for i=1:l
   refsumpos(i) = amplitude(i)*cos(2*pi*harmonics(i)*t0 + (phase(i)))'; % Y_harmonic 
end

a = p1*p2*p3*p4;

d0= 20;                 %damping in the shaft

for i=1:l
        c = real((p1-2*pi*harmonics(i)*1i)*(p2-2*pi*harmonics(i)*1i)*(p3-2*pi*harmonics(i)*1i)*(p4-2*pi*harmonics(i)*1i));
        d = imag((p1-2*pi*harmonics(i)*1i)*(p2-2*pi*harmonics(i)*1i)*(p3-2*pi*harmonics(i)*1i)*(p4-2*pi*harmonics(i)*1i));
        amp = amplitude(i);
        ef = harmonics(i);
        fai = phase(i);
        outpos(i) = a/(c^2+d^2)*amp*(cos(2*pi*ef*t0+fai)*(k*c+2*pi*ef*d*d0)/k-sin(2*pi*ef*t0+fai)*(2*pi*ef*d0*c-k*d)/k); %output signal for the i-th input 
        
end

sumoutpos = simplify(sum(outpos,1));

err_squared = (sum(refsumpos) - sum(sumoutpos)).^2;
ISE_value = int(err_squared, t0, t(1), t(end));
ISEbyk = vpa(partfrac(ISE_value), 4) ;
syms k_inv;
ISEbyk_inv = subs(ISEbyk, k, 1/k_inv );
%coe = coeffs(ISEbyk_inv, k_inv);
ksolution = solve(ISEbyk_inv==ISE, k_inv);
k_lower = 1/max(ksolution);
k_upper = 1/min(ksolution);

end


%  rysquare(j)=(pos(round((length(t)-1)/t(end)*t0+1))-sumoutpos)^2;
%end
 
%alex
%int((sumoutpos-pos((length(t)-1)/t(end)*t0+1))^2,t0,0,1.8)
%

%rysquareintegrate = sum(rysquare)/secpoints*t(end);
%vpa(partfrac(rysquareintegrate),6)

% sumoutpos = simplify(sum(outpos));  %the sun of all i signals
% vpa(partfrac(sumoutpos),6) ;        % get position output versus k
% vpa(subs(sumoutpos,1.8e5),4)  ;     % get position output with chosen k 
% pos((length(t)-1)/t(end)*t0+1) ;    % get the ref input at t0,say r(t0)