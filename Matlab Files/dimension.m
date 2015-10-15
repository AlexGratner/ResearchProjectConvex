function dimen_factor = dimension (kt, rm, lm, ll, kl, n, p1, p2, p3, p4, t, pos)

numHarm = 90;                               % Number of harmonies
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
    [amp_dynamic ang_dynamic] = dynamic_tf (rm, lm, ll, kl, n, p1, p2, p3, p4, harmonics(i));
        values(i,:) = amp_dynamic*amplitude(i)*cos(2*pi*harmonics(i)*t + (phase(i) + ang_dynamic))'; % Y_harmonic
end
dynamicsum = sum(values,1); 
dynamic_trms = sqrt(mean((kt*dynamicsum).^2));

Cm = 91350;
static_trms = Cm*lm*rm^2.5;
dimen_factor = (dynamic_trms/static_trms)^(3/3.5);
end