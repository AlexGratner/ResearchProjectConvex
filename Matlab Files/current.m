clc; clear all; close all;

load('roosload.mat');
numHarm = 30;                                % Number of harmonies
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

% for i=i:l
%     [amp_static(i) ang_static(i)] = staic_tf (0.03, 0.0795, 0.05, 3000, 84, 150, 150, 150, 150, harmonics(i));
% end
for i=1:l
    [amp_static ang_static] = static_tf (0.03, 0.0795, 0.05, 3000, 84, 150, 150, 150, 150, harmonics(i));
        valuescl(i,:) = amp_static*amplitude(i)*cos(2*pi*harmonics(i)*t + (phase(i) + ang_static))'; % Y_harmonic
end
% for i=1:l
%         valuescl(i,:) = amplitude(i)*cos(2*pi*harmonics(i)*t + (phase(i)))'; % Y_harmonic
% end


staticsum = sum(valuescl,1);                       % Sum of reference harmonics
static_max = max(staticsum);

for i=1:l
    [amp_dynamic ang_dynamic] = dynamic_tf (0.03, 0.0795, 0.05, 3000, 84, 150, 150, 150, 150, harmonics(i));
        values(i,:) = amp_dynamic*amplitude(i)*cos(2*pi*harmonics(i)*t + (phase(i) + ang_dynamic))'; % Y_harmonic
end
dynamicsum = sum(values,1); 
dynamic_max = max(dynamicsum);

dimension = dynamic_max/static_max;
plot(t, dynamicsum,'b'); 
hold on; 
plot(t, staticsum,'g'); 