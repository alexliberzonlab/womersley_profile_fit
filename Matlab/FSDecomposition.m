% function [KQ KQ0 KQmag KQphase Q Qf Qraw T j j32 k n nf nt nt1 t t1 w]=...
%     FSDecomposeCheck(EX3)
function [KQ KQ0 KQmag KQphase Q Qf Qraw T j k n nf nt nt1 t1 w]=...
    FSDecomposeCheck(EX3,nExp,T,nf,q)

%clear all; close all; clc;

% FSDecomposeCheck.m
%
%   Makes a Fourier series decomposition of a waveform, from a tab-delimited test file
%       named 'femoralraw.txt', containing equally spaced (time, flow) data
%		  Plots Fourier series with orginal data for quality control.
%
%   Use this routine to prepare data for a Womersley solution for real flow waveform
%       (ie as a preprocessing step prior to running WomerslyQWave.m)

% read flow rate waveform from file


dt = linspace(0,T,100);
dt1=dt(1:(length(dt)-1));
Qraw=EX3(nExp,1:(length(EX3)-1));%flow rate[ml/s]
 				
KQ = ones(1,nf)'; w = ones(1,nf)';% initialize arrays
nt = size(dt);  nt = max(nt);

Qf = fft(Qraw);% Calculate discrete Fourier Series expansion of flow waveform

% Recombine FFT coefficients to get Fourier Series Coefficients; compute frequencies
KQ0 = real(Qf(1))/nt;

j = sqrt(-1);

for n = 1:nf
   KQ(n) = 2*Qf(n+1)/nt;
   w(n)	= 2*pi*n/T;
end
                                                               
KQmag = [KQ0; abs(KQ)];
KQphase = [0.0; angle(KQ)];

t1 =[0:0.001*T:(T-(dt1(2)-dt1(1)))]';
nt1 = size(t1);  nt1 =max(nt1);                                                                                                                                                                                                                             max(nt1);

Q = KQ0*ones(1,nt1)';

for n = 1:nf
   for k = 1:nt1
      Q(k) = Q(k) + real(KQ(n)*exp(j*w(n)*t1(k)));
   end
end

%fprintf(1,['Mean volume flow rate = ' num2str(mean(Q)) '(ml/s) \n']);
t2=T*[0 0.125 0.25 0.375 0.5 0.625 0.75];
figure;
plot(dt1,Qraw,t1,Q);hold on; plot(t2,q(nExp,:),'LineStyle','none','Marker','o'); xlabel('Time (s)'); 
ylabel('Q (ml/s)'); 
legend('Original','Fourier Series');
title(sprintf('set %g' ,nExp));