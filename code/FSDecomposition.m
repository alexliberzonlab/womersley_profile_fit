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

if nargin == 0,
Qraw =[-7.7183,-8.2383,-8.6444,-8.8797,-9.6337,-10.5957,-11.8705,-10.0942,-6.2839,-1.1857,2.6043,4.4323,6.1785,7.8211,9.1311,9.9138,10.3447,10.4011,10.2807,9.8951,8.0597,5.6717,2.5232,1.3301,1.4405,1.9094,1.8145,0.8738,0.7055,0.7343,0.7788,0.7495,0.6711,-0.4796,-1.6541,-2.8643,-3.4902,-4.1714,-5.6581,-6.8024];

T = 1;
nf = 10;

q = Qraw;

else


Qraw=EX3(nExp,1:(length(EX3)-1));%flow rate[ml/s]

end
dt = linspace(0,T,length(Qraw));
dt1=dt(1:(length(dt)));



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

disp(length(dt1))
disp(length(Qraw))
disp(length(t1))
disp(length(Q))

figure;
plot(dt1,Qraw,t1,Q);
hold on;
% apparently q(nExp, :) was the experimental data
% and EX3 was the reconstructed Q spline

% plot(t2,q(nExp,:),'LineStyle','none','Marker','o'); 
xlabel('Time (s)'); 
ylabel('Q (ml/s)'); 
legend('Original','Fourier Series');

saveas(gcf,'../results/dt_Q.png')