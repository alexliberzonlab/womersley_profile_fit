% WomersleyQWave.m
%
%		Calculates Womersley solution for fully developed, periodic laminar flow in a straight tube
%
%       Calculates velocity profiles, and pressure, wall shear stress, and centerline velocity waveforms
%       Also produces an animation of the flow profile.
%
%       Instructions:  Prior to executing this routine (WomersleyQWave), run
%       FSDecomposeCheck, which computes Fourier Series amplitudes and phases
%       from a flow waveform, in file femoralraw.txt, and stores them in
%       array KQ.
%
% 		Developed for Biofluid Mechanics, Mechanical Engineering and Mechanics, Drexel University,
%			Philadlephia PA, 19104.
%
%		Copyright David Wootton,  March, 2003
%
clear all; close all;
clc;
load('./scaled_velocity_ex2.mat')
load('./qspline_ex2.mat')
load('./q_ex2.mat')
load('C:\Documents and Settings\Dikla\My Documents\MATLAB\PIV 2th ex27.4-29.4\inverse womersley model\trigger_points.mat')
ind=[];

for nExp = 1:6

% For the pulsatile wave forms
%[KQ KQ0 KQmag KQphase Q Qf Qraw T j j32 k n nf nt nt1 t t1 w]=...
%    FSDecomposeCheck();
% For the half-sine wave forms

period=[2 1 4 2 1 2];
nf=32;% Number of frequencies for resynthesized flow wave.
T=period(nExp);	% period, in seconds
trigger=trigger_points(nExp).normalized_time;
[KQ KQ0 KQmag KQphase Q Qf Qraw T j k n nf nt nt1 t1 w]=...
    FSDecomposition(qspline,nExp,T,nf,q.(sprintf('set%g',nExp)),trigger);

a =   1.8/2;		% radius, in cm, femoral artery = 0.53, aorta = 1.11, coronary = 0.186
% KQ is from workspace, as is KQ0 (DC component of flow)
nu = 0.0321;	% kinematic viscosity, in cm2/s
rho = 1.097;	% density, g/ml

%miu=0.0351 dynamic viscosity [g/cm*s]

w0 = 2*pi/T;			    % fundamental radian frequency
alpha0 = a*sqrt(w0/nu);	    % Womersley parameter, fundamental frequency

ny = 101;		% # points in profile (over diameter)


nt = 100;		% # time points to calculate over period - 7 points


nf = size(KQ);  nf = nf(1);		% # of nonzero frequencies
temp2 = ones(1,nt);  % temp2 = temp2(1,:);
mean_u = mean(Q)*100 /(6*(pi*a*a));% calculate mean velocity [cm/s]
Re = mean_u*(2*a)/nu;% calculate Reynolds number

y = [-1: 2/(ny-1): 1]';				% column vector of y's


t = 0: T/nt : T; t = t(1:nt);	% row vector of times
trigger_time=T*trigger;

for p=1:length(trigger_time)
    ind(p)=max(find(t>=(trigger_time(p)-0.02) &t<=(trigger_time(p)+0.02)));
end

fh = @(x,xdata,T,nu,a,KQ,KQ0,nf,rho)(womersley_profile_v5(x,xdata,T,nu,a,KQ,KQ0,nf,rho));


options = optimset('Display','iter');
options = optimset('MaxFunEvals',2000);
options = optimset('TolFun',1e-7);
options = optimset('TolX',1e-7);

lb = [0];
ub = [T];


t1 = zeros(length(trigger),1);

 figure;hold on;

for i = 1:length(trigger)


    temp=d_ex2(nExp).(sprintf('t%g',i-1));
    
    if(~isempty(temp))
  
    y1 = temp(:,1)/9; %r/R
    u1 = 100*temp(:,2);
    
    [tmp,resnorm] = lsqcurvefit(fh,t(ind(i)),y1,u1,lb,ub,options,T,nu,a,KQ,KQ0,nf,rho);%finds tje time stamp for which the compatibilty between model and measured profiles is maximal.
    %ts(i,:) = tmp;% model time stamp for a given set 
    col=['b';'k';'r';'g';'y';'c';'m'];
    plot(y1,u1,'.','Color',sprintf('%c',col(i)));hold on;
    plot(y,fh(tmp,y,T,nu,a,KQ,KQ0,nf,rho),'--','Color',sprintf('%c',col(i)));
    title(sprintf('set %g alpha=%g R=%g',nExp,alpha0,a));
  
    hold on;
grid on;

% velocity_profiles(nExp).(sprintf('t%g',i-1))=[y1 fh(tmp,y1,T,nu,a,KQ,KQ0,nf,rho)];

%velocity_profiles(nExp).(sprintf('t%g',i-1))=[y1 fh(tmp,y1,T,nu,a,KQ,KQ0,nf,rho)];
% [tau] = Copy_of_womersley_v5(tmp,y1,T,nu,a,KQ,KQ0,nf,rho);%calculating wall shear stress 
% Tau_wall(nExp,i)=tau;%Wall Shear Stress (dyn/cm^2)
 ts(nExp,i)=tmp;
    end
end

%velocity_profiles(nExp).phases=(360*t1)/T';


end
%save model_velocity.mat velocity_profiles
%save model_WSS.mat Tau
% w0 = 2*pi/T;			    % fundamental radian frequency
% alpha0 = a*sqrt(w0/nu);	    % Womersley parameter, fundamental frequency
% w = w0*[1:nf];			    % array of radian frequencies
% alpha = a*sqrt(w/nu);	    % array of Womersley parameters
% K = 0*KQ;					% initialize pressure gradient array
% Q = KQ0*temp2;				% initialize flow rate array with DC flow component
% Pp = -8*nu*Q/(pi*a^4);	% initialize pressure gradient array with DC component
% Tau_steady = -a*Pp/2;	% initialize wall shear stress array with DC component
% Tau = Tau_steady;
% j = sqrt(-1);
% j32 = j^1.5;
%
% u = (2*KQ0/(pi*a^2))*((1 - y.*y)*temp2);
%
% for n = 1:nf		% Sum over all frequencies
%    K(n) = KQ(n)*j*w(n)/(pi*a^2*(1 - 2*besselj(1,j32*alpha(n))/(j32*alpha(n)*besselj(0,j32*alpha(n)))));
%    Ktau(n) = rho*K(n)*a*besselj(1,j32*alpha(n))/(j32*alpha(n)*besselj(0,j32*alpha(n)));
%    for k = 1:nt		% Calculate for each point in time
%       Q(k) = Q(k)+real(KQ(n)*exp(j*w(n)*t(k)));
%       Pp(k) = Pp(k)+real(K(n)*exp(j*w(n)*t(k)));
%       Tau(k) = Tau(k)+real(Ktau(n)*exp(j*w(n)*t(k)));
%       for m = 1:ny		% Calculate for each spatial location y
%          u(m,k) = u(m,k) + real((K(n)*a^2/(nu*alpha(n)^2*i))*(1-besselj(0,j32*alpha(n)*y(m))/besselj(0,j32*alpha(n)))*exp(i*w(n)*t(k)));
%       end	% y loop
%    end	% t loop
% end	% f loop


%k=[1 13 23 38 51 63 76];
%col=['b';'k';'r';'g';'y';'c';'m'];

%subplot(2,4,nExp);
%for j=1:7
%subplot(2,4,j);
%for i=1:7
    %b = fh(t1(i),y,T,nu,a,KQ,KQ0,nf,rho);
    %plot(y,b,'LineStyle','--','Color',sprintf('%c',col(i)));hold on;
    %temp=d(nExp).(sprintf('t%g',i-1));
    %plot(temp(:,1)/max(temp(:,1)),100*temp(:,2),'LineStyle','none','Marker','o','Color',sprintf('%c',col(i)));
    %hold on;
   
%end

%xlabel('r/R');
%ylabel('u [cm/s]');
%title(sprintf('set %g alpha=%g R=%g',nExp,alpha0,a));
%grid on;
%axis([-1 1 -2 15]);

%end
%umin = min(min(u));  umax = max(max(u));



%Make an animated velocity profile plot
%{
clear M;
for n = 1:nt
   h1 = figure(3);
   plot(u(:,n),y); axis([1.1*umin 1.1*umax -1 1]);
   xlabel(['U(y,t=' num2str((n-1)*T/nt) ')']); ylabel('y'); title('Velocity Profiles');
   M(n) = getframe(h1);
end

moviereps = 30;	% number of repetions of the movie
speedratio = 0.5;	% ratio of movie speed to actual speed
%movie(M,moviereps,speedratio*nt/T);
movie2avi(M,'velocity_profile_femoral.avi','quality',100) % record the animation as AVI files
%}

%subplot(1,2,2);
%Plot selected velocity profiles on a static plot
%t = [0: T/nt : T];
%hold on;
%for i = 1:10:100
    %plot(y,fh(t(i),y,T,nu,a,KQ,KQ0,nf,rho));
%end
%legend('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%');
%title(sprintf('Velocity Profiles set%g, alpha=%g',nExp,alpha0)); xlabel('U(y;percent of period) [cm/s]'); ylabel('r/R');
%grid on;
%axis([-1 1 -2 20]);

%hold off
%figure;
%Plot time waveforms of flow, pressure gradient, wall shear stress, and centerline velocity
%subplot(2,2,1); plot(t,Q); xlabel('time (s)'); ylabel('Flow Rate (ml/s)');title(sprintf('set %g',nExp));
%subplot(2,2,2); plot(t,Tau,t,Tau_steady); xlabel('time (s)'); ylabel('Wall Shear Stress (dyn/cm^2)');title(sprintf('set %g',nExp));
%subplot(2,2,3); plot(t,Pp); xlabel('time (s)'); ylabel('Pressure Gradient (dyn/cm^3)');title(sprintf('set %g',nExp));
%subplot(2,2,4); plot(t,u((ny+1)/2,:)); xlabel('time (s)'); ylabel('Centerline Velocity (cm/s)');title(sprintf('set %g',nExp));



% save model_time_stamp.mat ts