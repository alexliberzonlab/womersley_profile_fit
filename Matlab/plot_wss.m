%Plot WSS model vs. measurements

clear all; 


load('tau_wall_measured_ex3_last.mat')
load('model_shear_stress_ex3_last.mat')

T=[4 4 4 2 2 2 1.2 1.2];

figure;

for n=1:8
    
t=[0 0.125 0.25 0.375 0.5 0.625 0.75];
dt=model_velocity(n).t;
tau=model_velocity(n).tau;

wss_measured=tau_wall_measured(:,n);

subplot(3,3,n);
plot(dt/T(n),tau(1,:),'LineStyle','--','Color','b');hold on; 
plot(t,wss_measured,'LineStyle','none','Marker','.','Color','r');
% axis([0 1 -1 5]);
ylabel('dyne/cm2');
xlabel('t/T');
title(sprintf('set %g',n));
end 