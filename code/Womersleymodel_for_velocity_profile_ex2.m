
%Calculating velocity profile data for a given dt,dy


clear all; close all;
clc;

load('./q_spline_ex2_last.mat')
load('./q_ex2_last.mat')
load('./trigger_points.mat')




for nExp =1:6

% For the pulsatile wave forms
%[KQ KQ0 KQmag KQphase Q Qf Qraw T j j32 k n nf nt nt1 t t1 w]=...
%    FSDecomposeCheck();
% For the half-sine wave forms

period=[2 1 4 2 1 2 ];
nf=32;% Number of frequencies for resynthesized flow wave.
T=period(nExp);% period, in seconds
trigger=trigger_points(nExp).normalized_time;

[KQ KQ0 KQmag KQphase Q Qf Qraw T j k n nf nt nt1 t1 w]=...
    FSDecomposition(qspline,nExp,T,nf,q.(sprintf('set%g',nExp)),trigger);



a =   1.7/2;		% radius, in cm, femoral artery = 0.53, aorta = 1.11, coronary = 0.186
% % KQ is from workspace, as is KQ0 (DC component of flow)
nu = 0.0321;	% kinematic viscosity, in cm2/s
rho = 1.097;	% density, g/ml

%miu=0.0351 dynamic viscosity [g/cm*s]
dt=0.01; %define time interval between data points.
dy=0.01;%define radial spacing 

[t y u]=calculate_velocity_profile_from_model1_ex2(dt,T,dy,nu,rho,nf,KQ,KQ0,a);

%saving the profiles data into a structure.

model_velocity(nExp).u=u;
model_velocity(nExp).y=y;
model_velocity(nExp).t=t;

 

end

save model_velocity_profiles_ex2.mat model_velocity




