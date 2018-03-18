
% This code apply a spline fit to the flow rate values calculated from the
% measured velocity profiles. The objective is to create a  vector with
% multiple data points that can be decomposed into fourier coefficient in
% the following analysis process. The number of phase in icreased from 7 to
% 100. 
%Input:q.mat- a 8x7 matrix of instantaneous flow rate values. when rows-ex. num, cols:phase. 
%Output:qspline. mat- a 8x100 matrix. same arrengement. 

clear all; %close all; clc; 

load('C:\Documents and Settings\Dikla\My Documents\MATLAB\PIV 2th ex27.4-29.4\inverse womersley model\q_ex2.mat')
load('C:\Documents and Settings\Dikla\My Documents\MATLAB\PIV 2th ex27.4-29.4\inverse womersley model\trigger_points.mat')

T=[2 1 4 2 1 2];

%qspline=zeros(length(q),100);
figure;


for i=1:6
%Because the number and phase of the triggering points in this experiment are not uniform,the
%triggering points data ws saved in a special struct trigger_points.

t=T(i)*trigger_points(i).normalized_time;

%Because the low number of data points the period was dubbeled to impose periodic conditions

t_n=[t (t+T(i))];

%Instantaneous flow rate values

temp=q.(sprintf('set%g',i));

%replicating values
temp_n=repmat(temp,1,2);
tt = linspace(0,T(i),100);%0:0.01:(T(i)-0.01);
tt_n=[tt(1:end-1) (tt+T(i))];
pp = csaps(t_n,temp_n,1);%Best order of fit is 6 % 
qq = fnval(pp,tt_n);

subplot(2,3,i)
plot(t/T(i),temp,'o','Color','r');hold on;
plot(tt/T(i),qq(1:length(tt)),'--');
title(sprintf('set%g',i));
axis([0 1 -5 40]);
qspline(i,:)=qq(1:length(tt));
t_int(i,:)=tt;
end

% save qspline_ex2.mat qspline
% save qspline_time_interval.mat t_int
