
% This code apply a spline fit to the flow rate values calculated from the
% measured velocity profiles. The objective is to create a  vector with
% multiple data points that can be decomposed into fourier coefficient in
% the following analysis process. The number of phase in icreased from 7 to
% 100. 
%Input:q.mat- a 8x7 matrix of instantaneous flow rate values. when rows-ex. num, cols:phase. 
%Output:qspline. mat- a 8x100 matrix. same arrengement. 

clear all; %close all; clc; 

load('C:\Documents and Settings\Dikla\My Documents\MATLAB\ex 3 3.8.2010\piv\analytical model\q1.mat')

T=[4 4 4 2 2 2 1.2 1.2];

%qspline=zeros(length(q),100);
figure;

n=[4 6 4 6 6 4 4 6];%Splin oreder.Determined by the best match to the data points as seen by the user. 
for i=1:8

t=T(i)*[0 0.125 0.25 0.375 0.5 0.625 0.75];
temp=q(i,1:7);
temp_n=repmat(temp,1,2);
t_n=[t (t+T(i))];

tt = linspace(0,T(i),100);%0:0.01:(T(i)-0.01);
tt_n=[tt(1:end-1) (tt+T(i))];

pp = csaps(t_n,temp_n,1);%Best order of fit is 6 % 
qq = ppval(pp,tt_n);
subplot(2,4,i)
plot(t/T(i),temp,'o','Color','r');hold on;
plot(tt/T(i),qq(1:length(tt)),'--')
title(sprintf('set%g ',i));
qspline(i,:)=qq(1:length(tt));
end

save qspline1.mat qspline