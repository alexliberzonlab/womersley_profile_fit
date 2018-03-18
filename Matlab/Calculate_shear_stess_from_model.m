%Calculating shear stress  from womersley model velocity profiles using spline fitting .
% Shear stresses are calculated as dynamic viscosity*grad(u)/grad(y)[dyne/cm^2],
% miu=0.0351g/cms.
%Cubic spline smoothinf fitting is applied to the data to calculate tau_f.





clear all; close all;


load('model_velocity_profiles_ex3_last.mat')

col=['b';'k';'r';'g';'y';'c';'m'];

T=[4 4 4 2 2 2 1.2 1.2];

for n=1:8;
    
    y=model_velocity(n).y;
    u=model_velocity(n).u;
    t=model_velocity(n).t;
    
    
    ymax=0.85;%cm
    
    y=y*ymax;
    
    phase_num=length(t);
    
    du=zeros(length(y),phase_num);
    tau=zeros(length(y),phase_num);
    
    
    for i=1:1:phase_num
        
        du(:,i)=gradient(u(:,i));
        dy=gradient(y);
        
        tau(:,i)=0.0351*du(:,i)./dy;
        
        
    end
    
    
    tau_f=zeros(length(y),length(t));
    
    
    for k=1:length(y)
        values=csaps(1:length(tau(k,:)),tau(k,:),.01,1:length(tau(k,:)));
        tau_f(k,:)=values;
        
    end
    
    
    model_velocity(n).tau=tau_f;
    
    
end

% save model_shear_stress_ex3_last.mat model_velocity



