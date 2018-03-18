% figure 5 Gharib

% Womersley model results: velocity profiles at
% 200 points in r and 400 points in time


load('model_velocity_profiles_ex3_last.mat')
load('./qspline_ex3_last.mat')
load('./q_ex3_last.mat')


t2 = [0 0.125 0.25 0.375 0.5 0.625 0.75]; % normalized phase

% plot WSS for the 8 experiments from run 3 (9-16)

% colors = colorGray(8); % for gray scale print
colors = distinguishable_colors(8,'w');


for i = 1:8
    figure
    [r,c] = size(model_velocity(i).u);
    dudy = zeros(size(c,1));
    q1 = zeros(size(c,1));
    for j = 1:c
        tmp = gradient(model_velocity(i).u(:,j),model_velocity(i).y);
        dudy(j) = 0.0351*tmp(1); % first, most left point
        q1(j) = trapz(model_velocity(i).y, squeeze(model_velocity(i).u(:,j)));        
    end
    [ax,h1,h2] = plotyy(model_velocity(i).t/model_velocity(i).t(end),dudy,...
        model_velocity(i).t/model_velocity(i).t(end),q1);
    
    axes(ax(2));
    hold on
    % errorbar(t2, q(i,:), 2*ones(length(q(i,:)),1),'o');
    errorbar(t2, q(i,:), .1 * max(q(i,:)) * ones(7,1),'o');
    
    model_velocity(i).tau = dudy;
    model_velocity(i).q = q1; 
    
    
    xlabel('t/T')
    ylabel(ax(1),'WSS, [dyne/cm^2]')
    ylabel(ax(2), 'Q [ml/s]' )
    title(int2str(i+8));
end

% save model_velocity_ex3 model_velocity





