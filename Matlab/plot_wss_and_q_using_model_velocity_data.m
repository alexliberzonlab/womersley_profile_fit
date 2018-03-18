% Womersley model results: velocity profiles at 
% 200 points in r and 400 points in time

%% Run 3

N = 8;

load('model_velocity_profiles_ex3_last.mat')

% plot WSS for the 8 experiments from run 3 (9-16)

% colors = {'r','g','b','k','m','y','c','b--'};

% colors = colorGray(8); % for gray scale print
colors = distinguishable_colors(N,'w');

figure, hold on
for i = 1:N
    [r,c] = size(model_velocity(i).u);
    dudy = zeros(size(c,1));
    for j = 1:c
        tmp = gradient(model_velocity(i).u(:,j),model_velocity(i).y);
        dudy(j) = 0.0351*tmp(1); % first, most left point
    end
    plot(model_velocity(i).t/model_velocity(i).t(end),dudy,'color',colors(i,:));
end

xlabel('t/T')
ylabel('WSS, [dyne/cm^2]')
legend({'9','10','11','12','13','14','15','16'})


% load('./qspline_ex3_last.mat')
load('./q_ex3_last.mat')
t2 = [0 0.125 0.25 0.375 0.5 0.625 0.75]; % normalized phase 

figure, hold on
for i = 1:N
    [r,c] = size(model_velocity(i).u);
    q1 = zeros(size(c,1));
    for j = 1:c
        q1(j) = trapz(model_velocity(i).y, squeeze(model_velocity(i).u(:,j)));
    end
    
    % plot(linspace(0,1,length(qspline(i,:))), qspline(i,:),'--');
    errorbar(t2,q(i,:),2*ones(length(q(i,:)),1),'color',colors(i,:),'Marker','o','LineStyle','none');
    plot(model_velocity(i).t/model_velocity(i).t(end),q1,'color',colors(i,:));
end

xlabel('t/T')
ylabel('Q [ml/s]' )
% legend({'9','10','11','12','13','14','15','16'})




%% Run 2

% Womersley model results: velocity profiles at 
% 200 points in r and 200 points in time

N = 6;

load('model_velocity_profiles_ex2.mat')

% plot WSS for the 6 experiments from run 2 (3-8)

% colors = {'r','g','b','k','m','y','c','b--'};

% colors = colorGray(8); % for gray scale print
colors = distinguishable_colors(N,'w');

figure, hold on
for i = 1:N
    [r,c] = size(model_velocity(i).u);
    dudy = zeros(size(c,1));
    for j = 1:c
        tmp = gradient(model_velocity(i).u(:,j),model_velocity(i).y);
        dudy(j) = 0.0351*tmp(1); % first, most left point
    end
    plot(model_velocity(i).t/model_velocity(i).t(end),dudy,'color',colors(i,:));
end

xlabel('t/T')
ylabel('WSS, [dyne/cm^2]')
% legend({'9','10','11','12','13','14','15','16'})
legend({'3','4','5','6','7','8'});


% load('./q_spline_ex2_last.mat')
load('./q_ex2_last.mat')
load('./trigger_points.mat')
period=[2 1 4 2 1 2 ];

exp = fieldnames(q); % names of the exps are set1, set2, ....

% t2 = [0 0.125 0.25 0.375 0.5 0.625 0.75]; % normalized phase 

figure, hold on
for i = 1:N
    [r,c] = size(model_velocity(i).u);
    q1 = zeros(size(c,1));
    for j = 1:c
        q1(j) = trapz(model_velocity(i).y, squeeze(model_velocity(i).u(:,j)));
    end
    
    % plot(linspace(0,1,length(qspline(i,:))), qspline(i,:),'--');
    t2 = trigger_points(i).normalized_time; 
    q2 = q.(exp{i}); 
    
    
    errorbar(t2,q2,2*ones(length(q2),1),'color',colors(i,:),'Marker','o','LineStyle','none');
    plot(model_velocity(i).t/model_velocity(i).t(end),q1,'color',colors(i,:));
end

xlabel('t/T')
ylabel('Q [ml/s]' )
legend({'3','4','5','6','7','8'});



