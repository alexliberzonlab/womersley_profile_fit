% let's try to create a modified figure 6b
% Ratio of Max negative to Maximum positive WSS


alpha0 =  [9,  13,  6,  9, 13,  9,  6, 6, 6,    8, 8, 8,    11, 11];
markers = {'o','s','<','o','s','o','<','<','<','p','p','p','d','d'};

PI = [2.4, 1.6 3.3 2.6 1.7 1.5 4 3.9 8.9 2.4 3.2 6.4 1.8 4.8];
Re = [633, 732, 490, 477, 496, 592, 190, 189, 69, 329, 195, 80, 335, 90];


%  T = [2 1 4 2 1 2 4 4 4 2 2 2 1.2 1.2 1.2];
% a =  1.7/2; % radius, in cm, femoral artery = 0.53, aorta = 1.11, coronary = 0.186
% nu = 0.0321;	% kinematic viscosity, in cm2/s
% rho = 1.097;	% density, g/ml
% 
% w0 = 2*pi./T;			    % fundamental radian frequency
% alpha0 = round(a*sqrt(w0/nu));	    % Womersley parameter, fundamental frequency

% ratio = zeros(length(PI),1);
% 
% N = 6;
% load('model_velocity_profiles_ex2.mat')
% 
% for i = 1:N
%     [r,c] = size(model_velocity(i).u);
%     dudy = zeros(size(c,1));
%     for j = 1:c
%         tmp = gradient(model_velocity(i).u(:,j),model_velocity(i).y);
%         dudy(j) = 0.0351*tmp(1); % this is WSS
%     end
%     % ratio(i) = length(find(dudy < 0 ))/length(dudy); see figure 6a
%     ratio(i) = min(dudy)/max(dudy);  %see figure 6b
% end
% 
% N = 8;
% load('model_velocity_profiles_ex3_last.mat')
% 
% figure, hold on
% for i = 1:N
%     [r,c] = size(model_velocity(i).u);
%     dudy = zeros(size(c,1));
%     for j = 1:c
%         tmp = gradient(model_velocity(i).u(:,j),model_velocity(i).y);
%         dudy(j) = 0.0351*tmp(1); % first, most left point
%     end
%     ratio(i+6) = min(dudy)/max(dudy);  %see figure 6b
% end
% 
% 
% save figure6b


%% 
load figure6b

figure
hold on
colors = distinguishable_colors(length(unique(PI)));

[~,k] = sort(PI);

for i = 1:length(k)
    ind = find(unique(PI) ==  PI(k(i)));
    plot(Re(k(i)),ratio(k(i)),'bo','MarkerSize',8,'MarkerFaceColor',colors(ind,:),'DisplayName',num2str(PI(k(i))));
end
legend('toggle');
xlabel('Re -> Cardiac Output')
ylabel('Min/Max')



figure
hold on

for i = 1:length(k)
    ind = find(unique(PI) ==  PI(k(i)));
    plot(alpha0(k(i)),ratio(k(i)),'bo','MarkerSize',8,'MarkerFaceColor',colors(ind,:),'DisplayName',num2str(PI(k(i))));
end
hl = legend('toggle');
set(get(hl,'Title'),'String','PI');
xlabel('\alpha_0 -> Heart beat')
ylabel('Min/Max')



figure
hold on
colors = distinguishable_colors(length(ratio));

[~,k] = sort(Re);

for i = 1:length(k)
    plot(PI(k(i)),ratio(k(i)),markers{k(i)},'MarkerSize',8,'MarkerFaceColor',colors(k(i),:),'DisplayName',num2str(Re(k(i))));
end

hl = legend('toggle');
set(get(hl,'Title'),'String','Re');
ylabel('Min/Max')
xlabel('Pulsativity index')
    


figure
hold on
colors = distinguishable_colors(length(unique(alpha0)));

[~,k] = sort(alpha0);

for i = 1:length(k)
    ind = find(unique(alpha0) ==  alpha0(k(i)));
    plot(PI(k(i)),ratio(k(i)),'bo','MarkerSize',8,'MarkerFaceColor',colors(ind,:),'DisplayName',num2str(alpha0(k(i))));
end
hl = legend('toggle');
set(get(hl,'Title'),'String','\alpha_0','FontSize',16);
xlabel('Pulsativity index')
ylabel('Min/Max')




figure
plot(PI,ratio,'o');
