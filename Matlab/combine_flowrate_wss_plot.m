% Manually removed the Womersley model data from the original Dikla figures
% then made the extract_data function to get the data of the measured
% and interpolated values in separate MAT files
% now we can plot it. 

% this file creates:
% fig5_gharib_3-8.fig
% fig5_gharib_9-16.fig
%

wss = load('/Users/alex/Dropbox/iditdiklaalexpaper/WSS_3-8.mat');
q = load('/Users/alex/Dropbox/iditdiklaalexpaper/Q_3-8.mat');



figure
for i = 1:6
    subplot(2,3,i);
    [ax,h1,h2] = plotyy(q.data(i).meas_x{1},q.data(i).meas_y{1},...
        wss.data(i).meas_x{1},wss.data(i).meas_y{1});
    set(ax(1),'ylim',[-.5 3],'ytickmode','auto');
    if i == 1 || i == 4, ylabel(ax(1),'Q, lpm'); end
    set(ax(2),'ylim',[-2 4],'ytickmode','auto');
    set(h2,'LineStyle','--');
    if i == 3 || i == 6, ylabel(ax(2),'WSS, dyne/cm^2'); end 
    xlabel(ax(1),'t/T')
    
end
    
    
    
wss = load('/Users/alex/Dropbox/iditdiklaalexpaper/WSS_9-16.mat');
q = load('/Users/alex/Dropbox/iditdiklaalexpaper/Q_9-16.mat');



figure
for i = 1:8
    subplot(3,3,i);
    [ax,h1,h2] = plotyy(q.data(i).meas_x{1},q.data(i).meas_y{1},...
        wss.data(i).meas_x{1},wss.data(i).meas_y{1});
    set(ax(1),'ylim',[-.5 1],'ytickmode','auto');
    if i == 1 || i == 4, ylabel(ax(1),'Q, lpm'); end
    set(ax(2),'ylim',[-2 4],'ytickmode','auto');
    set(h2,'LineStyle','--');
    if i == 3 || i == 6, ylabel(ax(2),'WSS, dyne/cm^2'); end 
    xlabel(ax(1),'t/T')
    
end
    


    

