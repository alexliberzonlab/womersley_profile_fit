function [data] = extract_data_from_flowrate_figure(hf)
% saves all the data in the figure of flow rate vs time
% in some structure
% reproduction test
%
% Example:
% 
% uiopen('/Users/alex/Dropbox/iditdiklaalexpaper/Q_3-8.fig',1)
% %uiopen('/Users/alex/Dropbox/iditdiklaalexpaper/Q_9-16.fig',1)
% data = extract_data_from_flowrate_figure(gcf)
% figure
% for i = 1:length(data)
%    subplot(3,3,i);
%     hold on
%        plot(data(i).meas_x{1},data(i).meas_y{1},'b--');
%        plot(data(i).meas_x{2},data(i).meas_y{2},'ro');
%        ylim([-.5 3]);
%        set(gca,'ytick',[0 3]);
%
%    hold off
% end
%


if nargin == 0, help('extract_data_from_flowrate_figure'), return, end

data = struct('meas_x',[],'meas_y',[]);



ha = findobj(hf,'type','axes');

data = repmat(data,1,numel(ha));


for i = 1:numel(ha)
    % each plot is supposedly 2 data sets:
    % measured dots and interpolated values
    ch = get(ha(i),'children');
    if length(ch) ~= 2, error('wrong number of plots'); end
    for j = 1:2
        data(i).meas_x{j} = get(ch(j),'xdata');
        data(i).meas_y{j} = get(ch(j),'ydata');
    end
end




        