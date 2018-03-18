function rescale_flow_rate(ha)
ch = get(ha,'children');

for i = 1:numel(ch)
    y = get(ch(i),'ydata');
    set(ch(i),'ydata',y*.06); % ml/sec to lpm
end
ylabel('Q [lpm]')
set(ha,'ylim',[-.5 3])
set(ha,'ytick',[0 3])
set(ha,'yticklabelmode','auto')


