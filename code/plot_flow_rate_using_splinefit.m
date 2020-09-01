load('/Users/alex/Dropbox/iditdiklaalexpaper/Measured Data/q_ex2_last.mat')
% t = 0:.2:.8;
% x = [t,1];

f = fieldnames(q);

figure, hold on

for i = 1:length(f)
    tmp = q.(f{i});
    y = [tmp,tmp(1)]; % complete the cycle
    t = 0:0.2:(length(tmp)-1)*0.2;
    x = [t,1]; 
    plot(x,y,'ro');
    pp = splinefit(x,y,3,3,'p');
    fnplt(pp)
end


load('/Users/alex/Dropbox/iditdiklaalexpaper/Measured Data/q_ex3_last.mat')
t = linspace(0,1,8);

figure, hold on
for i = 1:8
    
    plot(t,[q(i,:),q(i,1)],'ro');
    pp = splinefit(t,[q(i,:),q(i,1)],5,5,'p');
    fnplt(pp)
    
end
