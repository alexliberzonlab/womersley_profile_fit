%% Steady Flow,  Rigid Tube
clc, clear
% parameters
dP=10;
mu=0.01;
ro=1;
L=63;
R=0.5;

% Define spatial mesh
r=-R:0.01:R;

% Calculate  velocity
u=-dP/4/mu/L.*(r.^2-R^2);


% Plot
plot(r,u)
title ('Steady Flow,  Rigid Tube')
xlabel ('Radius (r)')
ylabel ('Velocity (u_z)')



%%  pulsatile Flow,  Rigid Tube

clear, clc, clf

% parameters definitions
mu=0.01; % fluid viscosity
ro=1;       % fluid density
ks=-10/62;      % ks=dp/dx=-Dp/L
R=0.5;           % tube radios
BPM=7;    % frequency [BPM]
T=60/BPM;  w=2*pi*BPM/60;  % cycle piriod


% Define spatial and time mesh
r=linspace(-R,R,50);
N=24;         % N number of time profiles in plot
t=linspace(0,T,N); 

% periodic pressure
% coeficients
alfa=sqrt(ro*w/mu)*R;            % Womersley Number
Lambda=((1i-1)/sqrt(2))*alfa;    % complex frequency parameter 4.5.5 p. 76
xi=Lambda/R.*r;                  % complex radial parameter

% calculating U

Jx=besselj(0,xi);
JL=besselj(0,Lambda);
us=-ks*R^2/4/mu;
Ur=us*(-4./Lambda.^2.*(1-Jx./JL));

ut=exp(1i.*w.*t);
UU=real(ut'*Ur);

% plot  settings

% mov = VideoWriter('rigid.avi');

figure, hold on
% plotting
for k=1:length(t)
    
    plot(r',UU(k,:))
        
CM=jet;
CM=CM(1:2:64,:);
set(0,'DefaultAxesColorOrder',CM)
title (['Rigid Tube , Pulsatile Flow,   \alpha=',num2str(alfa)])
xlabel ('Radius (r)')
ylabel ('Velocity (u_x)')

if BPM>1
    C=1.2/BPM;
else
    C=1;
end

ylim ([-us  us]*C)  %limits - the steady flow peak velocity
hold off

 hold all   % remark if want seperate profiles in movie

    % V(k) = getframe(figure(1));
    % mov = writeVideo(mov,V(k));
end
% mov=close(mov);
hold off
saveas(gcf,sprintf('../results/%d.png',k));



%% Qrel (Q(t)max/Qsteady)
clear, clc
 
alfa=linspace(0.0001,10);        %alfa
Lambda=1i^1.5.*alfa;                         %Lambda
J1=besselj(1,Lambda);
J0=besselj(0,Lambda);
g=2.*J1./Lambda./J0;
Qrel=real(-8./Lambda.^2.*(1-g));
 
plot (alfa,Qrel)
title ('Q_{puls}/Q_{steady}  , Rigid Tube')
xlabel ('alfa')
ylabel ('Qrel')
saveas(gcf,'../results/alpha_Qrel.png');


%%  pulsatile Flow,  elastic Tube


clear, clc%, clf

% parameters definitions
mu=0.01; % fluid viscosity
ro=1;       % fluid density
ni=0.75;      %Tube Poison's coeficient
% E=9E3;     %Tube Module of Elasticity - what are the units? 
E = 12E12; % 12E6; % Tygon is 0.005 - 0.008 GPa in paper we state 12Mpa
h=0.1;           %Tube's thickness
ks=-10/62;      % ks=dp/dx=-Dp/L
R=.5;           % tube radius, Alex, what are the units? cm? 
BPM=7;    % frequency [BPM]
T=60/BPM;  w=2*pi*BPM/60;  % cycle piriod
L=62;           %tube's length

% Define spatial and time mesh
r=linspace(-R,R,100);
t=linspace(0,T,24);

x=20;

% periodic pressure
% coeficients
c0=E*h/2/R/ro;                         %wave's speed
alfa=sqrt(ro*w/mu)*R;           % Womersley Number
Lambda=(1i^1.5)*alfa;     % complex frequency parameter
xi=Lambda/R.*r;                            % complex radial parameter


% wave parameters
c=c0*(1+1./-exp(alfa)); % <-- why this appears twice? Alex
Jx=besselj(0,xi);
J0=besselj(0,Lambda);
J1=besselj(1,Lambda);
g=2.*J1./Lambda./J0;
%c=c0*(1+1.1./-exp(alfa)+(0.3*sin(alfa).*exp(-alfa)*1i));  % <-- why this appears twice? Alex                                           %no dumping
z=2*c0^2./(1-ni^2)./c^2;
G=(2+z.*(2*ni-1))./z./(2*ni-g);


% calculating U (axial velocity)
us=-ks*R^2/4/mu;
Urx= us*(-4./Lambda.^2.*(1-G.*Jx./J0));
ut=exp(1i.*w.*(t-x./c));
UU=real(ut'*Urx);

% calculating V (radial velocity)
us=-ks*R^2/4/mu;
Vrx= us*(2.*R.*w./1i./Lambda.^2./c*(r./R-G.*2.*J1./J0./Lambda));
vt=exp(1i.*w.*(t-x./c));
VV=real(vt'*Vrx);

% plot  settings

% mov = avifile('elastic.avi');

figure, hold on
% plotting
for k=1:length(t)
    
    plot(r',UU(k,:))
    % plot(r',VV(k,:))
        
CM=jet;
CM=CM(1:2:64,:);
set(0,'DefaultAxesColorOrder',CM)
title (['Elastic Tube, E=', num2str(E),' x=', num2str(x),' , Pulsatile Flow,   \alpha=',num2str(alfa)])
xlabel ('Radius (r)')
% ylabel ('Velocity (u_x)')
ylabel ('Velocity (v_x)')
if BPM>1
    C=1.9/BPM;
else
    C=1;
end
C=max(max(UU))
% ylim ([-us  us]*C)  %limits - the steady flow peak velocity

% hold off

hold all   % remark if want seperate profiles in movie

    % V(k) = getframe(figure(1));
    % mov = addframe(mov,V(k));
end
% mov=close(mov);
hold off
saveas(gcf,sprintf('../results/elastic%d.png',k));




%% G(lambda)
clear, clc, clf
 
alfa=linspace(0.0001,10);        %alfa
Lambda=i^1.5.*alfa;                         %Lambda
ni=0.8;
J0=besselj(0,Lambda);
J1=besselj(1,Lambda);
g=2.*J1./Lambda./J0;
c0=3000;                                             %no dumping
c=c0*(1+1.1./-exp(alfa)+(0.3*sin(alfa).*exp(-alfa)*i));
z=2*c0^2./(1-ni^2)./c.^2;
G=(2+z.*(2*ni-1))./z./(2*ni-g);

realG=real(G);
imgG=imag(G);
 
plot (alfa,realG, '-r', alfa, imgG, ':b')
title ('G(alpha)')
xlabel ('alfa')
ylabel ('G')
ylim ([-0.7 1.7])
saveas(gcf,'../results/alfa_realG.png');

%% c(alpha)
clear, clc, clf
 
alfa=linspace(0.0001,10);        %alfa

 c0=3000;   
C=c0*(1+1./-exp(alfa)+(0.7*sin(alfa).*exp(-alfa)*i));
plot (alfa,real(C)/c0,  '-r', alfa, imag(C)/c0, ':b')
title ('C(alpha)')
xlabel ('alfa')
ylabel ('C')
saveas(gcf,'../results/alpha_real_C.png');

% ylim ([-0.7 1.7])