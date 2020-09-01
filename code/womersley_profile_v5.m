function F = womersley_profile_v5(x,xdata,T,nu,a,KQ,KQ0,nf,rho)


t =x;

y = xdata;

ny = length(y)-1;

temp2 = 1;


% nt = length(t);




w0 = 2*pi/T;			    % fundamental radian frequency
alpha0 = a*sqrt(w0/nu);	    % Womersley parameter, fundamental frequency
w = w0*[1:nf];			    % array of radian frequencies
alpha = a*sqrt(w/nu);	    % array of Womersley parameters
[K,Ktau] = deal(zeros(size(KQ)));					% initialize pressure gradient array
Q = KQ0*temp2;				% initialize flow rate array with DC flow component
Pp = -8*nu*Q/(pi*a^4);	% initialize pressure gradient array with DC component
Tau_steady = -a*Pp/2;	% initialize wall shear stress array with DC component
Tau = Tau_steady;
j = sqrt(-1);
j32 = j^1.5;



u = (2*KQ0/(pi*a^2))*((1 - y.*y)*temp2);

for n = 1:nf		% Sum over all frequencies
   K(n) = KQ(n)*j*w(n)/(pi*a^2*(1 - 2*besselj(1,j32*alpha(n))/(j32*alpha(n)*besselj(0,j32*alpha(n)))));
   Ktau(n) = rho*K(n)*a*besselj(1,j32*alpha(n))/(j32*alpha(n)*besselj(0,j32*alpha(n)));
   
 for k = 1:length(t)		% Calculate for each point in time
      Q(k) = Q(k)+real(KQ(n)*exp(j*w(n)*t(k)));
      Pp(k) = Pp(k)+real(K(n)*exp(j*w(n)*t(k)));
      Tau(k) = Tau(k)+real(Ktau(n)*exp(j*w(n)*t(k)));
   for m = 1:ny		% Calculate for each spatial location y
         u(m,k) = u(m,k) + real((K(n)*a^2/(nu*alpha(n)^2*1i))*(1-besselj(0,j32*alpha(n)*y(m))/besselj(0,j32*alpha(n)))*exp(1i*w(n)*t(k)));
      end	% y loop
   end	% t loop
end	% f loop



F = u(:,k);


% keyboard
