function res = source(z,holdon);

% This program finds and plots the pressure distribution 
% on an airfoil by representing the surface as a finite 
% number of constant strength source panels.
% (Neumann Boundary condition V.n = 0)
% Airfoils are taken from the Naca 4 digits library.

% Input data
% z is an array containing the airfoil panels coordinates

chord = 1; 
Vzero = 1;

nbp = max(size(z))-1;

% Turn z into an array of complex number
z = z(:,1)+i*z(:,2);

% Collocation points
m =(z(1:nbp)+z(2:nbp+1))/2;

% Panel angle
alfa(1:nbp) = imag(log(z(1:nbp)-z(2:nbp+1)));

% Influence matrix

% normal and tangent vectors to the panel
vecpanel = -diff(z);
tangent = vecpanel./abs(vecpanel);
normal = tangent*exp(i*pi/2);

% rhs term
RHSi = -Vzero*real(normal);

% Matrix 

mii = m*ones(1,nbp);
zjj = ones(nbp,1)*z(1:nbp).';
zjjp1 = ones(nbp,1)*z(2:nbp+1).';

up = 1/(2*pi)*log(abs(zjjp1-mii)./abs(zjj-mii));
angle = (imag(log(zjj-mii))-imag(log(zjjp1-mii))) ;
angle = mod(angle-pi,2*pi)-pi;
wp = angle/(2*pi);

% diagonal terms
up = up-diag(diag(up));
wp = wp-diag(diag(wp))+.5*eye(nbp);

% Velocity in global cs
ca = ones(nbp,1)*cos(-alfa);
sa = ones(nbp,1)*sin(-alfa);

u = up.*ca+wp.*sa;
w = -up.*sa+wp.*ca;

normal  = normal*ones(1,nbp);
tangent = tangent*ones(1,nbp);

% matrix coefficient = normal velocity
Aij = u.*real(normal)+w.*imag(normal);

% for later cp computation : tangential velocity  
Bij = u.*real(tangent)+w.*imag(tangent);

% Solve

sigma = Aij\RHSi;

% Compute velocity and cp

velocity = Bij*sigma;

velocity = velocity'+cos(alfa);

% Velocity at nodes
Qtj = (velocity(1:nbp-1) + velocity(2:nbp))./2;

dz12 = abs(z(1)-z(2));
dz23 = abs(z(3)-z(2));
dznbp1 = abs(z(nbp)-z(1));
dznbpnbpm1 = abs(z(nbp)-z(nbp-1));

Qtj1 = Qtj(1)+dz12.*(Qtj(1)-Qtj(2))./dz23;
Qtj2 = Qtj(nbp-1)+dznbp1.*(Qtj(nbp-1)-Qtj(nbp-2))./dznbpnbpm1;

v(1) = (Qtj1-Qtj2)./2;
v(2:nbp) = Qtj;
v(nbp+1)= - v(1);

cp = 1-velocity.^2/Vzero.^2;

% Graphics

cpmax = max(cp); cpmin = min(cp);
if holdon==0
  zplot = real(z)-i*(cpmax-cpmin)*imag(z);
  plot(zplot,'k');hold on; 
  axis('ij'); 
  axis([-0.1 1.1 cpmin cpmax]); 
end;
plot(real(m),cp,'r'); hold on;

% lift coefficient

res(1) = 0;
res(2) = 0;
res(3:nbp+3) = -v;
res(nbp+4) = 0;
res(nbp+5) = cpmin;
res(nbp+6) = cpmax;
