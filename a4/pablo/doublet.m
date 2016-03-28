function res = doublet(z,alfa,holdon);

% This program finds and plots the pressure distribution 
% on an airfoil by representing the surface as a finite 
% number of constant strength doublet panels.
% (Dirichlet Boundary condition)
% Airfoils are taken from the Naca 4 digits library.

% Input data
% z is an array containing the airfoil panels coordinates
% alfa is the angle of attack expressed in degrees

nbp = max(size(z))-1;
chord = 1;
Vzero = 1;

nbpo2 = nbp/2;

% Infinite wake point 

d1n = dist(z,nbp,nbp+1);
d12 = dist(z,2,1);
xP = (z(nbp,1)-z(nbp+1,1))./d1n + (z(2,1)-z(1,1))./d12;
yP = (z(nbp,2)-z(nbp+1,2))./d1n + (z(2,2)-z(1,2))./d12;
dPo = sqrt(xP.^2+yP.^2);
zi(1) = ((z(1,1)+z(nbp+1,1))/2-chord.*2000.*xP./dPo);
zi(2) = ((z(1,2)+z(nbp+1,2))/2-chord.*2000.*yP./dPo);

% Collocation points coordinates

m =(z(1:nbp,:)+z(2:nbp+1,:))./2;

% Upstream velocity potential

alfar = pi.*alfa./180;
Qinf(1) = Vzero*cos(alfar);
Qinf(2) = Vzero*sin(alfar); 
Phiinfi = (Qinf(1).*m(:,1)+Qinf(2)*m(:,2));

% Airfoil velocity potential

miu = Rsolver(nbp,z,m,zi,Phiinfi);

% Velocity

zm = m(2:nbp,:) - m(1:nbp-1,:);
deltaL = sqrt(zm(:,1).^2+zm(:,2).^2);
Qtj = (miu(2:nbp)-miu(1:nbp-1))./deltaL;

Qtj1 = Qtj(1,1)+dist(z,1,2).*(Qtj(1,1)-Qtj(2,1))./dist(z,2,3);
Qtj2 = Qtj(nbp-1,1)+dist(z,nbp,nbp+1).*(Qtj(nbp-1,1)-Qtj(nbp-2,1))./dist(z,nbp-1,nbp);

Vtj(1,1) = (Qtj1-Qtj2)./2;
Vtj(2:nbp,1) = Qtj;
Vtj(nbp+1,1)= - Vtj(1,1);

% Pressure coefficient
normQinf = Qinf(1).^2+Qinf(2).^2;
Cpj = 1- Vtj.^2/normQinf;

% Graphics

cpjmax = max(Cpj); cpjmin = min(Cpj);
if holdon==0
  zplot = z(:,1)-i*(cpjmax-cpjmin)*z(:,2);
  plot(zplot,'k');hold on; 
  axis('ij'); 
  axis([-0.1 1.1 cpjmin cpjmax]); 
end
plot(z(1:nbp+1,1),Cpj,'b'); hold on;

% Lift coefficient

F = 0;
Cm = 0;
cmle =0;

for j=1:nbp
  fj =  (Cpj(j)+Cpj(j+1))*.5*(z(j+1)-z(j))*exp(i*pi/2);
  F = F + fj;
  cmj(j) = real(fj)*imag(z(j+1)+z(j))/2 - imag(fj)*(real(z(j+1)+z(j))/2-chord/4);
  Cm = Cm + cmj(j);
  cmle = cmle + real(fj)*imag(z(j+1)+z(j))/2 - imag(fj)*real(z(j+1)+z(j))/2;
end

Cl = imag(F)*cos(alfar) - real(F)*sin(alfar);
fname = sprintf('q5cresult_%d.mat',1);
save(fname,'z','nbp', 'Cpj', 'Cl');
xcp = -cmle/imag(F);

if Cl<0.001
   Cl = 0;
   Cm = 0;
   xcp = 0;
else
   Cl = floor(10000*Cl)/10000;
   Cm = floor(10000*Cm)/10000;
   xcp = floor(100*xcp)/100;
end;

res(1) = Cl;
res(2) = Cm;
res(3:nbp+3) = Vtj(:,1);
res(nbp+4) = xcp;
res(nbp+5) = cpjmin;
res(nbp+6) = cpjmax;
