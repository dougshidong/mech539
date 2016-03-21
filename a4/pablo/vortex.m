function res = vortex(za,alfa,holdon);

% This program finds and plots the pressure distribution
% on an airfoil by representing the surface as a finite
% number of linear strength vortex panels.
% (Neumann Boundary condition V.n = 0)
% Airfoils are taken from the Naca 4 digits library.

% Input data
% za is an array containing the airfoil panels coordinates
% alfa is the angle of attack expressed in degrees

nbp = max(size(za))-1; % number of panels
chord = 1;
Vzero = 1;
alfar = pi.*alfa./180;

% Turn z into an array of complex number
z = za(:,1)+i*za(:,2);

% Change z to clockwise
z = z(nbp+1:-1:1);

% Collocation points
m =(z(1:nbp)+z(2:nbp+1))/2;

% Panel angle
th = imag(log(z(2:nbp+1)-z(1:nbp)));

% Free stream normal velocity component

RHSi(1:nbp,1) = cos(alfar)*sin(th(1:nbp))-sin(alfar)*cos(th(1:nbp));

% Influence matrix

% convert collocation pt to panel cs
xzt = m*ones(1,nbp)-ones(nbp,1)*z(1:nbp).';
xt = real(xzt);
zt = imag(xzt);

xz2t = diff(z);
x2t = real(xz2t);
z2t = imag(xz2t);

cth = ones(nbp,1)*cos(th).';
sth = ones(nbp,1)*sin(th).';

X = xt.*cth+zt.*sth;
Z = -xt.*sth+zt.*cth;
X2 = x2t.*cos(th)+z2t.*sin(th);

% compute r1,r2 and th2-th1

mii = m*ones(1,nbp);
zjj = ones(nbp,1)*z(1:nbp).';
zjjp1 = ones(nbp,1)*z(2:nbp+1).';

r1 = abs(zjj-mii);
r2 = abs(zjjp1-mii);

angle = imag(log((zjjp1-mii)./(zjj-mii))) ;
angle = mod(angle-pi,2*pi)-pi;
tmp = X2(:);
X2mat = ones(nbp,1)*tmp';

th2mth1 = angle./(2*pi*X2mat);
RR = log(r2./r1)./(2*pi*X2mat);
u2l = (Z.*RR+X.*th2mth1);
u1l = -(u2l-X2mat.*th2mth1);
cnst = 1/(2*pi);
TMP = 1/(2*pi)-Z.*th2mth1;
w1l = -TMP+(X2mat-X).*RR;
w2l = TMP+X.*RR;

tmp = diag(u1l) + 0.5*( diag(X)-X2  )./X2;
u1l = u1l - diag( tmp );
tmp = diag(u2l) - 0.5*diag(X)./X2;
u2l = u2l - diag(tmp);
tmp = diag(w1l);
w1l = w1l - diag(tmp + 1/(2*pi));
tmp = diag(w2l);
w2l = w2l - diag(tmp - 1/(2*pi));

% Velocity in global cs

ca = ones(nbp,1)*cos(-th)';
sa = ones(nbp,1)*sin(-th)';

u1 = u1l.*ca+w1l.*sa;
u2 = u2l.*ca+w2l.*sa;
w1 = -u1l.*sa+w1l.*ca;
w2 = -u2l.*sa+w2l.*ca;

% Influence matrix coefficient

CA = cos(th)*ones(1,nbp);
SA = sin(th)*ones(1,nbp);
Aij = zeros(nbp,nbp+1);
Bij = Aij;
Aij(:,1:nbp) = -u1.*SA + w1.*CA;
Bij(:,1:nbp) =  u1.*CA + w1.*SA;
Aij(:,2:nbp+1) = Aij(:,2:nbp+1) - u2.*SA + w2.*CA;
Bij(:,2:nbp+1) = Bij(:,2:nbp+1) + u2.*CA + w2.*SA;

% Add a wake panel with a constant-strength vortex

% Infinite wake point 

d1n = dist(za,nbp,nbp+1);
d12 = dist(za,2,1);
xP = (za(nbp,1)-za(nbp+1,1))./d1n + (za(2,1)-za(1,1))./d12;
yP = (za(nbp,2)-za(nbp+1,2))./d1n + (za(2,2)-za(1,2))./d12;
dPo = sqrt(xP.^2+yP.^2);
zi1 = ((za(1,1)+za(nbp+1,1))/2-chord.*2000.*xP./dPo);
zi2 = ((za(1,2)+za(nbp+1,2))/2-chord.*2000.*yP./dPo);

zi = zi1 + i*zi2;
zte = z(1);

d1 = zte*ones(nbp,1) - m;
d2 = zi*ones(nbp,1) - m;

angle = imag(log(d1./d2)) ;
angle = mod(angle-pi,2*pi)-pi;

r1or2 = abs(d1)./abs(d2);

u = 1/(2*pi)*angle;
w = -1/(2*pi)*log(r1or2);

% transfer to global cs

ca = cos(-th);
sa = sin(-th);

ug = u.*ca+w.*sa;
wg = -u.*sa+w.*ca;

% find the tangential component : 
CA = cos(th);
SA = sin(th);

Aw = -ug.*SA + wg.*CA;

Aij = [Aij Aw];

% Kutta condition

Aij(nbp+1,1)     = 1;
Aij(nbp+2,nbp+1) = 1;

RHSi(nbp+1,1) = 0;
RHSi(nbp+2,1) = 0;

% Solve

gamma = Aij\RHSi;

gamma = gamma(1:nbp+1,1);

% Compute velocity

vel =  Bij*gamma;
velocity = vel + cos(alfar)*cos(th)+sin(alfar)*sin(th);

% turn z back to anti-clockwise
z = z(nbp+1:-1:1); 

% turn velocity to anti-clockwise
velocity = -velocity(nbp:-1:1);

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

% compute cp
cp = 1-velocity.^2/Vzero.^2;

% Graphics

cpmax = max(cp); cpmin = min(cp);
if holdon==0
  zplot = real(z)-i*(cpmax-cpmin)*imag(z);
  plot(zplot,'k');hold on; 
  axis('ij'); 
  axis([-0.1 1.1 cpmin cpmax]);  
end;

plot(real(m),cp,'g'); hold on;

% lift coefficient

Fx = 0;
Fy = 0;
cm = 0;
cmle = 0;

for jj=1:nbp
  fxj = -cp(jj)*imag(z(jj+1)-z(jj));
  fyj = cp(jj)*real(z(jj+1)-z(jj));
  Fx = Fx + fxj;
  Fy = Fy + fyj;
  cm = cm + fxj*imag(m(jj)) - fyj*(real(m(jj))-chord/4);
  cmle = cmle + fxj*imag(m(jj)) - fyj*real(m(jj));
end

cl = Fy*cos(alfar) - Fx*sin(alfar);

xcp = -cmle/Fy;

if abs(cl)<0.001
    cl = 0;
    cm = 0;
    xcp = 0;
else
    cl = floor(10000*cl)/10000;
    cm = floor(10000*cm)/10000;
end;

xcp = floor(100*xcp)/100;

res(1) = cl;
res(2) = cm;
res(3:nbp+3) = v;
res(nbp+4) = xcp;
res(nbp+5) = cpmin;
res(nbp+6) = cpmax;
