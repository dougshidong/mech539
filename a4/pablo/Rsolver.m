function miou = Rsolver(nbp,z,m,zi,Phiinfi);

% Airfoil matrix

X1=m(:,1)*ones(1,nbp);
X2=m(:,2)*ones(1,nbp);
A1 = ones(nbp,1)*z(2:nbp+1,1)' - X1;
B1 = ones(nbp,1)*z(1:nbp,1)' - X1;
A2 = ones(nbp,1)*z(2:nbp+1,2)' - X2;
B2 = ones(nbp,1)*z(1:nbp,2)' - X2;

nAtnB = sqrt((A1.^2+A2.^2).*(B1.^2+B2.^2));

det = A1.*B2-A2.*B1;

cosbeta = (A1.*B1+A2.*B2)./(nAtnB);

cosbeta = cosbeta - diag(diag(cosbeta));

Cairfoil = -sign(det)./(2.*pi).*acos(cosbeta);

Cairfoil = real(Cairfoil);

Cairfoil = Cairfoil-diag(diag(Cairfoil))+.5*eye(nbp);

% Wake contribution

A1 = z(1,1) - m(:,1);
B1 = zi(1,1) - m(:,1);
A2 = z(1,2) - m(:,2);
B2 = zi(1,2) - m(:,2);
  
det = A1.*B2 - A2.*B1;

nAtnB = sqrt((A1.^2+A2.^2).*(B1.^2+B2.^2));

cosbeta = (A1.*B1+A2.*B2)./(nAtnB);

Cwake = sign(det)./(2.*pi).*acos(cosbeta);

% Substitution of the kutta condition 

Cairfoil(:,1)   = Cairfoil(:,1)   - Cwake;
Cairfoil(:,nbp) = Cairfoil(:,nbp) + Cwake;

% Computation of the solution

miou = Cairfoil\Phiinfi;
