function z  = naca4(x,Extra);

% Decoding arguments
eps = x(1)./100;
p = x(2)./10;
to = x(3)./100;

nbpo2 = Extra(1);
c = 1;

% x distribution   
beta = (0:(pi./nbpo2):pi).'; 
xc = c.*(1-.5.*(1-cos(beta)));   

% Thickness distribution = f(x)
thdis=5.*to.*c.*(0.2969.*sqrt(xc./c)-0.126.*xc./c-0.3537.*(xc./c).^2 +0.2843.*(xc./c).^3-0.1015.*(xc./c).^4);

% Camberline = f(x)
if p~=0 & eps~=0 
 I1=find(xc(1:nbpo2+1)<=p*c);
 I2=find(xc(1:nbpo2+1)>p*c);
 camberline(I1,1) = (eps.*xc(I1))./(p.^2).*(2.*p-xc(I1)./c);
 camberline(I2,1) = (eps.*(c-xc(I2)))./(1-p).^2.*(1+xc(I2)./c-2.*p);
end;

% Airfoil = camberline and thickness
if p==0 | eps==0
    xupper =  xc;
    yupper =  thdis;
    xlower =  xc;
    ylower = -thdis;
else
    theta(I1,1) = atan(2.*eps./p.*(-xc(I1)./(c.*p)+1));
    theta(I2,1) = atan(2.*eps./(1-p.^2).*(p-(xc(I2)./c)));
    xupper = xc         - thdis.*sin(theta);
    yupper = camberline + thdis.*cos(theta);
    xlower = xc         + thdis.*sin(theta);
    ylower = camberline - thdis.*cos(theta); 
end;

z = [xupper yupper ; xlower(nbpo2:-1:1,1) ylower(nbpo2:-1:1,1)];
