function res = solvebl(alfa,Re,z,n,ue,plotcp,side);

nbp2 = 200;  % bl discretization

% Arc length
ss(1,1) = 0;
for ii=2:n
  ss(ii,1) = ss(ii-1,1)+dist(z,ii,ii-1);
end;

sTE = ss(n,1);

% compute the boundary layer up to xe

Coff = 0.98;

nm = floor(0.5*n);
se = spline(z(nm:n,1),ss(nm:n),Coff); 

s = 0:se/nbp2:se;

% Velocity

spues = spline(ss,ue);
ues = ppval(spues,s);

% detect if there is some ue < 0
% which means a prb with the spline interpolation

in = find(ues <0);

if ~isempty(in)

  % Add some more points near the LE
  ue2(1) = ue(1);
  ue2(2) = 0.5*(ue(2)+ue(1));
  ue2(3) = ue(2);
  ue2(4) = 0.5*(ue(3)+ue(2));
  ue2(5) = ue(3);
  ue2(6) = 0.5*(ue(4)+ue(3));
  ue2(7:n+3) = ue(4:n);

  ss2(1) = ss(1);
  ss2(2) = 0.5*(ss(2)+ss(1));
  ss2(3) = ss(2);
  ss2(4) = 0.5*(ss(3)+ss(2));
  ss2(5) = ss(3);
  ss2(6) = 0.5*(ss(4)+ss(3));
  ss2(7:n+3) = ss(4:n);

  % re-spline
  spues = spline(ss2,ue2);
  ues = ppval(spues,s);

end;

% x coordinate

spx = spline(ss,z(:,1));

%plot(s,ues,'k');
%hold on;
%plot(ss,ue,'o');
%pause

ue = ues;

n = nbp2+1;

% x coordinate

spx = spline(ss,z(:,1));

% velocity gradient at nodes

v1= ue(1);  v2 = ue(2);  v3 = ue(3);
x1= s(1); x2 = s(2); x3 = s(3);
gamma = 1./(x3-x2)*((v3-v1)./(x3-x1)-(v2-v1)./(x2-x1));
fac = (v2-v1)./(x2-x1);
dueds(1) = gamma*(x1-x2) + fac;

if dueds(1) < 0
  dueds(1) = (v2-v1)/(x2-x1);
end;

v1=ue(1:n-2); v2=ue(2:n-1); v3=ue(3:n);
x1=s(1:n-2);x2=s(2:n-1);x3=s(3:n);
gamma = 1./(x3-x2).*((v3-v1)./(x3-x1)-(v2-v1)./(x2-x1));
fac = (v2-v1)./(x2-x1);
dueds(2:n-1) = gamma.*(x2-x1) + fac;

v1= ue(n-2);  v2 = ue(n-1);  v3 = ue(n);
x1= s(n-2); x2 = s(n-1); x3 = s(n);
gamma = 1./(x3-x2)*((v3-v1)./(x3-x1)-(v2-v1)./(x2-x1));
fac = (v2-v1)./(x2-x1);
dueds(n) = gamma*(2*x3-x1-x2) + fac;

%--------Laminar boundary layer

lsep = 0; trans=0; endofsurf=0;

theta(1) = sqrt(0.075/(Re*dueds(1)));
i = 1;

while lsep ==0 & trans ==0 & endofsurf ==0

  lambda = theta(i).^2*dueds(i)*Re;

  % test for laminar separation
  if lambda < -0.09 
    lsep = 1;
    itrans = i;
    break; 
  end;

  H(i) = fH(lambda);
  L = fL(lambda);

  cf(i) = 2*L./(Re*theta(i));
  if i>1, cf(i) = cf(i)./ue(i); end;
  i = i+1;

  % test for end of surface
  if i> n endofsurf = 1; itrans = n; break; end;  

  K = 0.45/Re;
  xm = (s(i)+s(i-1))/2;
  dx = (s(i)-s(i-1));
  coeff = sqrt(3/5);

  f1 = ppval(spues,xm-coeff*dx/2); f1 = f1^5;
  f2 = ppval(spues,xm);            f2 = f2^5;
  f3 = ppval(spues,xm+coeff*dx/2); f3 = f3^5;

  dth2ue6 = K*dx/18*(5*f1+8*f2+5*f3);
  theta(i) = sqrt((theta(i-1).^2*ue(i-1).^6 + dth2ue6)./ue(i).^6);

  % test for transition
  rex = Re*s(i)*ue(i);
  ret = Re*theta(i)*ue(i);
  retmax = 1.174*(rex^0.46+22400*rex^(-0.54));
  if ret>retmax 
    trans = 1; 
    itrans = i;
  end;

end;

%-------- Transition 

transorlamsep = 0;
transloc = 1;
tsep = 0;

if itrans < n

  if trans == 1

     uei = ue(i); thi = theta(i); si = s(i); duedsi = dueds(i);
     ueim1 = ue(i-1); thim1 = theta(i-1); sim1 = s(i-1); duedsim1 = dueds(i-1);

     % Find f(x) at i and i-1

     fxi = ret - retmax;  % already computed 
 
     rex = Re*sim1*ueim1;
     ret = Re*thim1*ueim1;
     retmax = 1.174*(rex^0.46+22400*rex^(-0.54));  
     fxim1 = ret - retmax;

     % Fit a linear function and find the root

     st = sim1 - fxim1/((fxi-fxim1)/(si-sim1));

     transorlamsep = 1;
     transloc = 100*ppval(spx,st);

     % Find the value of theta, and H at st using thwaites

     uet = ppval(spues,st);

     v1=ueim1; v2=uet; v3=uei;
     x1=sim1;  x2=st;  x3=si;
     gamma = 1./(x3-x2).*((v3-v1)./(x3-x1)-(v2-v1)./(x2-x1));
     fac = (v2-v1)./(x2-x1);
     duedst = gamma.*(x2-x1) + fac;

     xm = (st+sim1)/2;
     dx = (st-sim1);

     f1 = ppval(spues,xm-coeff*dx/2); f1 = f1^5;
     f2 = ppval(spues,xm);            f2 = f2^5;
     f3 = ppval(spues,xm+coeff*dx/2); f3 = f3^5;

     dth2ue6 = K*dx/18*(5*f1+8*f2+5*f3);
     thetat  = sqrt((thim1.^2*ueim1.^6 + dth2ue6)./uei.^6);

     lambdat = thetat.^2*duedst*Re;
     Ht = fH(lambdat);

     if Ht < 1.1
        Ht = 1.2;
     end;

     if Ht > 2  % to avoid turbulent separation just after transition
        Ht = 2;
     end;

     % Find the value of theta, and H at i using head

     y(1) = thetat;
     y(2) = H1ofH(Ht);
     dx = s(i) - st;

     y = runge(dx,y,Re,uet,duedst,ue(i),dueds(i));

     theta(i) = y(1);
     H(i) = HofH1(y(2));
     rtheta = Re*ue(i)*theta(i);
     cf(i) = cfturb(rtheta,H(i));     

elseif lsep == 1

     uei = ue(i); thi = theta(i); si = s(i); duedsi = dueds(i);
     ueim1 = ue(i-1); thim1 = theta(i-1); sim1 = s(i-1); duedsim1 = dueds(i-1);

     % Find f(x) at i and i-1

     fxi = thi.^2*duedsi*Re+0.09;
     fxim1 = thim1.^2*duedsim1*Re+0.09;
 
     % fit a linear function and find the root

     st = sim1 - fxim1/((fxi-fxim1)/(si-sim1));
 
     transorlamsep = 2;
     transloc = 100*ppval(spx,st);

     % Find the value of theta, and H at st using thwaites

     uet = ppval(spues,st);

     v1=ueim1; v2=uet; v3=uei;
     x1=sim1;  x2=st;  x3=si;
     gamma = 1./(x3-x2).*((v3-v1)./(x3-x1)-(v2-v1)./(x2-x1));
     fac = (v2-v1)./(x2-x1);
     duedst = gamma.*(x2-x1) + fac;

     xm = (st+sim1)/2;
     dx = (st-sim1);

     f1 = ppval(spues,xm-coeff*dx/2); f1 = f1^5;
     f2 = ppval(spues,xm);            f2 = f2^5;
     f3 = ppval(spues,xm+coeff*dx/2); f3 = f3^5;

     dth2ue6 = K*dx/18*(5*f1+8*f2+5*f3);
     thetat  = sqrt((thim1.^2*ueim1.^6 + dth2ue6)./uei.^6);

     lambdat = thetat.^2*duedst*Re;
     Ht = fH(lambdat);

     if Ht < 1.1
        Ht = 1.2;
     end;

     if Ht > 2  % to avoid turbulent separation just after transition
        Ht = 2;
     end;

     % Find the value of theta, and H at i using head

     y(1) = thetat;
     y(2) = H1ofH(Ht);
     dx = s(i) - st;

     y = runge(dx,y,Re,uet,duedst,ue(i),dueds(i));

     theta(i) = y(1);
     H(i) = HofH1(y(2));
     rtheta = Re*ue(i)*theta(i);
     cf(i) = cfturb(rtheta,H(i));     
     
end;

%--------TURBULENT BL

  tsep = 0;
  
  i = i+1;
  
  while endofsurf == 0 & tsep ==0;

    y = runge(s(i)-s(i-1),y,Re,ue(i-1),dueds(i-1),ue(i),dueds(i));
 
    theta(i) = y(1);
    H(i) = HofH1(y(2));

    if H(i) == 3 % which is actually a flag 
        tsep = 100*ppval(spx,s(i));
        i = i-1;
    end;

    rtheta = Re*ue(i)*theta(i);
    cf(i) = cfturb(rtheta,H(i));

    i = i+1;

    if i>n endofsurf = 1; end;

  end;

end;

if plotcp == 1

    deltas = H(1:i-1).*theta(1:i-1);

    if side ==1
  
      figure;plot(s(1:i-1),theta(1:i-1));grid;
      h = title('Upper Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Momentum Thickness Theta');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),deltas);grid;
      h = title('Upper Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Displacement Thickness delta star');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),H(1:i-1));grid;
      h = title('Upper Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Shape Factor H');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),cf(1:i-1));grid;
      h = title('Upper Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Skin Friction Coefficient Cf');set(h,'Fontsize',[14]);

    elseif side == 2

      figure;plot(s(1:i-1),theta(1:i-1));grid;
      h = title('Lower Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Momentum Thickness Theta');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),deltas);grid;
      h = title('Lower Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Displacement Thickness delta star');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),H(1:i-1));grid;
      h = title('Lower Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Shape Factor H');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),cf(1:i-1));grid;
      h = title('Lower Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Skin Friction Coefficient Cf');set(h,'Fontsize',[14]);
      fname = sprintf('q4result3_%d.mat', alfa);
      cflower=cf;
      slower=s;
      save(fname,'slower','cflower','-append');
    end;

end;

res(1) = theta(i-1);
res(2) = H(i-1);
res(3) = ue(i-1);
res(4) = transorlamsep;  % =0 is fully laminar, =1 if transition, =2 if laminar separation
res(5) = floor(100*transloc)/100;
res(6) = floor(100*tsep)/100; % if =0, means that there is no turbulent separation
