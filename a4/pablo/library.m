function z = library(nbpo2,select);

filename = ([select,'.DAT']);
load (filename); 
coord = eval(select);

% In the data file : 
% first line = nb points upper S / nb points lower S
% then upper surface starting from LE down to TE (both included)
% then lower surface starting from LE down to TE (both included)

% Reading the data from the file

nbup = coord(1,1); nblo = coord(1,2);
nbp = nbup+nblo-1;

xdat = coord(nbup+1:-1:2,1);        
ydat = coord(nbup+1:-1:2,2);
xdat(nbup+1:nbp) = coord(nbup+3:nbup+nblo+1,1);
ydat(nbup+1:nbp) = coord(nbup+3:nbup+nblo+1,2);

% build the arc length parameter s

s = [0,cumsum(sqrt(diff(xdat').^2+diff(ydat').^2))];

% locate xmin
global PPX;
ppx = spline(s,xdat);
PPX = ppx;
ppy = spline(s,ydat);

% Test wether the data is smooth
% ypl = ppval(ppy,s);
% plot(diff(ypl)./diff(s));
%  pause

xpl = ppval(ppx,linspace(s(1),s(nbp),4*nbpo2));

% find min x (defined as LE)
sle = fminbnd('splf',s(1),s(nbp),[0,1e-6]);

% where is the LE ?
i1 = min(find( s >= sle));
% so bracketed by i1-1 & i1
tt = 0.1*(s(i1)-s(i1-1));

% selection of the rigth LE and creation of splines 
% for upper and lower surfaces with proper LE

if sle - s(i1-1) < tt
%  disp('very close to i1-1');
  ile = i1-1;
  xdatup = xdat(ile:-1:1);  ydatup = ydat(ile:-1:1);
  xdatlo = xdat(ile:nbp);   ydatlo = ydat(ile:nbp);
elseif s(i1) - sle < tt
  ile = i1;
%  disp('very close to i1');
  xdatup = xdat(ile:-1:1);  ydatup = ydat(ile:-1:1);
  xdatlo = xdat(ile:nbp);   ydatlo = ydat(ile:nbp);
else
%  disp('add');
  % add l.e.point to data
  ile = i1;
  xdat2 = xdat(1:ile-1);
  xdat2(ile) = ppval(ppx,sle);
  xdat2(ile+1:nbp+1) = xdat(ile:nbp);
  ydat2 = ydat(1:ile-1);
  ydat2(ile) = ppval(ppy,sle);
  ydat2(ile+1:nbp+1) = ydat(ile:nbp);
  xdatup = xdat2(ile:-1:1);  ydatup = ydat2(ile:-1:1);
  xdatlo = xdat2(ile:nbp+1); ydatlo = ydat2(ile:nbp+1);
end;

% Set the LE at x=0 and y = 0 and scale the airfoil
% so that the TE is at x = 1

chord = xdatup(ile) - xdatup(1);
xle = xdatup(1);
yle = ydatup(1);

xdatup = xdatup-xle;     xdatlo = xdatlo-xle;
xdatup = xdatup/chord;  xdatlo = xdatlo/chord;

ydatup = ydatup-yle;     ydatlo = ydatlo-yle;
ydatup = ydatup/chord;  ydatlo = ydatlo/chord;

% Parametrize the new arc length

sup = [0,cumsum(sqrt(diff(xdatup').^2+diff(ydatup').^2))];
slo = [0,cumsum(sqrt(diff(xdatlo').^2+diff(ydatlo').^2))];

% x distribution   
beta = (0:(pi./nbpo2):pi); 
xc = 0.5*(1-cos(beta)).';   

% Upper surface

% take the square root instead
Vxc = sqrt(xc);

% solve Vx(s) = Vxc 
sVxc = interp1(sqrt(xdatup),sup,Vxc,'spline'); 

% find ycup
ycup = spline(sup,ydatup,sVxc);

% Lower surface

% solve Vx(s) = Vxc 
sVxc = interp1(sqrt(xdatlo),slo,Vxc,'spline'); 

% find yclo
yclo = spline(slo,ydatlo,sVxc);

% Airfoil coordinates

z = [xc(nbpo2+1:-1:1) ycup(nbpo2+1:-1:1) ; xc(2:nbpo2+1) yclo(2:nbpo2+1)];

% treat the blunt TE case

gap = z(1,2)-z(2*nbpo2+1,2);

if gap > 1e-16
  disp('Blunt trailing edge !');
end;

