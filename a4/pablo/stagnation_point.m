function isp = stagnation_point(Qtj);

isp = 2; 

if Qtj(isp)>0
   while Qtj(isp)>0  
       isp = isp+1;
   end;
elseif Qtj(isp)<0
   while Qtj(isp)<0  
       isp = isp+1;
   end;
end;

ispu = isp -1;
ispl = isp;

if abs(Qtj(ispu)) < abs(Qtj(ispl)) 
   isp = ispu;
else
   isp = ispl;
end; 
