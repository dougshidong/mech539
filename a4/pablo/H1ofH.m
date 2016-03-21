function y=H1ofH(H);

if H <1.1
   disp('H < 1.1 !  -> H1 = 16');
   y = 16;
else   
  if H <= 1.6       
    y = 3.3 + 0.8234*(H-1.1).^(-1.287);  
  else
    y = 3.3 + 1.5501*(H-0.6778).^(-3.064);
  end;
end;
