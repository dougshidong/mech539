%% a
close all;
clear;
let = 'a';
st = 4;
en = 26;
in = 2;


param = [st:in:en];
i=0;
for alfa = param;
    i=i+1;
    fname = sprintf('./q5results/q5%cresult_%d.mat', let, alfa);
    load(fname);
    cdlist(i) = cd;
end

h = figure('Position', [0, 0, 750, 500]);
hold all
plot(param/100, cdlist,'-o');
grid on
ylabel('$C_D$')
title('Drag vs Thickness');
xlabel('Thickness $\left[ \% \right]$')
export_fig ('../report/figures/q5dthick.pdf')

clear
let = 'b';
st = 2;
en = 8;
in = 1;

param = [st:in:en];
i=0;
for alfa = param;
    i=i+1;
    fname = sprintf('./q5results/q5%cresult_%d.mat', let, alfa);
    load(fname);
    cdlist(i) = cd;
end

h = figure('Position', [0, 0, 750, 500]);
hold all
plot(param/100, cdlist,'-o');
grid on
ylabel('$C_D$')
title('Drag vs Max Camber');
xlabel('Max Camber $\left[ \% \right]$')
export_fig ('../report/figures/q5dmaxcam.pdf')

clear
let = 'c';
st = 2;
en = 8;
in = 1;

param = [st:in:en];
i=0;
for alfa = param;
    i=i+1;
    fname = sprintf('./q5results/q5%cresult_%d.mat', let, alfa);
    load(fname);
    cdlist(i) = cd;
end

h = figure('Position', [0, 0, 750, 500]);
hold all
plot(param/10, cdlist,'-o');
grid on
title('Drag vs Camber Location');
xlabel('Camber Location $\left[ \% \right]$')
export_fig ('../report/figures/q5dcamloc.pdf')