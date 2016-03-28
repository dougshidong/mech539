close all;
clear;

set(0,'DefaultTextInterpreter','Latex');

h = figure('Position', [0, 0, 750, 500]);
hold on;
for alfa = 0:2:6
    fname = sprintf('./q3results/q3result%d.mat', alfa);
    load(fname);
    plot(z(1:nbp+1,1),Cpj,'-');
end
grid on
legend('\alpha = 0','\alpha = 2', ...
       '\alpha = 4', '\alpha = 6');
set(gca,'YDir','Reverse');
title('Pressure Distribution for Varying AoA');
xlabel('x')
ylabel('$C_p$')
export_fig ('../report/figures/q3pressure.pdf')
hold off

cllist = zeros(10,1);
cdlist = zeros(10,1);
for alfa = 0:10
    fname = sprintf('./q3results/q3result%d.mat', alfa);
    load(fname);
    cllist(alfa + 1) = Cl;
    cdlist(alfa + 1) = cd;
end

h = figure('Position', [0, 0, 750, 500]);
plot(cdlist,cllist,'-o');
grid on
title('NACA0012 Drag Polar');
xlabel('$C_D$')
ylabel('$C_L$')
export_fig ('../report/figures/q3dp.pdf')