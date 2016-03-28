close all;

cl_list = zeros(9,1);
cd_list = zeros(9,1);
set(0,'DefaultTextInterpreter','Latex');

h = figure('Position', [0, 0, 750, 500]);
hold on;
i = 0;
for nbp = 20:10:100
    i = i + 1;
    fname = sprintf('./q1results/pressure%d.mat', nbp);
    load(fname);
    plot(z(1:nbp+1,1),Cpj);
    cl_list(i) = Cl;
    cd_list(i) = cd;
end
grid on
legend('20 panels','30 panels','40 panels','50 panels','60 panels',...
       '70 panels','80 panels','90 panels','100 panels');
set(gca,'YDir','Reverse');
title('Pressure Distribution for Various Grid Sizes');
xlabel('x')
ylabel('$C_p$')
export_fig ('../report/figures/q1pressure.pdf')
hold off

h = figure('Position', [0, 0, 1000, 300]);
plot(20:10:100, cl_list, '-o');
title('Lift Convergence');
xlabel('Number of Panels')
ylabel('$C_L$')
export_fig ('../report/figures/q1lift.pdf')



h = figure('Position', [0, 0, 1000, 300]);
plot(20:10:100, cd_list, '-o');
title('Drag Convergence');
xlabel('Number of Panels')
ylabel('$C_D$')
export_fig ('../report/figures/q1drag.pdf')
