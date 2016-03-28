close all;
clear;
set(0,'DefaultTextInterpreter','Latex');

cllist1 = zeros(13,1);
cdlist1 = zeros(13,1);
cllist2 = zeros(13,1);
cdlist2 = zeros(13,1);
cllist3 = zeros(13,1);
cdlist3 = zeros(13,1);
for alfa = -2:10
    fname = sprintf('./q4results/q4result1_%d.mat', alfa);
    load(fname);
    cllist1(alfa + 3) = Cl;
    cdlist1(alfa + 3) = cd;
    fname = sprintf('./q4results/q4result2_%d.mat', alfa);
    load(fname);
    cllist2(alfa + 3) = Cl;
    cdlist2(alfa + 3) = cd;
    fname = sprintf('./q4results/q4result3_%d.mat', alfa);
    load(fname);
    cllist3(alfa + 3) = Cl;
    cdlist3(alfa + 3) = cd;
    
end
h = figure('Position', [0, 0, 750, 500]);
hold all
plot(cdlist1,cllist1,'-o');
plot(cdlist2,cllist2,'-s');
plot(cdlist3,cllist3,'-v');
grid on
title('Drag Polar for Various Airfoils');
legend('NACA0012','GA(W)-1','DAE31','Location','best');
xlabel('$C_D$')
ylabel('$C_L$')
export_fig ('../report/figures/q4dp.pdf')


h = figure('Position', [0, 0, 750, 500]);
hold on;
for alfa = 1:3
    fname = sprintf('./q4results/q4result%d_4.mat', alfa);
    load(fname);
    plot(z(1:nbp+1,1),Cpj,'-');
end
grid on
legend('NACA0012','GA(W)-1','DAE31');
set(gca,'YDir','Reverse');
title('Pressure Distribution for $\alpha = 4^\circ$');
xlabel('x')
ylabel('$C_p$')
export_fig ('../report/figures/q4pressure.pdf')
hold off

co = [            0   4.4700e-01   7.4100e-01
   8.5000e-01   3.2500e-01   9.8000e-02
   9.2900e-01   6.9400e-01   1.2500e-01
   4.9400e-01   1.8400e-01   5.5600e-01
   4.6600e-01   6.7400e-01   1.8800e-01
   3.0100e-01   7.4500e-01   9.3300e-01
   6.3500e-01   7.8000e-02   1.8400e-01];
h = figure('Position', [0, 0, 750, 500]);
hold on;
for alfa = 1:3
    fname = sprintf('./q4results/q4result%d_4.mat', alfa);
    load(fname);
    plot(supper(1:size(cfupper,2)),cfupper,'-','Color', co(alfa,:));
end
grid on
legend('NACA0012','GA(W)-1','DAE31');
title('Skin Friction Coefficient on Upper Surface');
xlabel('Arc Length')
ylabel('$C_f$')

% Label Transition Point
sq = plot(NaN, NaN, 'sk');
fname = sprintf('./q4results/q4result1_4.mat');
load(fname);
plot(supper(37),cfupper(37),'s','Color', co(1,:));
fname = sprintf('./q4results/q4result2_4.mat');
load(fname);
plot(supper(46),cfupper(46),'s','Color', co(2,:));
fname = sprintf('./q4results/q4result3_4.mat');
load(fname);
plot(supper(89),cfupper(89),'s','Color', co(3,:));

legend('NACA0012','GA(W)-1','DAE31','Transition');

export_fig ('../report/figures/q4uppercf.pdf')
hold off


h = figure('Position', [0, 0, 750, 500]);
hold on;
for alfa = 1:3
    fname = sprintf('./q4results/q4result%d_4.mat', alfa);
    load(fname);
    plot(slower(1:size(cflower,2)),cflower,'-','Color', co(alfa,:));
end
%plot(slower(1:size(cflower,2)),cflower,'s','Color', co(2,:));
grid on
title('Skin Friction Coefficient on Lower Surface');
xlabel('Arc Length')
ylabel('$C_f$')
sq = plot(NaN, NaN, 'sk');
fname = sprintf('./q4results/q4result1_4.mat');
load(fname);
plot(slower(164),cflower(164),'s','Color', co(1,:));
fname = sprintf('./q4results/q4result2_4.mat');
load(fname);
plot(slower(118),cflower(118),'s','Color', co(2,:));
legend('NACA0012','GA(W)-1','DAE31','Transition');

export_fig ('../report/figures/q4lowercf.pdf')
hold off
















