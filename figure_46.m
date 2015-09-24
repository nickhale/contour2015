clc, close all

LW = 'LineWidth';
MFC = 'MarkerFaceColor';

grey = .7*[1 1 1];
h = plot(NaN);
blue = get(h, 'color');

lw1 = 1;
lw2 = 2;

r1 = 1;
r2 = 3;

t = linspace(0, 2*pi);
x = linspace(0, 1, 2);

t = linspace(0, 2*pi);
fill3(.5+0*t, r2*cos(t), r2*sin(t), grey, ...
    'facealpha', .5,'linestyle', ':', LW, lw1); hold on
t = linspace(0, 2*pi, 20);
fill3(.51+0*t, r1*cos(t), r1*sin(t), 'w', 'linestyle', ':', LW, lw1); hold on

t = linspace(0, 2*pi);
fill3(0*t, r2*cos(t), r2*sin(t), 'w', ...
    'facealpha', .5,'linestyle', ':', 'LineWidth', 2); hold on
t = linspace(0, 2*pi, 20);
fill3(0*t, r1*cos(t), r1*sin(t), 'w', 'linestyle', ':', LW, lw1); hold on

t = linspace(0, 2*pi);
fill3(.99+0*t, r2*cos(t), r2*sin(t), grey, ...
    'facealpha', .5, 'edgealpha', 0); hold on
t = linspace(0, 2*pi);
fill3(1+0*t, r1*cos(t), r1*sin(t), 'w', 'edgealpha', 0); hold on

t = linspace(0, 2*pi);
plot3(0*t, r1*cos(t), r1*sin(t), ':k', LW, lw1); hold on
plot3(1+0*t, r1*cos(t), r1*sin(t), 'k', LW, lw2);
plot3(1+0*t, r2*cos(t), r2*sin(t), 'k', LW, lw2);
t1 = 1.35;
plot3(x, r1*cos(t1)+0*x, r1*sin(t1)+0*x, ':k', LW, lw1)
plot3(x, r2*cos(t1)+0*x, r2*sin(t1)+0*x, 'k', LW, lw2)
t2 = 4.4;
plot3(x, r1*cos(t2)+0*x, r1*sin(t2)+0*x, ':k', LW, lw1)
plot3(x, r2*cos(t2)+0*x, r2*sin(t2)+0*x, 'k', LW, lw2)
t = linspace(t1, t2);
plot3(0*t, r2*cos(t), r2*sin(t), '-k', LW, lw2);
t = linspace(t2-2*pi,t1, 10);
plot3(0*t, r2*cos(t), r2*sin(t), ':k', LW, lw1);

set(gca, 'view', [25, 18])
view([23, 11])
axis off equal square
xlim([0 1])
ylim([-3 3])
zlim([-3 3])

% print -depsc2 cyl.eps

%%

figure(2)
cla, clc
t = linspace(0, 2*pi);
fill(r2*cos(t), r2*sin(t), grey, 'facealpha', .5, 'LineWidth', 2); hold on
fill(r1*cos(t), r1*sin(t), 'w', 'facealpha', 1, 'LineWidth', 2);
view([0, 90])
rr = linspace(0, 1);
plot(3.5*rr, 0*rr, ':k', LW, 2), hold on
t = pi/4;
plot(2*rr.*cos(t), 2*rr.*sin(t), '-k', LW, 2)
tt = linspace(0, t);
plot(1.5*cos(tt), 1.5*sin(tt), '-k', LW, 2);
z = 1.5*[cos(t), sin(t)];
plot(z(1)+[-.005 0 .2], z(2) + [-.2 0 -.05], '-k', LW, 2);
z = 2*[cos(t), sin(t)];
plot(z(1)+[-.0125 0 -.2], z(2) + [-.2 0 -.01], '-k', LW, 2);
axis equal square off
axis(1.05*[-3 3 -3 3])

% print -depsc2 cyl2.eps