clc, close all

LW = 'LineWidth';
MFC = 'MarkerFaceColor';
b = .2;

grey = .7*[1 1 1];
h = plot(NaN);
blue = get(h, 'color');

t = linspace(0,1, 2);
plot3(t, 0*t, 0*t, 'k', LW, 3); hold on
plot3(t, 0*t+b, 0*t, ':k', LW, 3);
plot3(t, 0*t, 0*t+b, 'k', LW, 3);
plot3(t, 0*t+b, 0*t+b, 'k', LW, 3);

plot3(0*t, t*b, 0*t, ':k', LW, 3)
plot3(0*t, t*b, 0*t+b, 'k', LW, 3)
plot3(0*t, 0*t, b*t, 'k', LW, 3)
plot3(0*t, 0*t+b, b*t, ':k', LW, 3)

plot3(1+0*t, t*b, 0*t, 'k', LW, 3)
plot3(1+0*t, t*b, 0*t+b, 'k', LW, 3)
plot3(1+0*t, 0*t, b*t, 'k', LW, 3)
plot3(1+0*t, 0*t+b, b*t, 'k', LW, 3)

plot3(.5*[1 1 1 1 1], [0 b b 0 0], [0 0 b b 0], ':k', LW, 1)
b2 = b/2;
fill3(.5*[1 1 1 1 1], [b2 b b b2 b2], [b2 b2 b b b2], grey, 'facealpha', .5)
h = fill3([1 1 1 1 1], [0 b b 0 0], [0 0 b b 0], grey, 'facealpha', .5);

set(gca, 'view', [25, 18])
set(gca, 'xtick', [0 1], 'ytick', [0 b], 'ztick', [0 b])
set(gca, 'xticklabel', {}, 'yticklabel', {}, 'zticklabel', {})

axis off equal
