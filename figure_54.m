clc, close all, clear all

LW = 'LineWidth'; lw = 3;
MFC = 'MarkerFaceColor';
d = pi/6;
ell = pi;
grey =.9*[1 1 1];

w = linspace(-10,10,1000);
wp = w + 1i*pi/2;
wm = w - 1i*pi/2;

k = 0;
N = 7;
for s = linspace(-pi/2, pi/2, N);
    k = k + 1;
    ww{k} = w + 1i*s + 1i*eps;
end
ww([1, (N+1)/2, end]) = [];
N = N - 3;

cols = get(0, 'DefaultAxesColorOrder');
blue = cols(1,:);
red = cols(2,:);

c1 = blue;
c2 = red;


figure(1)
c = 10;
fill([-10 10 10 -10 -10], pi/2*[-1 -1 1 1 -1], grey, ...
    'edgecolor', 'none'); hold on
for k = 1:N
    plot(ww{k}, '--k', LW, 1); hold on
end
h1 = plot(wp, LW, lw, 'color', blue); hold on
h2 = plot(wm, LW, lw, 'color', red); 
plot(w+1i*eps, '--k', LW, lw)
plot(0, pi/2, 'o', 'color', c1, MFC, c1);
plot(0, -pi/2, 'o', 'color', c2, MFC, 'w');
trim = 5;
% h2 = plot([0], [-2], '.w', 'markerfacecolor', 'none', 'markersize', eps);
h2 = plot([0], [-trim], '.w', 'markerfacecolor', 'none', 'markersize', eps);
h2 = plot([0], [trim], '.w', 'markerfacecolor', 'none', 'markersize', eps);
axis equal
axis(trim*[-1 1 -1 1])
hold off
set(gca, 'xtick', [])
set(gca, 'ytick', [])
axis off

print -depsc2 ~/Dropbox/work/2015/contour/paper/figures/sectorala

%%

figure(2)
v = @(w) .5*(1+1i*sinh(w));
fill(100*[-1 1 1 -1 -1], 100*[-1 -1 1 1 -1], grey, ...
    'edgecolor', 'none'); hold on
for k = 1:N
    plot(v(ww{k}), '--k', LW, 1); hold on
end
plot(v(wp), LW, lw, 'color', blue); hold on
plot(v(wm), LW, lw, 'color', red);
plot(v(w), '--k', LW, lw)
plot(0, 0, 'o', 'color', c1, MFC, c1);
plot(1, 0, 'o', 'color', c2, MFC, 'w');
axis equal
axis([-2 2 -2 2])
set(gca, 'xtick', [0,  1], 'xticklabel', {'',''})
% set(gca, 'xtick', [])
set(gca, 'ytick', [])
hold off
set(gca, 'box', 'off')
axis off
print -depsc2 ~/Dropbox/work/2015/contour/paper/figures/sectoralb

%%
figure(3)
z = @(w) (ell^2+pi^2)*v(w).^(1-d/pi) - ell^2;
c = 20;
b = (20-ell^2)*tan(d);
fill([-c c c -c -c -ell^2 -c], [-c -c c c b 0 -b], grey, ...
    'edgecolor', 'none'); hold on
for k = 1:N
    plot(z(ww{k}), '--k', LW, 1); hold on
end
plot(z(wp), LW, lw, 'color', blue); hold on
plot(z(wm), LW, lw, 'color', red); 
plot(z(w), '--k', LW, lw)
plot(-ell^2, 0, 'o', 'color', c1, MFC, c1);
plot(pi^2, 0, 'o', 'color', c2, MFC, 'w');
axis equal
axis([-20 20 -20 20])
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca, 'box', 'off')
axis off
circ = exp(1i*linspace(d,-d, 1000)+1i*pi);
ang = -ell^2+7*circ;
plot(ang, 'k', LW, 1)
plot(real(ang(end))+[-.5 0 .125], imag(ang(end))+[-.25 0 -.5], '-k', LW, 1)
% plot([-c, -ell^2], [0 0], '--k')
hold off
print -depsc2 ~/Dropbox/work/2015/contour/paper/figures/sectoralc

alignfigs

