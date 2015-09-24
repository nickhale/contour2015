
% Test PIAR for various h and N.

clc, close all, clear all

LW = 'LineWidth';         lw = 3;
MS = 'MarkerSize';        ms = 20;
FS = 'FontSize';          fs = 12;
MFC = 'MarkerFaceColor';  fc = 'k';

% Params:
b = .1;
x = .5;
ell = (pi/b);
ell2 = ell^2;
n = 15;

% a and h range:
hh = linspace(0, 1, 50);
aa = linspace(-pi/2, pi/2, 50);

% Test points:
y = linspace(0, b, 100);
y = y(2:end-1);


%% 'Exact' solution (many terms)

% Optimal contour:
nOpt = 100;
[z, w] = nodesandweights(nOpt, x, ell2);
% Solve resolvent systems
u = 0*y; 
for k = nOpt:-1:1
    zk = z(k);
    v = 1./zk*(1-cosh((.5*b-y)*sqrt(zk))./cosh(b*.5*sqrt(zk)));
    u = u + imag(v*w(k));
end
sol = u;

%% Solve for each h and a:
err = zeros(length(hh), length(aa));
for jj = 1:numel(aa);
    for kk = 1:numel(hh);
        [jj, kk]
        aStar = aa(jj);
        h = hh(kk);
        % Define optimal contour
        [z, w] = nodesandweights_advanced(n, x, ell2, h, aStar);
        % Solve resolvent systems
        u = 0*y; 
        for k = 1:numel(w)
            zk = z(k);
            v = 1./zk*(1-cosh((.5*b-y)*sqrt(zk))./cosh(b*.5*sqrt(zk)));
            v(isnan(v)) = 0;
            u = u + imag(v*w(k));
        end
        err(kk,jj) = norm((u - sol)./sol, inf);
    end
end


%%

[aaa, hhh] = meshgrid(aa, hh);
err(err < 1e-14) = 1e-14;
[~, ~] = contourf(aaa, hhh, log10(err), [-14:2:0], 'LineWidth', 2);
h = colorbar;
set(h, 'limits', [-14,0]);
ylim([0, 1])
grid off, axis on
set(gca, 'fontsize', 14)
caxis([-14,0])
set(gca, 'xtick', [-pi/2, 0, pi/2], 'xticklabel', {'-\pi/2', '0', '\pi/2'})

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

figure(2)

a = linspace(-pi/2, pi/2);
h = linspace(0, 1);
[aa, hh] = meshgrid(a, h);
TE = exp(-.5*(1-x)*sqrt(pi^2+ell2)*exp(n*hh/2).*sin(aa/2+pi/4));
DEm = exp(-2*pi*(pi/2+aa)./hh);
DEp = exp(-2*pi*(pi/2-aa)./hh-(1-x)*ell);
E = DEm + DEp + TE;

%%

contourf(aa, hh, log10((TE + DEp + DEm)/norm(sol,inf)), -14:2:0, 'LineWidth', 2);
colorbar
caxis([-14 0])
set(gcf, 'position', get(gcf, 'position') + [0 -560  0 0 ]);
set(gca, 'fontsize', 14)
caxis([-14,0])
set(gca, 'xtick', [-pi/2, 0, pi/2], 'xticklabel', {'-\pi/2', '0', '\pi/2'})

%%

h = @(a) 2/n*wapr(sqrt(2)*pi^2*n/((1-x)*sqrt(pi^2+ell2)));
% Original
a0 = 0;
h0 = h(a0);
% Nonzero alpha
aStar = ell*(1-x)*h0/(4*pi);

%%
ms = 8;
figure(1)
hold on
plot(a0, h0, 'ow', MFC, 'w', MS, ms)
plot(aStar, h0, '^w', MFC, 'w', MS, ms)
hold off

figure(2)
hold on
plot(a0, h0, 'ow', MFC, 'w', MS, ms)
plot(aStar, h0, '^w', MFC, 'w', MS, ms)
hold off

%%


figure(1)
eval(['print -dpng ../paper/figures/contourPlotAH_x=', num2str(100*x), '_n=', int2str(n), '_rect'])
figure(2)
eval(['print -dpng ../paper/figures/contourPlotAH_x=', num2str(100*x), '_n=', int2str(n), '_rect_theory'])

