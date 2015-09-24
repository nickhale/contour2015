
clc, close all

LW = 'LineWidth';         lw = 3;
MS = 'MarkerSize';        ms = 20;
FS = 'FontSize';          fs = 12;
MFC = 'MarkerFaceColor';  fc = 'w';

% Params:
b = 1;
x = .5;
ell2 = (pi/b)^2;

% n and h range:
nn = 1:40;
hh = linspace(0, 1, 50);

% Test points:
y = linspace(0, b, 100);
y = y(2:end-1);

%% 'Exact' solution (many terms in contour series)
n = 100;
% Optimal contour:
[z, w] = nodesandweights(n, x, ell2);
%  Solve resolvent systems
u = 0*y; 
for k = n:-1:1
    zk = z(k);
    v = 1./zk*(1-cosh((.5*b-y)*sqrt(zk))./cosh(b*.5*sqrt(zk)));
    v(isinf(v) | isnan(v)) = 0;
    u = u + imag(v*w(k));
end
sol = u;

%% Solve for each h and n:
err = zeros(length(hh), length(nn));
for jj = 1:numel(nn);
    for kk = 1:numel(hh);
        disp([jj, kk])
        % Define optimal contour
        [z, w] = nodesandweights_advanced(nn(jj), x, ell2, hh(kk));
        % Solve resolvent systems
        u = 0*y; 
        for k = 1:numel(w)
            zk = z(k);
            v = 1./zk*(1-cosh((.5*b-y)*sqrt(zk))./cosh(b*.5*sqrt(zk)));
            v(isinf(v) | isnan(v)) = 0;
            u = u + imag(v*w(k));
        end
        err(kk,jj) = norm((u - sol)./sol, inf);
    end
end

%% Once more with the optimal contour:
errOpt = zeros(size(nn));
hOpt = zeros(size(nn));
for jj = 1:numel(nn);
    disp(jj)
    % Define optimal contour
    [z, w, hOpt(jj)] = nodesandweights(nn(jj), x, ell2);
    % Solve resolvent systems
    u = 0*y; 
    for k = 1:numel(w)
        zk = z(k);
        v = 1./zk*(1-cosh((.5*b-y)*sqrt(zk))./cosh(b*.5*sqrt(zk)));
        v(isinf(v) | isnan(v)) = 0;
        u = u + imag(v*w(k));
    end
    errOpt(jj) = norm((u - sol)./sol,inf);
end

%%

figure(1)

[nnn, hhh] = meshgrid(nn, hh);

% To clean up noisey contours for large h.
err(hhh>.81) = max(err(hhh>.81), 1e-4);
err(err<1e-16) = 1e-16;
err(end) = 1e-24;

contourf(nnn, hhh, log10(err), -18:2:0, 'LineWidth', 2); hold on
contourf(nnn, hhh, log10(err), -16:2:0, 'LineWidth', 2);
plot3(nn, hOpt, 1e6+log10(errOpt), '--w', 'LineWidth', 4);
h = colorbar;
set(h, 'limits', [-16, 0]);
grid off, axis on
set(gca, 'fontsize', 14)
caxis([-16, 0])
axis([0 40 0, 1])
view(0, 90)

% Text:
su = .025;
sr = 0.25;
if ( x == .5 )
    dc = 'k'; ds = 'o'; ms = 5;
    p1 = [9.5, .787];
    text(p1(1), p1(2), '-6', FS, fs, 'color', fc);
    plot3(p1(1)+sr, p1(2)+su, 1e6, ds, 'color', dc, MFC, dc, MS, ms)
    p2 = [13.7, .595];
    text(p2(1), p2(2), '-8', FS, fs, 'color', fc);
    plot3(p2(1)+sr, p2(2)+su, 1e6, ds, 'color', dc, MFC, dc, MS, ms)
    p3 = [19.5, .448];
    text(p3(1), p3(2), '-10', FS, fs, 'color', fc);
    plot3(p3(1)+sr, p3(2)+su, 1e6, ds, 'color', dc, MFC, dc, MS, ms)
    p4 = [25, .365];
    text(p4(1), p4(2), '-12', FS, fs, 'color', fc);
    plot3(p4(1)+sr, p4(2)+su, 1e6, ds, 'color', dc, MFC, dc, MS, ms)
    p5 = [30.5, .308];
    text(p5(1), p5(2), '-14', FS, fs, 'color', fc);
    plot3(p5(1)+sr, p5(2)+su, 1e6, ds, 'color', dc, MFC, dc, MS, ms)
    p6 = [45-8.75, .263];
    text(p6(1), p6(2), '-16', FS, fs, 'color', fc);
    plot3(p6(1)+sr, p6(2)+su, 1e6, ds, 'color', dc, MFC, dc, MS, ms)
end
hold off

eval(['print -dpng ../paper/figures/contourPlotNH_x=', num2str(100*x)])

%%

figure(2)

sol = norm(sol, inf);
TE = exp(-(1-x)*sqrt(pi^2+ell2)*exp(nnn.*hhh/2)/(2*sqrt(2)));
DE = exp(-pi^2./hhh) + exp(-pi^2./hhh - (1-x)*sqrt(ell2)/sqrt(2));
err2 = (DE + TE)./sol;
err2(err2<1e-16) = 1e-16;
[c, ~] = contourf(nnn, hhh, log10(err2), -16:2:0, 'LineWidth', 2); hold on
plot3(nn, hOpt, 1e6+log10(errOpt), '--w', 'linewidth', 4); hold off

h = colorbar;
set(h, 'limits', [-16, 0]);
grid off, axis on
set(gca, 'fontsize', 14)
caxis([-16, 0])
axis([0 40 0, 1])
view(0, 90)

eval(['print -dpng ../paper/figures/contourPlotNH_x=', num2str(100*x),'_theory'])

%%

alignfigs