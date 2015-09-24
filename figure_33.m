close all

LW = 'LineWidth';         lw = 3;
MS = 'MarkerSize';        ms = 20;
MFC = 'MarkerFaceColor';  mfc = 'k';
FC = 'color';             fc = 'w';
FS = 'FontSize';          fs = 12;

% Params:
b = .1;
x = .5;
ell2 = (pi/b)^2;

% n and ell^2 range:
nn = 1:40;
ll2 = linspace(0, 1.5*ell2, 150);

% Test points:
y = linspace(0, b, 100); y = y(2:end-1);

%% 'Exact' solution (many terms)

n = 70;
% Optimal contour
[z, w] = nodesandweights(n, x, ell2);
u = 0*y; 
%  Solve resolvent systems
for k = n:-1:1
    zk = z(k);
    v = 1./zk*(1-cosh((.5*b-y)*sqrt(zk))./cosh(b*.5*sqrt(zk)));
    u = u + imag(v*w(k));
end
sol = u;

%%

hh = [];
for jj = 1:numel(nn);
    for kk = 1:numel(ll2);
        [jj, kk]
        n = nn(jj);
        % Optimal contour
        [z, w, h] = nodesandweights(n, x,  ll2(kk));
        %  Solve resolvent systems
        u = 0*y; 
        for k = n:-1:1
            zk = z(k);
            v = 1./zk*(1-cosh((.5*b-y)*sqrt(zk))./cosh(b*.5*sqrt(zk)));
            u = u + imag(v*w(k));
        end
        hh(kk,jj) = h;
        err(kk,jj) = norm((u - sol)./sol, inf);
    end
end

%%

for jj = 1:numel(nn);
    n = nn(jj);
    % Optimal step size
    [z, w] = nodesandweights(n, x, ell2);    
    % Solve resolvent systems
    u = 0*y; 
    for k = n:-1:1
        zk = z(k);
        v = 1./zk*(1-cosh((.5*b-y)*sqrt(zk))./cosh(b*.5*sqrt(zk)));
        u = u + double(v*w(k));
    end
    u = 2*h/pi*imag(u);

    errOpt(kk,jj) = norm((u - sol)./sol, inf);
end


%%

sol = norm(sol, inf);
err(~err | err < 1e-14) = 1e-14;

[c, ~] = contourf(nn, ll2, log10(err), -14:2:0, 'lineWidth', 2); hold on
h = colorbar;
set(h, 'limits', [-14 0]);
caxis([-14,0])
plot3(nn, ell2+0*nn, 1e10+log10(errOpt), '--w', 'linewidth', 3);

su = 20;
sr = -0.4;
if ( x == .5 )
    ds = '.'; dc = 'k';
    p1 = [4.45, 1400];
    text(p1(1), p1(2), '-2', FS, fs, FC, fc);
    plot(p1(1)+sr, p1(2)+su, ds, FC, dc, MFC, mfc, MS, ms)
    p2 = [10.8, 1360];
    text(p2(1), p2(2), '-4', FS, fs, FC, fc);
    plot(p2(1)+sr, p2(2)+su, ds, FC, dc, MFC, mfc, MS, ms)
    p3 = [16, 1300];
    text(p3(1), p3(2), '-6', FS, fs, FC, fc);
    plot(p3(1)+sr, p3(2)+su, ds, FC, dc, MFC, mfc, MS, ms)
    p4 = [20.25, 1245];
    text(p4(1), p4(2), '-8', FS, fs, FC, fc);
    plot(p4(1)+sr, p4(2)+su, ds, FC, dc, MFC, mfc, MS, ms)
    p5 = [23.6, 1195];
    text(p5(1), p5(2), '-10', FS, fs, FC, fc);
    plot(p5(1)+sr, p5(2)+su, ds, FC, dc, MFC, mfc, MS, ms)
    p6 = [26, 1145];
    text(p6(1), p6(2), '-12', FS, fs, FC, fc);
    plot(p6(1)+sr, p6(2)+su, ds, FC, dc, MFC, mfc, MS, ms)
end

set(gca, FS, 14)
view(0,90);
xlim([0 40])
grid off, shg

%%

% eval(['print -dpng ../paper/figures/contourPlotM_x=', num2str(100*x)])

%%

figure(2)

[NNN, mmm] = meshgrid(nn, ll2);

zz = .5*(pi^2-mmm) + 1i*(pi^2+mmm).*sinh(NNN.*hh);
TE = abs(sin(x*sqrt(zz))./sin(sqrt(zz))./(zz+mmm).*(pi^2+mmm).*cosh(NNN.*hh));
TE(isnan(TE)) = 0;

zInv = @(w,m) asinh(((w+m)./((m+pi^2)/2)-1)/1i);
ddd = imag(zInv(-ell2, mmm));
DE = exp(-2*pi*ddd./hh);
DEm = exp(-pi^2./hh);
DEp = exp(-(1-x)*sqrt(mmm/2)).*exp(-2*pi*ddd./hh);
errTheory = (DEp + DEm + TE)/sol;
errTheory(~errTheory | errTheory < 1e-14) = 1e-14;

[c, ~] = contourf(nn, ll2, log10(errTheory), -14:2:0, 'lineWidth', 2); hold on
h = colorbar;
set(h, 'limits', [-14 0]);
caxis([-14,0])
plot3(nn, ell2+0*nn, 1e10+log10(errOpt), '--w', 'linewidth', 3); hold off

set(gca, FS, 14)
view(0,90);
xlim([0 40])
grid off

%%

% eval(['print -dpng ../paper/figures/contourPlotM_x=', num2str(100*x), '_theory'])

%%

% alignfigs
