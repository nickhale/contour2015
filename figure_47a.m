%%
% Parameters:
r1 = 1; r2 = 3;             % Annulus radii
N = 21;                     % Number of quadrature nodes
x = 0.5;
LW = 'LineWidth'; lw = 2;
m = 31;                     % Discretisation in r and t
u1 = @(r, t) (r-r1).*(r2-r).*(1-sin(t));
%% 

tic

% Differential operators:
[t, D2t] = fourdif(m-2, 2);
[r, D2r, D1r, phip] = cheb2bc(m, [1 0 0 ; 1 0 0]);

% Scaling r from [-1 1] to [r1 r2]:
scl = (r2-r1)/2;
r = scl*r+(r2+r1)/2;
D1r = 1/scl*D1r;
D2r = 1/scl^2*D2r;
R = diag(r); R2 = diag(r.^2);
Ri = diag(1./r); Ri2 = diag(1./r.^2);
A = D2t;
B = (R*D1r + R2*D2r).';

% RHS
[rr, tt] = meshgrid(r, t);
U1 = u1(rr,tt)*R2; % Scaled RHS.

% Compute Schur factorization of A.
[Z1, A] = schur(A, 'complex');
% Transform the righthand side.
U1 = Z1'*U1;

% Define optimal contour and quadrature nodes/weights:
ell2 = 0;
[z, w, h] = nodesandweights(N, x, ell2);

% Solve resolvent systems
I = speye(m-2, m-2);
v = cell(N,1);
for k = 1:N
    v{k} = lyap(A, B - z(k)*R2, U1);
end

xx = linspace(0, x, 50);
u = cell(numel(xx),1);
z0 = zeros(m-1,1);
th = ((0:N-1).' + 0.5)*h;                % Equally spaced.
mu = ell2 + pi^2;
z = 0.5*mu*(1+1i*sinh(th)) - ell2;       % Mapped nodes.
dzdt = 0.5i*mu*cosh(th);
U = [];
for j = 1:numel(xx)
    xj = xx(j);
    E = sin(xj*sqrt(z))./sin(sqrt(z));
    w = (h/pi)*dzdt.*E;                     % Weights.
    w(isnan(w)) = 0;
    
    uj = 0;
    for k = 1:N
        uj = uj + w(k)*v{k};
    end
    uj = imag(Z1*uj);
    % Append BCs
    uj(end+1,:) = uj(1,:);    %#ok<SAGROW>
    uj = [z0, uj, z0];        %#ok<AGROW>
    u{j} = uj;
    U(:,j) = uj(:,(m+1)/2);
end

toc

% Plotting
[X,Y,Z] = cylinder(2*ones(length(xx),1), m-2);
h1 = surf(X,Y,Z, 'edgealpha', 0);
set(h1, 'CData', U.');
set(h1, 'ZData', get(h1, 'ZData')/2);
axis square, shading interp
tmp = get(h1, 'XData');
set(h1, 'XData', get(h1, 'ZData'));
set(h1, 'ZData', tmp);
caxis(gca, [0 1]), colorbar
view([23, 11])
set(gca, 'xtick', [0 .5 1])
set(gca, 'ytick', [-3 3])
set(gca, 'ztick', [-3 3])
xlim([0 1])
ylim([-3 3])
zlim([-3 3])

hold on
xx = linspace(0, 1, 100);
th = .9*pi/2;
plot3(xx, 3*cos(th)+0*xx, 3*sin(th)+0*xx, 'k', 'LineWidth', lw);
plot3(xx, cos(th)+0*xx, sin(th)+0*xx, ':k', 'LineWidth', lw);
th = .95*3*pi/2;
plot3(xx, 3*cos(th)+0*xx, 3*sin(th)+0*xx, 'k', 'LineWidth', lw);
plot3(xx, cos(th)+0*xx, sin(th)+0*xx, ':k', 'LineWidth', lw);
plot3(0*xx+1, cos(2*pi*xx), sin(2*pi*xx), 'k', LW, 3)
plot3(0*xx+1, 3*cos(2*pi*xx), 3*sin(2*pi*xx), 'k', LW, 3)
plot3(0*xx, cos(2*pi*xx), sin(2*pi*xx), 'k', LW, 3)
plot3(0*xx, 3*cos(2*pi*xx), 3*sin(2*pi*xx), 'k', LW, 3)
% plot3(.5+0*xx, cos(2*pi*xx), sin(2*pi*xx), ':k', LW, 3)
grey = .7*[1 1 1];
fill3([0*xx 0*xx] + 1, [3*cos(2*pi*xx) cos(2*pi*xx)], ...
    [3*sin(2*pi*xx) sin(2*pi*xx)], grey, 'facealpha', .5, 'edgealpha', 0)
hold off

toc

pause(1)
% print -dpng ../paper/figures/cyl_plot2
