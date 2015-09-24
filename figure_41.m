% Solve the PIAS problem using Hess/Schur/US.
% Nick Hale, 27 May 2015.

clc, close all

x = .5;                            
n = 100;                           % n = number of resolvent solves  
mm = ceil(logspace(1, 3, 50));    % m = dimension of DM-2

t_direct(mm(end)) = 0;
t_hess(mm(end)) = 0;
t_schur(mm(end)) = 0;
t_us(mm(end)) = 0;

for m = [mm(1) mm]
    disp(m)

    %% Collocation:

    % A = diffmat(m+2, 2, [0 1]); % }
    % A([1,end],:) = [];          % } <-- Chebfun
    % A(:,[1,end],:) = [];        % }
    % y = chebpts(m+2);           % }

    [y, A] = cheb2bc(m+2, [1 0 0 ; 1 0 0]); % } <-- DMSUITE
    A = 4*A;                                % }

    I = eye(m);

    %% Define contour
    [z, w] = nodesandweights(n, x, 0);

    %% Direct

    u  = zeros(m,1); 
    u0 = ones(m,1);
    tic
    for k = 1:n
        v = (z(k)*I - A) \ u0;
        u = u + v*w(k);
    end
    u = imag(u);
    t_direct(m) = toc;
%     u05 = u(m/2+.5);

    %% Using Hess:

    u   = zeros(m, 1); 
    u0  = ones(m, 1);
    tic
    [P, H] = hess(A);
    Pu0 = P'*u0;
    for k = 1:n
        v = ( (z(k)*I - H) \ Pu0 );
        u = u + v*w(k);
    end
    u = imag(P*u);
    t_hess(m) = toc;
%     u05 = u(m/2+.5)

%% Using Schur complex:

    u   = zeros(m, 1); 
    u0  = ones(m, 1);
    tic
    [U, T] = schur(A, 'complex');
    Uu0 = U'*u0;
    for k = 1:n
        v = ( (z(k)*I - T) \ Uu0 );
        u = u + v*w(k);
    end
    u = imag(U*u);
    t_schur(m) = toc;
%     u05 = u(m/2+.5)

%% Using US:

    D2 = 4*ultraS.diffmat(m, 2);
    S = ultraS.convertmat(m, 0, 1);
    u   = zeros(m, 1); 
    u0  = [ 1 ; zeros(m-1, 1) ]; 
    rhs = S*u0;

    % Boundary conditions:
    B = ones(2, m); B(1,2:2:end) = -1;
    rhs = [ 0 ; 0 ; rhs(1:m-2) ];

    A = sparse(m, m);
    A(1:2,1:2) = [1 -1 ; 1 1];
    A(3:m,:) = -D2(1:m-2,:);
    S(3:m,:) = S(1:m-2,:); S(1:2,:) = 0;

    C = speye(2);
    U = zeros(m,2); U(1:2,1:2) = C;
    V = B; V(1:2,1:2) = 0;
    b = [rhs, U];
    tic
    for k = 1:n
        tmp = (A + z(k)*S)\b;
        Vtmp = V*tmp;
        v = tmp(:,1) - tmp(:,2:3)*((C + Vtmp(:,2:3))\(Vtmp(:,1)));
        u = u + v*w(k);
    end
    u = imag(u);
    t_us(m) = toc;
    
    % Evaluate at y = .5; 
%     u05 = u; u05(2:2:end) = []; u05(2:2:end) = -u05(2:2:end); u05 = sum(u05)

end


%%
ca, clc

cols = get(0, 'DefaultAxesColorOrder');
style = {'-s', '-^', '-v', '-d'};
MFC = 'MarkerFaceColor';
LW = 'LineWidth'; lw = 3;
mm2 = [mm(1:5:end) mm(end)];
MS = 'MarkerSize'; 

h1 = loglog(NaN, NaN, style{1}, 'color', cols(1,:), MFC, 'w', LW, 2); hold on
loglog(mm, t_direct(mm), LW, lw, 'color', cols(1,:)); hold on
loglog(mm2, t_direct(mm2), style{1}(2), LW, 2, MS, 8, 'color', cols(1,:), MFC, 'w');


h2 = loglog(NaN, NaN, style{2}, 'color', cols(2,:), MFC, 'w', LW, 2);
loglog(mm, t_hess(mm), LW, lw, 'color', cols(2,:)); hold on
loglog(mm2, t_hess(mm2), style{2}(2), LW, 2, MS, 8, 'color', cols(2,:), MFC, 'w');

h3 = loglog(NaN, NaN, style{3}, 'color', cols(5,:), MFC, 'w', LW, 2);
loglog(mm, t_schur(mm), LW, lw, 'color', cols(5,:)); hold on
loglog(mm2, t_schur(mm2), style{3}(2), LW, 2, MS, 8, 'color', cols(5,:), MFC, 'w');

h4 = loglog(NaN, NaN, style{4}, 'color', cols(3,:), MFC, 'w', LW, 2);
loglog(mm, t_us(mm), LW, lw, 'color', cols(3,:)); hold on
loglog(mm2, t_us(mm2), style{4}(2), LW, 2, MS, 8, 'color', cols(3,:), MFC, 'w');

if ( n == 1 )
    legend([h1 h2 h3 h4], 'Direct', 'Hessenberg', 'Schur', 'Ultraspherical', 'location', 'NW')
end

grid on
xlim([8 1200])
ylim([5e-5 2e0])

% set(gca, 'FontSize', 16);
% pause(1)
% str = ['~/Dropbox/work/2015/contour/paper/figures/pias_time_' int2str(n)];
% print('-depsc2', str)
% savefig(str)

return
%%
ca
n = 1;
str = ['~/Dropbox/work/2015/contour/paper/figures/pias_time_' int2str(n)];
openfig(str)
ylim([5e-5 2e1])
print('-depsc2', str)

get(gca)
% savefig(str)