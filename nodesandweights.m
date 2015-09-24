function [z, w, h] = nodesandweights(n, x, ell2)
%NODESANDWEIGHTS  Optimal nodes and weights for contour integral method.
% Inputs:   n  - Number of points
%           x  - Height of required slice
%           ell2 - Largest (most positive) eigenvalue of A is at -ell2
% Ouputs:   z  - Quadrature nodes
%           w  - Quadrature weights
%           h  - Optimal step size
% For wapr.m see www.mathworks.com/matlabcentral/fileexchange/3644/
if ( nargin < 3 ), ell2 = 0; end                % Default to ell2 = 0
h = 2/n*wapr(sqrt(2/(ell2+pi^2))*pi^2*n/(1-x)); % Optimal step size: Eqn (2.12)
t = ((0:n-1).'+0.5)*h;                          % Equispaced t_k (midpoint rule)
z = (pi^2-ell2)/2+1i*(pi^2+ell2)/2*sinh(t);     % Mapped nodes z(t_k): Eqn (2.1)
dzdt = 1i*(pi^2+ell2)/2*cosh(t);                % z'(t_k)
E = sin(x*sqrt(z))./sin(sqrt(z));               % E(x;z(t_k)): Eqn (1.1)
indx = imag(sqrt(z)) > .95*asinh(realmax);      % Deal with overflow ...
E(indx) = exp(1i*(1-x)*sqrt(z(indx)));          %  when im(z) >> 1.
w = (h/pi)*dzdt.*E;                             % Weights: Eqn (2.5)
end