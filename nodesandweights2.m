function [z, w, h] = contourIntegral2(n, x, l2, h, a)
%CONTOURINTEGRAL  Nodes and weights for contour integral method.
% Inputs:   n  - Number of points.
%           x  - Height of required slice.
%           l2 - Largest eigenvalue of A.
% Ouputs:   z  - Quadrature nodes
%           w  - Quadrature weights.
%           h  - Optimal step size.
% For wapr.m see www.mathworks.com/matlabcentral/fileexchange/3644/
if ( nargin < 2 ), x = 0.5; end                % Default to x = 0.5.
if ( nargin < 3 ), l2 = 0;  end                % Default to l2 = 0.
if ( nargin < 4 || isempty(h))
%     h = pi/sqrt(n)-log(sqrt(pi^2+l2))/n;
    h = pi/sqrt(n);
end
if ( nargin < 5 ), a = 0; end                  % Default to a = 0.
t = ((0:n-1).'+0.5)*h;                         % Equally spaced.
z = (pi^2-l2)/2+1i*(pi^2+l2)/2*sinh(t+1i*a);   % Mapped nodes.
dzdt = 1i*(pi^2+l2)/2*cosh(t+1i*a);            % dz/dt
E = sin(x*sqrt(z))./sin(sqrt(z));              % E(x;z(t)) from (?.?)
idx = imag(sqrt(z)) > asinh(realmax);          % Deal with overflow ...
E(idx) = exp(1i*(1-x)*sqrt(z(idx)));           %  when im(z) >> 1.
w = (h/pi)*dzdt.*E;                            % Weights.
w(isnan(w)) = 0;
end
