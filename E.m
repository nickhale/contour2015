function y = E(x, z)
% E(X, Z) evaluates SIN(X*SQRT(Z))/SIN(SQRT(Z)) accurately for all values of X
% in [0, 1] and large magnitude (complex-valued) Z. In particular it avoids both
% the risk of overflow and rounding error as z ~ 1.

if ( x == 1 || x == 0 )
    % Trivial cases:
    y = x + 0*z;
    
elseif ( x < .9 )
    
    % Easy case:
    y = sin(x*sqrt(z))./sin(sqrt(z));
    
else
    
    % Tricky case (i.e., z close to 1):
%     y = cos((1-x)*sqrt(z)) - cot(sqrt(z)).*sin((1-x)*sqrt(z));
    y = sin(x*sqrt(z))./sin(sqrt(z));

%     yv = double(sin(vpa(x)*sqrt(vpa(z)))./sin(sqrt(vpa(z))));
%     norm(y - yv, inf)   
    
    idx = imag(sqrt(z)) > asinh(realmax);
    y(idx) = exp(1i*(1-x)*sqrt(z(idx)));

%     norm(y - yv, inf)
    
end

    idx = isnan(y);
    sz = sign(imag(z(idx))).*sqrt(z(idx));
    y(idx) = exp(1i*(1-x)*sz);

end