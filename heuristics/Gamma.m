% Gamma(x, y, rho) = Pr[X <= Phi^{-1}(x) and Y <= Phi^{-1}(y)]
% where (X,Y) are jointly normal random variables with mean 0, variance 1, and covariance rho
%
% The arguments can be matrices, in which case they must all be of the
% same dimension.
%
% Requires:
%   0 <= x <= 1
%   0 <= y <= 1
%  -1 <= rho <= 1
%
function G = Gamma(x, y, rho)
   x = trunc(x, 0, 1);
   y = trunc(y, 0, 1);
   rho = trunc(rho, -1, 1);
   t1 = trunc(normalinv(x), -1e10, 1e10);
   t2 = trunc(normalinv(y), -1e10, 1e10);
   G = zeros(size(x));
   for i=1:size(x,1)
       for j = 1:size(x,2)
           G(i,j) = bvnl(t1(i,j), t2(i,j), rho(i,j));
       end;
   end;


function X = trunc(X, lo, hi)
   X(find(X < lo)) = lo;
   X(find(X > hi)) = hi;
