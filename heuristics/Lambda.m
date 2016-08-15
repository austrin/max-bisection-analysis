% Lambda(r1, r2, trho) = 2Gamma((1-r1)/2, (1-r2)/2, trho) + (r1+r2)/2
%
% The arguments can be matrices, in which case they must all be of the
% same dimension.
%
% Requires:
%  -1 <= r1 <= 1
%  -1 <= r2 <= 1
%  -1 <= trho <= 1
function L = Lambda(r1, r2, trho)
  L = 2*Gamma((1-r1)/2, (1-r2)/2, trho) + (r1+r2)/2;
