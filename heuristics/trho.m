% trho(mu1, mu2, rho) = (rho-mu1*mu2)/sqrt((1-mu1^2)(1-mu2^2))
%
% The arguments can be matrices, in which case they must all be of the
% same dimension.
%
% Requires:
%   -1 <= mu1 <= 1
%   -1 <= mu2 <= 1
%   -1 <= rho <= 1
%   -mu1-mu2-rho <= 1
%   -mu1+mu2+rho <= 1
%    mu1+mu2-rho <= 1
%    mu1-mu2+rho <= 1
%
function tr = trho(mu1, mu2, rho)
  d1 = (1-mu1.^2);
  d2 = (1-mu2.^2);
  tr = (rho - mu1.*mu2) ./ sqrt(d1.*d2);
  tr(find(d1 <= 0)) = 0;
  tr(find(d2 <= 0)) = 0;
  tr(find(tr >= 1)) = 1;
  tr(find(tr <= -1)) = -1;
