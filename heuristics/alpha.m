% alpha(mu1, mu2, rho, r1, r2) = (2(1-Lambda(r1, r2, trho)))/(1-rho)
% where trho = (rho-mu1*mu2)/sqrt((1-mu1^2)(1-mu2^2))
%
% Evaluates a (set of) roundings on a (set of) configurations.
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
%   -1 <= r1 <= 1
%   -1 <= r2 <= 1
function a = alpha(mu1, mu2, rho, r1, r2)

tr = trho(mu1, mu2, rho);
a = (2*(1-Lambda(r1, r2, tr))) ./ (1-rho);
a(find(rho >= 1)) = 1;
a(find(a >= 1000000)) = 1000000; % Steer clear of degenerate places