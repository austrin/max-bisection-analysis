% Phi^{-1}(x), where Phi(x) is the cdf of a standard N(0,1) Gaussian
function phiinv=normalinv(x)
    phiinv = sqrt(2)*erfinv(2*x-1);
