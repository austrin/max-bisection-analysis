This directory contains some heuristic tools for computing
approximation ratios, and finding the best rounding functions.

The code is written in Matlab.  It might work in Octave but this has
not been tested.

## Example usage:

1. To try out the linear algorithm:
```
>> format long
>> c = 0.864503;
>> alpha_linear(c)

ans =

   0.873682718504500

>> alpha(0.177, 0.177, -0.646, c*0.177, c*0.177)

ans =

   0.873682889846272
```

2. To try out the pairing algorithm:
```
>> c = 0.8056;
>> f = inline('1.618*max(0, x-0.478)');
>> alpha_pairing(c, f)

ans =

   0.877629684815347

>> alpha_cf(0.177, 0.177, -0.646, c, f)

ans =

   0.877988468767727
```

## Scripts

The scripts here are as follows.  Some of the scripts take additional
parameters specifying precision and verbosity, see the respective
scripts for details.

alpha*.m - the various alpha functions from the paper:

 * alpha.m: the "plain" alpha taking five parameters mu1, mu2, rho, r1, r2.
 * alpha_linear.m: the function alpha(c) from the linear algorithm,
   		   finding the worst configuration for a given c
 * alpha_I.m: the variant of alpha from the pairing algorithm where r1
              and r2 are replaced by two intervals I1 and I2
 * alpha_cf.m: the function alpha_{c,f}(mu1, mu2, rho) from the pairing algorithm.
 * alpha_pairing.m: the function alpha(c, f) from the pairing algorithm,
                    finding the worst configuration for a given c and f

gendata*.m - various functions to generate data to be plotted:

 * gendata_linear_limit.m: compute the approximation ratio on the two
                           worst configurations (and their mix) as a
                           function of c.
 * gendata_pairing_contour.m: compute the approximation ratio of
                              pairing algorithm along a fine grid on
                              the boundary of the polytope.

The rest are helper functions:

 * Phi.m - normal distribution cdf
 * normalinv.m - inverse of Phi
 * bnvl.m - bivariate normal distribution cdf
 * Gamma.m - as in paper
 * Lambda.m - as in paper
 * trho.m - \tilde{\rho}(mu1, mu2, rho) = (rho-mu1*mu2)/sqrt((1-mu1^2)(1-mu2^2))
            as in definition of alpha
 * random_configuration.m - returns a random configuration satisfying
                            the triangle inequalities
