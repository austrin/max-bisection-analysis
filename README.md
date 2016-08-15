# Analysis of Max Bisection Algorithms

This repository provides some programs to analyze pairing-based
approximation algorithms for Max Bisection introduced in the paper
[Better Balance by Being Biased: a 0.8776-Approximation for Max
Bisection](http://arxiv.org/abs/1205.0458).

Two sets of programs are provided:

1. [proof/](proof) contains a program that is used to formally
   prove the approximation ratio for a given choice of rounding
   algorithm.

2. [heuristics/](heuristics) contains programs used to heuristically
   compute approximation ratios (using local optimization), which is
   much faster than the prover and can be used to experiment with
   different rounding algorithms.

For more details, see the paper.
