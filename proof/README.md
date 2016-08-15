The prover program requires the Boost interval arithmetic and math
libraries, and the GNU MPFR library installed (in Ubuntu, the packages
libboost-dev and libmpfr4 can be used to provide these libraries).

The makefile produces the program prover, which checks that the
approximation ratio of the pairing algorithm on smooth configurations
is at least 0.87762.  The parameters used can be tweaked by editing
prover.cc.

Be aware that the prover uses around 1.5 GB of memory (with the
parameters configured here -- other parameters may result in much less
or much more memory), so if your machine is low on memory you will
probably suffer from excessive swapping.

Example usage:

```
terminal$ make
g++ -O3 -Wall -c -o alpha.o alpha.cc
g++ -O3 -Wall -c -o bvnl.o bvnl.cc
g++ -O3 -Wall -c -o configurations.o configurations.cc
g++ -O3 -Wall -c -o gaussian.o gaussian.cc
g++ -O3 -Wall -c -o search.o search.cc
g++ -O3 -Wall -o prover prover.cc alpha.o bvnl.o configurations.o gaussian.o search.o -lmpfr
terminal$ time ./prover
Searching in box ([-0.99999,0.99999], [-0.99999,0.99999], [-0.99999,0.99999])
Level 0: 1 cases. LB: 1. UB: 1. Bbox: ([-0.99999,0.99999], [-0.99999,0.99999], [-0.99999,0.99999])
Level 1: 8 cases. LB: 1. UB: 1. Bbox: ([-0.99999,0.99999], [-0.99999,0.99999], [-0.99999,0.99999])
Level 2: 38 cases. LB: 1. UB: 1. Bbox: ([-0.499995,0.99999], [-0.99999,0.99999], [-0.99999,0.99999])
Level 3: 168 cases. LB: 1. UB: 1. Bbox: ([-0.2499975,0.99999], [-0.99999,0.99999], [-0.99999,0.99999])
Level 4: 820 cases. LB: 0.87846614. UB: 1. Bbox: ([-0.12499875,0.99999], [-0.99999,0.99999], [-0.99999,0.99999])
Level 5: 3461 cases. LB: 0.87762028. UB: 0.96763867. Bbox: ([-0.062499375,0.99999], [-0.99999,0.99999], [-0.99999,0.99999])
Level 6: 12553 cases. LB: 0.87762028. UB: 0.9229909. Bbox: ([-0.031249688,0.99999], [-0.99999,0.99999], [-0.99999,0.99999])
Level 7: 51745 cases. LB: 0.87762. UB: 0.90035217. Bbox: ([-0.015624844,0.99999], [-0.99999,0.99999], [-0.99999,0.99999])
Level 8: 240079 cases. LB: 0.87762. UB: 0.88900159. Bbox: ([-0.0078124219,0.99999], [-0.99999,0.99999], [-0.99999,0.99999])
Level 9: 1117812 cases. LB: 0.87762. UB: 0.88332748. Bbox: ([-0.0039062109,0.99999], [-0.97655273,0.99999], [-0.97655273,0.99999])
Level 10: 4599038 cases. LB: 0.87762. UB: 0.88049061. Bbox: ([-0.0019531055,0.99999], [-0.94530305,0.99999], [-0.94530305,0.99999])
Level 11: 13971033 cases. LB: 0.87762. UB: 0.87907203. Bbox: ([-0.00097655273,0.99999], [-0.93163131,0.99999], [-0.93163131,0.99999])
Level 12: 16876426 cases. LB: 0.87762. UB: 0.8783629. Bbox: ([0,0.99999], [-0.92381889,0.99999], [-0.92381889,0.99999])
Level 13: 24662726 cases. LB: 0.87762. UB: 0.87800832. Bbox: ([0,0.99999], [-0.91747129,0.99999], [-0.91698302,0.99999])
Level 14: 28084949 cases. LB: 0.87762. UB: 0.87783104. Bbox: ([0,0.99999], [-0.90904853,0.99999], [-0.90868232,0.99999])
Level 15: 31226392 cases. LB: 0.87762. UB: 0.8777424. Bbox: ([0,0.99999], [-0.871024,0.99999], [-0.8708409,0.99999])
Level 16: 18690729 cases. LB: 0.87762. UB: 0.87769808. Bbox: ([0,0.99999], [-0.82033481,0.99999], [-0.82024326,0.99999])
Level 17: 432602 cases. LB: 0.87762. UB: 0.87769808. Bbox: ([0,0.99999], [-0.76935571,0.99999], [-0.76930993,0.99999])
Level 18: 507082 cases. LB: 0.87762. UB: 0.87769808. Bbox: ([0,0.99999], [-0.70256866,0.99999], [-0.70253814,0.99999])
Level 19: 65 cases. LB: 0.87762. UB: 0.87769808. Bbox: ([0,0.99999], [0,0.99998619], [0,0.99999])
Search complete. A total of 140477727 cases were inspected
LB: 0.87762. UB: 0.87769808.


real	24m3.252s
user	24m2.660s
sys	0m0.552s
