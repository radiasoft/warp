
#define MACHEPS 1.0e-14
#define LARGEPOS 1.0e+36
#define SMALLPOS 1.0e-36
c Length of diagnostic control arrays
#define NCONTROL 50
c Default length of arrays describing lattice
#define NELEMENT 100
c Number of particle per group
#define NPARPGRP 256
c Max number of ptcl "subsets" for scatter plots
#define NSUBSETS 3
c Max number of diagnostic windows
#define NWINDOWS 9
c Number of z moments
#define NUMZMMNT 28
c Max number of lab windows
#define MAXNUMLW 50
c Used only for data statements
#define TNWINM  2*NWINDOWS
c Used only for data statements
#define NWPNSP1 NWINDOWS + NSUBSETS + 1
#define NEVER   0
#define SELDOM  1
#define ALWAYS  2
#define ERR 1
#define OK 0
#define STDOUT 6
#define dvnz(X) sign(abs(X)+SMALLPOS,X)
