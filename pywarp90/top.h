
#define MACHEPS 1.0e-14
#define LARGEPOS 1.0e+36
#define SMALLPOS 1.0e-36
! Length of diagnostic control arrays
#define NCONTROL 50
! Default length of arrays describing lattice
#define NELEMENT 100
! Number of particle per group
#define NPARPGRP 256
! Max number of ptcl "subsets" for scatter plots
#define NSUBSETS 3
! Max number of diagnostic windows
#define NWINDOWS 9
! Number of z moments
#define NUMZMMNT 28
! Max number of lab windows
#define MAXNUMLW 50
! Used only for data statements
#define TNWINM  2*NWINDOWS
! Used only for data statements
#define NWPNSP1 NWINDOWS + NSUBSETS + 1
#define NEVER   0
#define SELDOM  1
#define ALWAYS  2
#define STDOUT 6
#define dvnz(X) sign(abs(X)+SMALLPOS,X)

! Define size of integers. Must be the same as the size of a long int in C.
#ifndef ISZ
#if defined ALPHA || defined T3E || defined J90 || defined X86_64 || defined IA64
#define ISZ 8
#else
#define ISZ 4
#endif
#endif

! Define size of integers used in the MPI libraries.
#ifndef MPIISZ
#define MPIISZ 4
#endif

! Define size of word.
#ifndef WORDSIZE
#if defined T3E || defined J90
#define WORDSIZE 64
#else
#define WORDSIZE 32
#endif
#endif
