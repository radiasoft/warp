#include "top.h"
c=============================================================================
c@(#) File w3d_load.F, version $Revision: 1.1 $, $Date: 2004/11/30 18:42:33 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is a source code file of package W3D of the WARP 
c  3d electrostatic PIC code, Cartesian geometry, for beam problems
c=============================================================================
c=============================================================================
      real(kind=8) function bessi0(x)
      real(kind=8):: x

c     bessi0(x) calculates the Bessel function I_0(x) for real x using 
c     approximate polynomial coefficients.  Because I_0 is smooth, this 
c     provides a good, efficient approximation using only a few terms and
c     together with recusion relations can be used to calculate higher order 
c     I_n(x) with integer n > 0.    
c     The routine is slightly modified from the one presented in:
c        Numerical Recipes The Art of Scientific Computing, 1st Edition 
c        Press, Flannery, Teukolsky, and Vetterling (1986)
 
      real(kind=8):: P1 =  1.0e0
      real(kind=8):: P2 =  3.5156229e0
      real(kind=8):: P3 =  3.0899424e0
      real(kind=8):: P4 =  1.2067492e0
      real(kind=8):: P5 =  0.2659732e0
      real(kind=8):: P6 =  0.360768e-1
      real(kind=8):: P7 =  0.45813e-2
      real(kind=8):: Q1 =  0.39894228e0
      real(kind=8):: Q2 =  0.1328592e-1
      real(kind=8):: Q3 =  0.225319e-2
      real(kind=8):: Q4 = -0.157565e-2 
      real(kind=8):: Q5 =  0.916281e-2
      real(kind=8):: Q6 = -0.2057706e-1
      real(kind=8):: Q7 =  0.2635537e-1
      real(kind=8):: Q8 = -0.1647633e-1
      real(kind=8):: Q9 =  0.392377e-2
      
      real(kind=8):: Y, AX

      if ( abs(x) .lt. 3.75) then
        Y = (x/3.75)**2
        bessi0 = P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      else
        AX=abs(x)
        Y=3.75/AX
        bessi0 = (exp(AX)/sqrt(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     &      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      endif
 
      return
      end
c=============================================================================
      real(kind=8) function bessi1(x)
      real(kind=8):: x

c     bessi1(x) calculates the Bessel function I_1(x) for real x using 
c     approximate polynomial coefficients.  Because I_1 is smooth, this 
c     provides a good, efficient approximation using only a few terms and
c     together with recusion relations can be used to calculate higher order 
c     I_n(x) with integer n > 0.    
c     The routine is slightly modified from the one presented in:
c        Numerical Recipes The Art of Scientific Computing, 1st Edition 
c        Press, Flannery, Teukolsky, and Vetterling (1986)

      real(kind=8):: P1 =  0.5e0
      real(kind=8):: P2 =  0.87890594e0
      real(kind=8):: P3 =  0.51498869e0
      real(kind=8):: P4 =  0.15084934e0
      real(kind=8):: P5 =  0.2658733e-1
      real(kind=8):: P6 =  0.301532e-2
      real(kind=8):: P7 =  0.32411e-3
      real(kind=8):: Q1 =  0.39894228e0 
      real(kind=8):: Q2 = -0.3988024e-1 
      real(kind=8):: Q3 = -0.362018e-2
      real(kind=8):: Q4 =  0.163801e-2 
      real(kind=8):: Q5 = -0.1031555e-1
      real(kind=8):: Q6 =  0.2282967e-1
      real(kind=8):: Q7 = -0.2895312e-1 
      real(kind=8):: Q8 =  0.1787654e-1
      real(kind=8):: Q9 = -0.420059e-2
      
      real(kind=8):: Y, AX

      if (abs(x) .lt. 3.75) then
        Y=(x/3.75)**2
        bessi1 = x*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      else
        AX=abs(x)
        Y=3.75/AX
        bessi1 = (exp(AX)/sqrt(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+
     &      Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
        if (x .lt. 0.) bessi1 = -bessi1
      endif

      return
      end
c=============================================================================
      real(kind=8) function bessi(n,x)
      integer(ISZ):: n 
      real(kind=8):: x

c     bessi(n,x) calculates the Bessel function I_n(x) for real x 
c     and integer n >= 2 using downward recursion operators seeded using 
c     calls to bessi0(x) for I_0(x).  
c     The routine is slightly modified from the one presented in:
c        Numerical Recipes The Art of Scientific Computing, 1st Edition 
c        Press, Flannery, Teukolsky, and Vetterling (1986)

      integer(ISZ):: J,M
      real(kind=8):: BI,BIM,BIP,TOX,bessi0

      integer(ISZ):: IACC=40
      real(kind=8):: BIGNO = 1.0e10
      real(kind=8):: BIGNI = 1.0e-10

      if ((n .lt. 2) .or. mod(n,1) .ne. 0.) then 
        call remark('bad argument n in bessi(n,x)')
        call kaboom(1) 
      endif 
      
      if (x .eq. 0.) then 
        bessi = 0. 
      else
        TOX=2.0/x
        BIP=0.0
        BI=1.0
        bessi=0.
        M=2*((N+INT(sqrt(float(IACC*N)))))
        do J=M,1,-1
          BIM=BIP+float(J)*TOX*BI
          BIP=BI
          BI=BIM
          if (abs(BI) .gt. BIGNO) then
            bessi=bessi*BIGNI
            BI=BI*BIGNI
            BIP=BIP*BIGNI
          endif 
          if (J .eq. N) bessi=BIP
        enddo
        bessi = bessi*bessi0(x)/BI
        if (x .lt. 0. .and. mod(n,2) .eq. 1) bessi = -bessi 
      endif 

      return
      end
c=============================================================================
