#include "top.h"
c=============================================================================
c@(#) File w3d_load.F, version $Revision: 1.3 $, $Date: 2004/12/02 23:55:28 $
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
c
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
c
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
c
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
        bessi = 0.
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
          if (J .eq. N) bessi = BIP
        enddo
        bessi = bessi*bessi0(x)/BI
        if (x .lt. 0. .and. mod(n,2) .eq. 1) bessi = -bessi 
      endif 

      return
      end
c=============================================================================
      real(kind=8) function zbrent(func,x1,x2,tol,itmax)

      integer(ISZ):: itmax 
      real(kind=8):: func,x1,x2,tol

c     zbrent caluclates the root of a func(x) that is bracketed between 
c     x1 and x2.  The root is found to accuracy tol using up to itmax 
c     iterations.  
c
c     The routine is slightly modified from the one presented in:
c        Numerical Recipes The Art of Scientific Computing, 1st Edition 
c        Press, Flannery, Teukolsky, and Vetterling (1986)

  
      integer(ISZ):: ITER
      real(kind=8):: A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM

      real(kind=8):: EPS = 3.e-8

      A=x1
      B=x2
      FA=func(A)
      FB=func(B)
      if(FB*FA .gt. 0.) then
        call remark('Root must be bracketed for zbrent.') 
        call kaboom(1)  
      endif 
      FC=FB
      do ITER=1,itmax
        if(FB*FC .gt. 0.) then
          C=A
          FC=FA
          D=B-A
          E=D
        endif
        if(abs(FC) .lt. abs(FB)) then
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        endif 
        TOL1=2.*EPS*abs(B)+0.5*tol
        XM=.5*(C-B)
        if(abs(XM) .le. TOL1 .or. FB .eq. 0.) then
          zbrent = B
          return
        endif
        if(abs(E) .ge. TOL1 .and. abs(FA) .gt. abs(FB)) then
          S=FB/FA
          if(A .eq. C) then
            P=2.*XM*S
            Q=1.-S
          else
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          endif
          if(P .gt. 0.) Q=-Q
          P=ABS(P)
          if(2.*P .lt. min(3.*XM*Q-abs(TOL1*Q),abs(E*Q))) then
            E=D
            D=P/Q
          else
            D=XM
            E=D
          endif
        else
          D=XM
          E=D
        endif
        A=B
        FA=FB
        if(abs(D) .gt. TOL1) then
          B=B+D
        else
          B=B+sign(TOL1,XM)
        endif
        FB=func(B)
      enddo
      zbrent = B
      call remark('zbrent exceeding maximum iterations.')
      call kaboom(1)

      return
      end
c=============================================================================
      real(kind=8) function watbagconstr(k1re)

      real(kind=8):: k1re

c    watbagconstr returns a scaled constraint equation assoicated with 
c    a matched waterbag equilbrium distribution function in a continuous 
c    focusing channel.  
c
c    k1re = k_1 * r_e       k_1 = wavenumber parameter 
c                           r_e = beam edge radius
c

      real(kind=8):: bessi0,bessi

      watbagconstr = bessi0(k1re)**2/bessi(2,k1re)**2 - 
     &  4./(k1re)*( 2*bessi(2,k1re)/bessi0(k1re) +
     &              k1re*bessi0(k1re)*bessi(3,k1re)/Bessi(2,k1re)**2
     &            )

      return 
      end 
c=============================================================================
      real(kind=8) function zbrent_test(x1,x2,tol,itmax)

      integer(ISZ):: itmax 
      real(kind=8):: func,x1,x2,tol

c     zbrent caluclates the root of a func(x) that is bracketed between 
c     x1 and x2.  The root is found to accuracy tol using up to itmax 
c     iterations.  
c
c     The routine is slightly modified from the one presented in:
c        Numerical Recipes The Art of Scientific Computing, 1st Edition 
c        Press, Flannery, Teukolsky, and Vetterling (1986)

  
      integer(ISZ):: ITER
      real(kind=8):: A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
      real(kind=8):: EPS = 3.e-8
      real(kind=8):: watbagconstr

      A=x1
      B=x2
      FA=watbagconstr(A)
      FB=watbagconstr(B)
      if(FB*FA .gt. 0.) then
        call remark('Root must be bracketed for zbrent.') 
        call kaboom(1)  
      endif 
      FC=FB
      do ITER=1,itmax
        if(FB*FC .gt. 0.) then
          C=A
          FC=FA
          D=B-A
          E=D
        endif
        if(abs(FC) .lt. abs(FB)) then
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        endif 
        TOL1=2.*EPS*abs(B)+0.5*tol
        XM=.5*(C-B)
        if(abs(XM) .le. TOL1 .or. FB .eq. 0.) then
          zbrent_test = B
          return
        endif
        if(abs(E) .ge. TOL1 .and. abs(FA) .gt. abs(FB)) then
          S=FB/FA
          if(A .eq. C) then
            P=2.*XM*S
            Q=1.-S
          else
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          endif
          if(P .gt. 0.) Q=-Q
          P=ABS(P)
          if(2.*P .lt. min(3.*XM*Q-abs(TOL1*Q),abs(E*Q))) then
            E=D
            D=P/Q
          else
            D=XM
            E=D
          endif
        else
          D=XM
          E=D
        endif
        A=B
        FA=FB
        if(abs(D) .gt. TOL1) then
          B=B+D
        else
          B=B+sign(TOL1,XM)
        endif
        FB=watbagconstr(B)
      enddo
      zbrent_test = B
      call remark('zbrent exceeding maximum iterations.')
      call kaboom(1)

      return
      end
c=============================================================================
      real(kind=8) function watbagconstr(k1re)

      real(kind=8):: k1re 

c    watbagconstr returns a scaled constraint equation assoicated with 
c    a matched waterbag equilbrium distribution function in a continuous 
c    focusing channel.  
c
c    k1re = k_1 * r_e       k_1 = wavenumber parameter 
c                           r_e = beam edge radius
c

      real(kind=8):: bessi0,bessi
  
      watbagconstr = bessi0(k1re)**2/bessi(2,k1re)**2 - 
     &  4./(k1re)*( 2*bessi(2,k1re)/bessi0(k1re) +
     &              k1re*bessi0(k1re)*bessi(3,k1re)/Bessi(2,k1re)**2
     &            )

      return 
      end 
c=============================================================================
      subroutine rk4(y,dydx,n,x,h,yout,derivs)

      integer(ISZ):: n 
      real(kind=8):: x,h
      real(kind=8):: y(n),dydx(n),yout(n) 
      external derivs 

c     4th order Runge-Kutta advance 
c       x          initial coordinate 
c       y(1:n)     initial y values at x 
c       dydx(1:n)  initial dy/dx derivaties at x 
c       h          step 
c       yout(1:n)  final y values at x + h
c                    y and yout need not be distinct  
c
c       external user supplied subroutine derivs(x,y,dydx) 
c       returns the derivatives dydx evaluated at x, y(x) 
c       where y and dydx have dimension nvar.     
c
c     The routine is slightly modified from the one presented in:
c        Numerical Recipes The Art of Scientific Computing, 1st Edition 
c        Press, Flannery, Teukolsky, and Vetterling (1986)

      integer(ISZ),parameter:: NMAX=10
      real(kind=8):: YT(NMAX),DYT(NMAX),DYM(NMAX)
      integer(ISZ):: I
      real(kind=8):: HH,H6,XH

      HH=H*0.5
      H6=H/6.
      XH=X+HH

c     First step
      do I=1,n
        YT(I)=y(I)+HH*dydx(I)
      enddo
c     Second step
      call derivs(XH,YT,DYT)
      do I=1,N
        YT(I)=Y(I)+HH*DYT(I)
      enddo 
c     Third step 
      call derivs(XH,YT,DYM)
      do I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
      enddo
c     Fourth step 
      call derivs(X+H,YT,DYT)
      do I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
      enddo

      return
      end
c=============================================================================
      subroutine rkdumb(y,x,nvar,nstep,x1,x2,y_work,dy_work,derivs)

      integer(ISZ):: nvar,nstep 
      real(kind=8):: x1,x2
      real(kind=8):: y(nvar,nstep),x(nstep),y_work(nvar),dy_work(nvar)
      external derivs 

c     4th order Runge-Kutta advance with uniform step size  
c
c       external user supplied subroutine derivs(x,y,dydx) 
c       returns the derivatives dydx evaluated at x, y(x) 
c       where y and dydx have dimension nvar.   
c
c     The routine is slightly modified from the one presented in:
c        Numerical Recipes The Art of Scientific Computing, 1st Edition 
c        Press, Flannery, Teukolsky, and Vetterling (1986)

      integer(ISZ):: I,K
      real(kind=8):: H,XX 

c     load starting values 
      do I=1,nvar
        y_work(I)=y(I,1)
      enddo 
      X(1)=X1
      XX=X1
      H=(X2-X1)/float(nstep)
c     take nstep 4th order Runge-Kutta steps using rk4
      do K=1,NSTEP
        call derivs(XX,y_work,dy_work)
        call rk4(y_work,dy_work,nvar,XX,H,y_work,DERIVS)
        if(XX+H .eq. XX) then
          call remark('Stepsize not significant in rkdumb') 
          call kaboom(1) 
        endif 
c       store intermediate steps 
        XX=XX+H
        x(K+1)=XX
        do I=1,NVAR
          y(I,K+1)=y_work(I)
        enddo
      enddo 

      return
      end
c=============================================================================
