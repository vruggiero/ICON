!
!+ Routines for calculating Bessel functions and the root of a function
!
! $Id$
!
MODULE mo_nr
!
! Description:
!   Routines for calculating Bessel functions and the root of a function.
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_2         2008/12/04 Andreas Rhodin
!  remove obsolete PAUSE statement
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Andreas Rhodin  MPIfM  2001  original source
! Harald Anlauf   DWD    2007  minor cleanups
!------------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  use mo_kind,      only: wp      ! working precision kind parameter
  use mo_exception, only: finish  ! abort routine
  implicit none
  !----------------
  ! Public entities
  !----------------
  public :: bessj0      ! Bessel function J_0
  public :: bessj1      ! Bessel function J_1
  public :: zbrent      ! Find zero of a function
  private

contains

  FUNCTION bessj0(x)
  !---------------------------------------------------
  ! Returns the Bessel function J_0(x) for any real X.
  !---------------------------------------------------
  real(wp)             :: bessj0
  real(wp), intent(in) :: x
    real(wp) ax,xx,z
    DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, &
                     &s1,s2,s3,s4,s5,s6,y
    SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4, &
    &s5,s6
    DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,          &
     &-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,  &
     &.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,            &
     &651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2, &
     &s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,           &
     &59272.64853d0,267.8532712d0,1.d0/
    if(abs(x).lt.8.)then
       y=x**2
       bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y* &
            (s4+y*(s5+y*s6)))))
    else
       ax=abs(x)
       z=8._wp/ax
       y=z**2
       xx=ax-.785398164_wp
       bessj0=sqrt(.636619772_wp/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*    &
            p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
    endif
    return
  END FUNCTION bessj0

  FUNCTION bessj1(x)
  !---------------------------------------------------
  ! Returns the Bessel function J_1(x) for any real X.
  !---------------------------------------------------
  real(wp)             :: bessj1
  real(wp), intent(in) :: x
      real(wp) ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, &
     &s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4, &
     &s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,             &
     &242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2, &
     &s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,          &
     &99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,              &
     &.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,    &
     &-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
    if(abs(x).lt.8._wp)then
       y=x**2
       bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+ &
            y*(s4+y*(s5+y*s6)))))
    else
       ax=abs(x)
       z=8._wp/ax
       y=z**2
       xx=ax-2.356194491_wp
       bessj1=sqrt(.636619772_wp/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*    &
            p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1._wp,x)
    endif
    return
  END FUNCTION bessj1

  FUNCTION zbrent(func,x1,x2,tol)
  !------------------------------------------------------------------------
  ! Using Brent's method, find the root of a function FUNC known to lie
  ! between X1 and X2. The root, returned as ZBRENT, will be refined until
  ! its accuracy is TOL.
  !
  ! Parameters: Maximum allowed number of iterations and machine
  ! floating-point precision.
  !------------------------------------------------------------------------
  DOUBLE PRECISION zbrent,tol,x1,x2,func
  EXTERNAL func
    INTEGER ITMAX
    DOUBLE PRECISION EPS
    PARAMETER (ITMAX=100,EPS=3.e-16) !EPS=3.e-8
    INTEGER iter
    DOUBLE PRECISION a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    a=x1
    b=x2
    fa=func(a)
    fb=func(b)
    if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) &
      call finish ('zbrent','root must be bracketed')
    c=b
    fc=fb
    do 11 iter=1,ITMAX
       if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
       endif
       if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2.*EPS*abs(b)+0.5*tol
       xm=.5*(c-b)
       if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
       endif
       if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
             p=2.*xm*s
             q=1.-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
             q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
             e=d
             d=p/q
          else
             d=xm
             e=d
          endif
       else
          d=xm
          e=d
       endif
       a=b
       fa=fb
       if(abs(d) .gt. tol1) then
          b=b+d
       else
          b=b+sign(tol1,xm)
       endif
       fb=func(b)
11  continue
    call finish ('zbrent','exceeding maximum iterations')
    zbrent=b
    return
  END FUNCTION zbrent

end module mo_nr
