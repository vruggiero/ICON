!
!+ Fast Fourier Transform
!
MODULE mo_fft
!
! Description:
!   Fast Fourier Transform
!
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
! V1_8         2009/12/09 Harald Anlauf
!  NEC: replace macro SX by __SX__
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  add private statement, verbosity flag
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
!------------------------------------------------------------------------------

use mo_kind,      only : wp
use mo_exception, only : finish
use mo_mpi_dace,  only : dace
  ! ---------------------------------------------------------------
  !
  ! module *mo_fft* - quantities needed for the fast *fourier transforms.
  !
  ! ---------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: inifft   ! perform factorisation for FFT
  PUBLIC :: fft991cy ! low level FFT routine
  PUBLIC :: trig     ! array passed to fft991cy
  PUBLIC :: nfax     ! array passed to fft991cy
  PUBLIC :: nfft     ! maximum length of calls to FFT
  PUBLIC :: verbose  ! verbosity flag

  REAL(WP), ALLOCATABLE :: trig(:)           ! array used by *fft991*.
  INTEGER               :: nfax(10) = 0      ! array used by *fft991*.
  LOGICAL               :: verbose = .false. ! verbosity flag

  INTEGER,PARAMETER :: nout = 6      ! Fortran unit for stdout

  !  nfft: maximum length of calls to *fft.

#if (defined CRAY) && (! defined _CRAYMPP)
#if (_MAXVL == 64)
  INTEGER, PARAMETER :: nfft = 64
#endif
#if (_MAXVL == 128)
  INTEGER, PARAMETER :: nfft = 128
#endif
#else /* !CRAY */
#if defined (__SX__)
  INTEGER, PARAMETER :: nfft = 256
#else
#if (defined __uxp__)
  INTEGER, PARAMETER :: nfft = 2048
#else
  INTEGER, PARAMETER :: nfft = 16
#endif
#endif
#endif

CONTAINS

  SUBROUTINE inifft(nlon,ier)

    ! Description:
    !
    ! Preset constants in *mo_fft*.
    !
    ! Method:
    !
    ! *inifft* is called from *setdyn*.
    !
    ! Externals:
    ! *set99*      set up trigonometric tables.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, December 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A.Rhodin, DWD May 2007, optional: ier
    !
    ! for more details see file AUTHORS
    !

    IMPLICIT NONE

    INTEGER           :: nlon
    INTEGER ,OPTIONAL :: ier

    !  Executable statements

    if (allocated(trig)) then
      if (size(trig)<nlon) deallocate(trig)
    endif
    IF (.NOT. ALLOCATED(trig)) ALLOCATE(trig(nlon))

  ! Define parameters depending on length of vectorregisters

    if (dace% lpio .and. verbose) then
      WRITE(nout,*) 'Maximum ffts to run in parallel:',nfft
      WRITE(nout,*)
    endif

    !-- 1. Set up trigonometric tables

    CALL set99(trig,nfax,nlon,ier)

    RETURN
  END SUBROUTINE inifft

  SUBROUTINE fft991cy(a,work,trigs,ifax,inc,jump,n,lot,isign)

    ! Description:
    !
    ! Calls fortran-versions of fft s.
    !
    ! Method:
    !
    ! Subroutine 'fft991cy' - multiple fast real periodic transform
    ! supersedes previous routine 'fft991cy'.
    !
    ! Real transform of length n performed by removing redundant
    ! operations from complex transform of length n.
    !
    ! a       is the array containing input & output data.
    ! work    is an area of size (n+1)*min(lot,nfft).
    ! trigs   is a previously prepared list of trig function values.
    ! ifax    is a previously prepared list of factors of n.
    ! inc     is the increment within each data 'vector'
    !         (e.g. inc=1 for consecutively stored data).
    ! jump    is the increment between the start of each data vector.
    ! n       is the length of the data vectors.
    ! lot     is the number of data vectors.
    ! isign = +1 for transform from spectral to gridpoint
    !       = -1 for transform from gridpoint to spectral.
    !
    ! ordering of coefficients:
    ! a(0),b(0),a(1),b(1),a(2),b(2),.,a(n/2),b(n/2)
    ! where b(0)=b(n/2)=0; (n+2) locations required.
    !
    ! ordering of data:
    ! x(0),x(1),x(2),.,x(n-1), 0 , 0 ; (n+2) locations required.
    !
    ! Vectorization is achieved on cray by doing the transforms
    ! in parallel.
    !
    ! n must be composed of factors 2,3 & 5 but does not have to be even.
    !
    ! definition of transforms:
    !
    ! isign=+1: x(j)=sum(k=0,.,n-1)(c(k)*exp(2*i*j*k*pi/n))
    ! where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)

    ! isign=-1: a(k)=(1/n)*sum(j=0,.,n-1)(x(j)*cos(2*j*k*pi/n))
    ! b(k)=-(1/n)*sum(j=0,.,n-1)(x(j)*sin(2*j*k*pi/n))

    ! calls fortran-versions of fft s  !!!
    ! dimension a(n),work(n),trigs(n),ifax(1)

    IMPLICIT NONE

    !  Scalar arguments
    INTEGER :: inc, isign, jump, lot, n

    !  Array arguments
    REAL(WP) :: a(*), trigs(*), work(*)
    INTEGER :: ifax(10)

    !  Local scalars:
    INTEGER :: i, ia, ibase, ierr, ifac, igo, ii, istart, ix, iz, j, jbase, jj, &
         &      k, la, nb, nblox, nfax, nvex, nx

    !  Intrinsic functions
    INTRINSIC MOD


    !  Executable statements

    IF (ifax(10)/=n) CALL set99(trigs,ifax,n)
    nfax = ifax(1)
    nx = n + 1
    IF (MOD(n,2)==1) nx = n
    nblox = 1 + (lot-1)/nfft
    nvex = lot - (nblox-1)*nfft
    IF (isign==-1) GO TO 50

    ! isign=+1, spectral to gridpoint transform

    istart = 1
    DO nb = 1, nblox
       ia = istart
       i = istart
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO j = 1, nvex
          a(i+inc) = 0.5_wp*a(i)
          i = i + jump
       END DO
       IF (MOD(n,2)==1) GO TO 10
       i = istart + n*inc
       DO j = 1, nvex
          a(i) = 0.5_wp*a(i)
          i = i + jump
       END DO
10     CONTINUE
       ia = istart + inc
       la = 1
       igo = + 1

       DO k = 1, nfax
          ifac = ifax(k+1)
          ierr = -1
          IF (igo==-1) GO TO 20
          CALL rpassm(a(ia),a(ia+la*inc),work(1),work(ifac*la+1),trigs,inc,1, &
               &          jump,nx,nvex,n,ifac,la,ierr)
          GO TO 30
20        CONTINUE
          CALL rpassm(work(1),work(la+1),a(ia),a(ia+ifac*la*inc),trigs,1,inc,nx, &
               &          jump,nvex,n,ifac,la,ierr)
30        CONTINUE
          IF (ierr/=0) GO TO 100
          la = ifac*la
          igo = -igo
          ia = istart
       END DO

       ! If necessary, copy results back to a

       IF (MOD(nfax,2)==0) GO TO 40
       ibase = 1
       jbase = ia
       DO jj = 1, nvex
          i = ibase
          j = jbase
          DO ii = 1, n
             a(j) = work(i)
             i = i + 1
             j = j + inc
          END DO
          ibase = ibase + nx
          jbase = jbase + jump
       END DO
40     CONTINUE

       ! Fill in zeros at end

       ix = istart + n*inc
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO j = 1, nvex
          a(ix) = 0.0_wp
          a(ix+inc) = 0.0_wp
          ix = ix + jump
       END DO

       istart = istart + nvex*jump
       nvex = nfft
    END DO
    RETURN

    ! isign=-1, gridpoint to spectral transform

50  CONTINUE
    istart = 1
    DO nb = 1, nblox
       ia = istart
       la = n
       igo = + 1

       DO k = 1, nfax
          ifac = ifax(nfax+2-k)
          la = la/ifac
          ierr = -1
          IF (igo==-1) GO TO 60
          CALL qpassm(a(ia),a(ia+ifac*la*inc),work(1),work(la+1),trigs,inc,1, &
               &          jump,nx,nvex,n,ifac,la,ierr)
          GO TO 70
60        CONTINUE
          CALL qpassm(work(1),work(ifac*la+1),a(ia),a(ia+la*inc),trigs,1,inc,nx, &
               &          jump,nvex,n,ifac,la,ierr)
70        CONTINUE
          IF (ierr/=0) GO TO 100
          igo = -igo
          ia = istart + inc
       END DO

       ! If necessary, copy results back to a

       IF (MOD(nfax,2)==0) GO TO 80
       ibase = 1
       jbase = ia
       DO jj = 1, nvex
          i = ibase
          j = jbase
          DO ii = 1, n
             a(j) = work(i)
             i = i + 1
             j = j + inc
          END DO
          ibase = ibase + nx
          jbase = jbase + jump
       END DO
80     CONTINUE

       ! Shift a(0) & fill in zero imag parts

       ix = istart
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO j = 1, nvex
          a(ix) = a(ix+inc)
          a(ix+inc) = 0.0_wp
          ix = ix + jump
       END DO
       IF (MOD(n,2)==1) GO TO 90
       iz = istart + (n+1)*inc
       DO j = 1, nvex
          a(iz) = 0.0_wp
          iz = iz + jump
       END DO
90     CONTINUE

       istart = istart + nvex*jump
       nvex = nfft
    END DO
    RETURN

    ! Error messages

100 CONTINUE

    if (dace% lpio) then
      SELECT CASE (ierr)
      CASE (:-1)
         WRITE (nout,'(A,I5,A)') ' Vector length =',nvex,', greater than nfft'
      CASE (0)
         WRITE (nout,'(A,I3,A)') ' Factor =',ifac,', not catered for'
      CASE (1:)
         WRITE (nout,'(A,I3,A)') ' Factor =',ifac,&
                                 ', only catered for if la*ifac=n'
      END SELECT
    endif

    RETURN
  END SUBROUTINE fft991cy

  SUBROUTINE qpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la,ierr)

    ! Description:
    !
    ! Performs one pass through data as part of
    ! multiple real fft (fourier analysis) routine.
    !
    ! Method:
    !
    ! a       is first real input vector
    !         equivalence b(1) with a(ifac*la*inc1+1)
    ! c       is first real output vector
    !         equivalence d(1) with c(la*inc2+1)
    ! trigs   is a precalculated list of sines & cosines
    ! inc1    is the addressing increment for a
    ! inc2    is the addressing increment for c
    ! inc3    is the increment between input vectors a
    ! inc4    is the increment between output vectors c
    ! lot     is the number of vectors
    ! n       is the length of the vectors
    ! ifac    is the current factor of n
    !         la = n/(product of factors used so far)
    ! ierr    is an error indicator:
    !         0 - pass completed without error
    !         1 - lot greater than nfft
    !         2 - ifac not catered for
    !         3 - ifac only catered for if la=n/ifac
    !

    IMPLICIT NONE

    !  Scalar arguments
    INTEGER :: ierr, ifac, inc1, inc2, inc3, inc4, la, lot, n

    !  Array arguments
    ! REAL(WP) :: a(n),b(n),c(n),d(n),trigs(n)
    REAL(WP) :: a(*), b(*), c(*), d(*), trigs(*)

    !  Local scalars:
    REAL(WP) :: a0, a1, a10, a11, a2, a20, a21, a3, a4, a5, a6, b0, b1, b10, b11, &
         &      b2, b20, b21, b3, b4, b5, b6, c1, c2, c3, c4, c5, qrt5, s1, s2, s3, s4, &
         &      s5, sin36, sin45, sin60, sin72, z, zqrt5, zsin36, zsin45, zsin60, &
         &      zsin72
    INTEGER :: i, ia, ib, ibad, ibase, ic, id, ie, if, ig, igo, ih, iink, ijk, &
         &      ijump, j, ja, jb, jbase, jc, jd, je, jf, jink, k, kb, kc, kd, ke, kf, &
         &      kstop, l, m

    !  Intrinsic functions
    INTRINSIC REAL, SQRT

    !  Data statements
    DATA sin36/0.587785252292473_wp/, sin72/0.951056516295154_wp/, &
         &      qrt5/0.559016994374947_wp/, sin60/0.866025403784437_wp/


    !  Executable statements

    m = n/ifac
    iink = la*inc1
    jink = la*inc2
    ijump = (ifac-1)*iink
    kstop = (n-ifac)/(2*ifac)

    ibad = 1
    IF (lot>nfft) GO TO 180
    ibase = 0
    jbase = 0
    igo = ifac - 1
    IF (igo==7) igo = 6
    ibad = 2
    IF (igo<1 .OR. igo>6) GO TO 180
    GO TO (10,40,70,100,130,160) igo

    ! Coding for factor 2

10  CONTINUE
    ia = 1
    ib = ia + iink
    ja = 1
    jb = ja + (2*m-la)*inc2

    IF (la==m) GO TO 30

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ia+i) + a(ib+i)
          c(jb+j) = a(ia+i) - a(ib+i)
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    ja = ja + jink
    jink = 2*jink
    jb = jb - jink
    ibase = ibase + ijump
    ijump = 2*ijump + iink
    IF (ja==jb) GO TO 20
    DO k = la, kstop, la
       kb = k + k
       c1 = trigs(kb+1)
       s1 = trigs(kb+2)
       jbase = 0
       DO l = 1, la
          i = ibase
          j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
          DO ijk = 1, lot
             c(ja+j) = a(ia+i) + (c1*a(ib+i)+s1*b(ib+i))
             c(jb+j) = a(ia+i) - (c1*a(ib+i)+s1*b(ib+i))
             d(ja+j) = (c1*b(ib+i)-s1*a(ib+i)) + b(ia+i)
             d(jb+j) = (c1*b(ib+i)-s1*a(ib+i)) - b(ia+i)
             i = i + inc3
             j = j + inc4
          END DO
          ibase = ibase + inc1
          jbase = jbase + inc2
       END DO
       ibase = ibase + ijump
       ja = ja + jink
       jb = jb - jink
    END DO
    IF (ja>jb) GO TO 170
20  CONTINUE
    jbase = 0
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ia+i)
          d(ja+j) = -a(ib+i)
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    GO TO 170
30  CONTINUE
    z = 1.0_wp/REAL(n,wp)
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = z*(a(ia+i)+a(ib+i))
          c(jb+j) = z*(a(ia+i)-a(ib+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    GO TO 170

    ! Coding for factor 3

40  CONTINUE
    ia = 1
    ib = ia + iink
    ic = ib + iink
    ja = 1
    jb = ja + (2*m-la)*inc2
    jc = jb

    IF (la==m) GO TO 60

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ia+i) + (a(ib+i)+a(ic+i))
          c(jb+j) = a(ia+i) - 0.5_wp*(a(ib+i)+a(ic+i))
          d(jb+j) = sin60*(a(ic+i)-a(ib+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    ja = ja + jink
    jink = 2*jink
    jb = jb + jink
    jc = jc - jink
    ibase = ibase + ijump
    ijump = 2*ijump + iink
    IF (ja==jc) GO TO 50
    DO k = la, kstop, la
       kb = k + k
       kc = kb + kb
       c1 = trigs(kb+1)
       s1 = trigs(kb+2)
       c2 = trigs(kc+1)
       s2 = trigs(kc+2)
       jbase = 0
       DO l = 1, la
          i = ibase
          j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
          DO ijk = 1, lot
             a1 = (c1*a(ib+i)+s1*b(ib+i)) + (c2*a(ic+i)+s2*b(ic+i))
             b1 = (c1*b(ib+i)-s1*a(ib+i)) + (c2*b(ic+i)-s2*a(ic+i))
             a2 = a(ia+i) - 0.5_wp*a1
             b2 = b(ia+i) - 0.5_wp*b1
             a3 = sin60*((c1*a(ib+i)+s1*b(ib+i))-(c2*a(ic+i)+s2*b(ic+i)))
             b3 = sin60*((c1*b(ib+i)-s1*a(ib+i))-(c2*b(ic+i)-s2*a(ic+i)))
             c(ja+j) = a(ia+i) + a1
             d(ja+j) = b(ia+i) + b1
             c(jb+j) = a2 + b3
             d(jb+j) = b2 - a3
             c(jc+j) = a2 - b3
             d(jc+j) = -(b2+a3)
             i = i + inc3
             j = j + inc4
          END DO
          ibase = ibase + inc1
          jbase = jbase + inc2
       END DO
       ibase = ibase + ijump
       ja = ja + jink
       jb = jb + jink
       jc = jc - jink
    END DO
    IF (ja>jc) GO TO 170
50  CONTINUE
    jbase = 0
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ia+i) + 0.5_wp*(a(ib+i)-a(ic+i))
          d(ja+j) = -sin60*(a(ib+i)+a(ic+i))
          c(jb+j) = a(ia+i) - (a(ib+i)-a(ic+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    GO TO 170
60  CONTINUE
    z = 1.0_wp/REAL(n,wp)
    zsin60 = z*sin60
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = z*(a(ia+i)+(a(ib+i)+a(ic+i)))
          c(jb+j) = z*(a(ia+i)-0.5_wp*(a(ib+i)+a(ic+i)))
          d(jb+j) = zsin60*(a(ic+i)-a(ib+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    GO TO 170

    ! Coding for factor 4

70  CONTINUE
    ia = 1
    ib = ia + iink
    ic = ib + iink
    id = ic + iink
    ja = 1
    jb = ja + (2*m-la)*inc2
    jc = jb + 2*m*inc2
    jd = jb

    IF (la==m) GO TO 90

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = (a(ia+i)+a(ic+i)) + (a(ib+i)+a(id+i))
          c(jc+j) = (a(ia+i)+a(ic+i)) - (a(ib+i)+a(id+i))
          c(jb+j) = a(ia+i) - a(ic+i)
          d(jb+j) = a(id+i) - a(ib+i)
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    ja = ja + jink
    jink = 2*jink
    jb = jb + jink
    jc = jc - jink
    jd = jd - jink
    ibase = ibase + ijump
    ijump = 2*ijump + iink
    IF (jb==jc) GO TO 80
    DO k = la, kstop, la
       kb = k + k
       kc = kb + kb
       kd = kc + kb
       c1 = trigs(kb+1)
       s1 = trigs(kb+2)
       c2 = trigs(kc+1)
       s2 = trigs(kc+2)
       c3 = trigs(kd+1)
       s3 = trigs(kd+2)
       jbase = 0
       DO l = 1, la
          i = ibase
          j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
          DO ijk = 1, lot
             a0 = a(ia+i) + (c2*a(ic+i)+s2*b(ic+i))
             a2 = a(ia+i) - (c2*a(ic+i)+s2*b(ic+i))
             a1 = (c1*a(ib+i)+s1*b(ib+i)) + (c3*a(id+i)+s3*b(id+i))
             a3 = (c1*a(ib+i)+s1*b(ib+i)) - (c3*a(id+i)+s3*b(id+i))
             b0 = b(ia+i) + (c2*b(ic+i)-s2*a(ic+i))
             b2 = b(ia+i) - (c2*b(ic+i)-s2*a(ic+i))
             b1 = (c1*b(ib+i)-s1*a(ib+i)) + (c3*b(id+i)-s3*a(id+i))
             b3 = (c1*b(ib+i)-s1*a(ib+i)) - (c3*b(id+i)-s3*a(id+i))
             c(ja+j) = a0 + a1
             c(jc+j) = a0 - a1
             d(ja+j) = b0 + b1
             d(jc+j) = b1 - b0
             c(jb+j) = a2 + b3
             c(jd+j) = a2 - b3
             d(jb+j) = b2 - a3
             d(jd+j) = -(b2+a3)
             i = i + inc3
             j = j + inc4
          END DO
          ibase = ibase + inc1
          jbase = jbase + inc2
       END DO
       ibase = ibase + ijump
       ja = ja + jink
       jb = jb + jink
       jc = jc - jink
       jd = jd - jink
    END DO
    IF (jb>jc) GO TO 170
80  CONTINUE
    sin45 = SQRT(0.5_wp)
    jbase = 0
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ia+i) + sin45*(a(ib+i)-a(id+i))
          c(jb+j) = a(ia+i) - sin45*(a(ib+i)-a(id+i))
          d(ja+j) = -a(ic+i) - sin45*(a(ib+i)+a(id+i))
          d(jb+j) = a(ic+i) - sin45*(a(ib+i)+a(id+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    GO TO 170
90  CONTINUE
    z = 1.0_wp/REAL(n,wp)
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = z*((a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i)))
          c(jc+j) = z*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
          c(jb+j) = z*(a(ia+i)-a(ic+i))
          d(jb+j) = z*(a(id+i)-a(ib+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    GO TO 170

    ! Coding for factor 5

100 CONTINUE
    ia = 1
    ib = ia + iink
    ic = ib + iink
    id = ic + iink
    ie = id + iink
    ja = 1
    jb = ja + (2*m-la)*inc2
    jc = jb + 2*m*inc2
    jd = jc
    je = jb

    IF (la==m) GO TO 120

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          a1 = a(ib+i) + a(ie+i)
          a3 = a(ib+i) - a(ie+i)
          a2 = a(ic+i) + a(id+i)
          a4 = a(ic+i) - a(id+i)
          a5 = a(ia+i) - 0.25_wp*(a1+a2)
          a6 = qrt5*(a1-a2)
          c(ja+j) = a(ia+i) + (a1+a2)
          c(jb+j) = a5 + a6
          c(jc+j) = a5 - a6
          d(jb+j) = -sin72*a3 - sin36*a4
          d(jc+j) = -sin36*a3 + sin72*a4
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    ja = ja + jink
    jink = 2*jink
    jb = jb + jink
    jc = jc + jink
    jd = jd - jink
    je = je - jink
    ibase = ibase + ijump
    ijump = 2*ijump + iink
    IF (jb==jd) GO TO 110
    DO k = la, kstop, la
       kb = k + k
       kc = kb + kb
       kd = kc + kb
       ke = kd + kb
       c1 = trigs(kb+1)
       s1 = trigs(kb+2)
       c2 = trigs(kc+1)
       s2 = trigs(kc+2)
       c3 = trigs(kd+1)
       s3 = trigs(kd+2)
       c4 = trigs(ke+1)
       s4 = trigs(ke+2)
       jbase = 0
       DO l = 1, la
          i = ibase
          j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
          DO ijk = 1, lot
             a1 = (c1*a(ib+i)+s1*b(ib+i)) + (c4*a(ie+i)+s4*b(ie+i))
             a3 = (c1*a(ib+i)+s1*b(ib+i)) - (c4*a(ie+i)+s4*b(ie+i))
             a2 = (c2*a(ic+i)+s2*b(ic+i)) + (c3*a(id+i)+s3*b(id+i))
             a4 = (c2*a(ic+i)+s2*b(ic+i)) - (c3*a(id+i)+s3*b(id+i))
             b1 = (c1*b(ib+i)-s1*a(ib+i)) + (c4*b(ie+i)-s4*a(ie+i))
             b3 = (c1*b(ib+i)-s1*a(ib+i)) - (c4*b(ie+i)-s4*a(ie+i))
             b2 = (c2*b(ic+i)-s2*a(ic+i)) + (c3*b(id+i)-s3*a(id+i))
             b4 = (c2*b(ic+i)-s2*a(ic+i)) - (c3*b(id+i)-s3*a(id+i))
             a5 = a(ia+i) - 0.25_wp*(a1+a2)
             a6 = qrt5*(a1-a2)
             b5 = b(ia+i) - 0.25_wp*(b1+b2)
             b6 = qrt5*(b1-b2)
             a10 = a5 + a6
             a20 = a5 - a6
             b10 = b5 + b6
             b20 = b5 - b6
             a11 = sin72*b3 + sin36*b4
             a21 = sin36*b3 - sin72*b4
             b11 = sin72*a3 + sin36*a4
             b21 = sin36*a3 - sin72*a4
             c(ja+j) = a(ia+i) + (a1+a2)
             c(jb+j) = a10 + a11
             c(je+j) = a10 - a11
             c(jc+j) = a20 + a21
             c(jd+j) = a20 - a21
             d(ja+j) = b(ia+i) + (b1+b2)
             d(jb+j) = b10 - b11
             d(je+j) = -(b10+b11)
             d(jc+j) = b20 - b21
             d(jd+j) = -(b20+b21)
             i = i + inc3
             j = j + inc4
          END DO
          ibase = ibase + inc1
          jbase = jbase + inc2
       END DO
       ibase = ibase + ijump
       ja = ja + jink
       jb = jb + jink
       jc = jc + jink
       jd = jd - jink
       je = je - jink
    END DO
    IF (jb>jd) GO TO 170
110 CONTINUE
    jbase = 0
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          a1 = a(ib+i) + a(ie+i)
          a3 = a(ib+i) - a(ie+i)
          a2 = a(ic+i) + a(id+i)
          a4 = a(ic+i) - a(id+i)
          a5 = a(ia+i) + 0.25_wp*(a3-a4)
          a6 = qrt5*(a3+a4)
          c(ja+j) = a5 + a6
          c(jb+j) = a5 - a6
          c(jc+j) = a(ia+i) - (a3-a4)
          d(ja+j) = -sin36*a1 - sin72*a2
          d(jb+j) = -sin72*a1 + sin36*a2
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    GO TO 170
120 CONTINUE
    z = 1.0_wp/REAL(n,wp)
    zqrt5 = z*qrt5
    zsin36 = z*sin36
    zsin72 = z*sin72
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          a1 = a(ib+i) + a(ie+i)
          a3 = a(ib+i) - a(ie+i)
          a2 = a(ic+i) + a(id+i)
          a4 = a(ic+i) - a(id+i)
          a5 = z*(a(ia+i)-0.25_wp*(a1+a2))
          a6 = zqrt5*(a1-a2)
          c(ja+j) = z*(a(ia+i)+(a1+a2))
          c(jb+j) = a5 + a6
          c(jc+j) = a5 - a6
          d(jb+j) = -zsin72*a3 - zsin36*a4
          d(jc+j) = -zsin36*a3 + zsin72*a4
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    GO TO 170

    ! Coding for factor 6

130 CONTINUE
    ia = 1
    ib = ia + iink
    ic = ib + iink
    id = ic + iink
    ie = id + iink
    if = ie + iink
    ja = 1
    jb = ja + (2*m-la)*inc2
    jc = jb + 2*m*inc2
    jd = jc + 2*m*inc2
    je = jc
    jf = jb

    IF (la==m) GO TO 150

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          a11 = (a(ic+i)+a(if+i)) + (a(ib+i)+a(ie+i))
          c(ja+j) = (a(ia+i)+a(id+i)) + a11
          c(jc+j) = (a(ia+i)+a(id+i)-0.5_wp*a11)
          d(jc+j) = sin60*((a(ic+i)+a(if+i))-(a(ib+i)+a(ie+i)))
          a11 = (a(ic+i)-a(if+i)) + (a(ie+i)-a(ib+i))
          c(jb+j) = (a(ia+i)-a(id+i)) - 0.5_wp*a11
          d(jb+j) = sin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(if+i)))
          c(jd+j) = (a(ia+i)-a(id+i)) + a11
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    ja = ja + jink
    jink = 2*jink
    jb = jb + jink
    jc = jc + jink
    jd = jd - jink
    je = je - jink
    jf = jf - jink
    ibase = ibase + ijump
    ijump = 2*ijump + iink
    IF (jc==jd) GO TO 140
    DO k = la, kstop, la
       kb = k + k
       kc = kb + kb
       kd = kc + kb
       ke = kd + kb
       kf = ke + kb
       c1 = trigs(kb+1)
       s1 = trigs(kb+2)
       c2 = trigs(kc+1)
       s2 = trigs(kc+2)
       c3 = trigs(kd+1)
       s3 = trigs(kd+2)
       c4 = trigs(ke+1)
       s4 = trigs(ke+2)
       c5 = trigs(kf+1)
       s5 = trigs(kf+2)
       jbase = 0
       DO l = 1, la
          i = ibase
          j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
          DO ijk = 1, lot
             a1 = c1*a(ib+i) + s1*b(ib+i)
             b1 = c1*b(ib+i) - s1*a(ib+i)
             a2 = c2*a(ic+i) + s2*b(ic+i)
             b2 = c2*b(ic+i) - s2*a(ic+i)
             a3 = c3*a(id+i) + s3*b(id+i)
             b3 = c3*b(id+i) - s3*a(id+i)
             a4 = c4*a(ie+i) + s4*b(ie+i)
             b4 = c4*b(ie+i) - s4*a(ie+i)
             a5 = c5*a(if+i) + s5*b(if+i)
             b5 = c5*b(if+i) - s5*a(if+i)
             a11 = (a2+a5) + (a1+a4)
             a20 = (a(ia+i)+a3) - 0.5_wp*a11
             a21 = sin60*((a2+a5)-(a1+a4))
             b11 = (b2+b5) + (b1+b4)
             b20 = (b(ia+i)+b3) - 0.5_wp*b11
             b21 = sin60*((b2+b5)-(b1+b4))
             c(ja+j) = (a(ia+i)+a3) + a11
             d(ja+j) = (b(ia+i)+b3) + b11
             c(jc+j) = a20 - b21
             d(jc+j) = a21 + b20
             c(je+j) = a20 + b21
             d(je+j) = a21 - b20
             a11 = (a2-a5) + (a4-a1)
             a20 = (a(ia+i)-a3) - 0.5_wp*a11
             a21 = sin60*((a4-a1)-(a2-a5))
             b11 = (b5-b2) - (b4-b1)
             b20 = (b3-b(ia+i)) - 0.5_wp*b11
             b21 = sin60*((b5-b2)+(b4-b1))
             c(jb+j) = a20 - b21
             d(jb+j) = a21 - b20
             c(jd+j) = a11 + (a(ia+i)-a3)
             d(jd+j) = b11 + (b3-b(ia+i))
             c(jf+j) = a20 + b21
             d(jf+j) = a21 + b20
             i = i + inc3
             j = j + inc4
          END DO
          ibase = ibase + inc1
          jbase = jbase + inc2
       END DO
       ibase = ibase + ijump
       ja = ja + jink
       jb = jb + jink
       jc = jc + jink
       jd = jd - jink
       je = je - jink
       jf = jf - jink
    END DO
    IF (jc>jd) GO TO 170
140 CONTINUE
    jbase = 0
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = (a(ia+i)+0.5_wp*(a(ic+i)-a(ie+i))) + sin60*(a(ib+i)-a(if+i))
          d(ja+j) = -(a(id+i)+0.5_wp*(a(ib+i)+a(if+i))) - sin60*(a(ic+i)+a(ie+i))
          c(jb+j) = a(ia+i) - (a(ic+i)-a(ie+i))
          d(jb+j) = a(id+i) - (a(ib+i)+a(if+i))
          c(jc+j) = (a(ia+i)+0.5_wp*(a(ic+i)-a(ie+i))) - sin60*(a(ib+i)-a(if+i))
          d(jc+j) = -(a(id+i)+0.5_wp*(a(ib+i)+a(if+i))) + sin60*(a(ic+i)+a(ie+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    GO TO 170
150 CONTINUE
    z = 1.0_wp/REAL(n,wp)
    zsin60 = z*sin60
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          a11 = (a(ic+i)+a(if+i)) + (a(ib+i)+a(ie+i))
          c(ja+j) = z*((a(ia+i)+a(id+i))+a11)
          c(jc+j) = z*((a(ia+i)+a(id+i))-0.5_wp*a11)
          d(jc+j) = zsin60*((a(ic+i)+a(if+i))-(a(ib+i)+a(ie+i)))
          a11 = (a(ic+i)-a(if+i)) + (a(ie+i)-a(ib+i))
          c(jb+j) = z*((a(ia+i)-a(id+i))-0.5_wp*a11)
          d(jb+j) = zsin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(if+i)))
          c(jd+j) = z*((a(ia+i)-a(id+i))+a11)
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    GO TO 170

    ! Coding for factor 8

160 CONTINUE
    ibad = 3
    IF (la/=m) GO TO 180
    ia = 1
    ib = ia + iink
    ic = ib + iink
    id = ic + iink
    ie = id + iink
    if = ie + iink
    ig = if + iink
    ih = ig + iink
    ja = 1
    jb = ja + la*inc2
    jc = jb + 2*m*inc2
    jd = jc + 2*m*inc2
    je = jd + 2*m*inc2
    z = 1.0_wp/REAL(n,wp)
    zsin45 = z*SQRT(0.5_wp)

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = z*(((a(ia+i)+a(ie+i))+(a(ic+i)+a(ig+i)))+((a(id+i)+ &
               &          a(ih+i))+(a(ib+i)+a(if+i))))
          c(je+j) = z*(((a(ia+i)+a(ie+i))+(a(ic+i)+a(ig+i)))-((a(id+i)+ &
               &          a(ih+i))+(a(ib+i)+a(if+i))))
          c(jc+j) = z*((a(ia+i)+a(ie+i))-(a(ic+i)+a(ig+i)))
          d(jc+j) = z*((a(id+i)+a(ih+i))-(a(ib+i)+a(if+i)))
          c(jb+j) = z*(a(ia+i)-a(ie+i)) + zsin45*((a(ih+i)-a(id+i))-(a(if+ &
               &          i)-a(ib+i)))
          c(jd+j) = z*(a(ia+i)-a(ie+i)) - zsin45*((a(ih+i)-a(id+i))-(a(if+ &
               &          i)-a(ib+i)))
          d(jb+j) = zsin45*((a(ih+i)-a(id+i))+(a(if+i)-a(ib+i))) + &
               &          z*(a(ig+i)-a(ic+i))
          d(jd+j) = zsin45*((a(ih+i)-a(id+i))+(a(if+i)-a(ib+i))) - &
               &          z*(a(ig+i)-a(ic+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    ! Return

170 CONTINUE
    ibad = 0
180 CONTINUE
    ierr = ibad
    RETURN
  END SUBROUTINE qpassm

  SUBROUTINE rpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la,ierr)
    ! Dimension a(n),b(n),c(n),d(n),trigs(n)

    IMPLICIT NONE

    !  Scalar arguments
    INTEGER :: ierr, ifac, inc1, inc2, inc3, inc4, la, lot, n

    !  Array arguments
    REAL(WP) :: a(*), b(*), c(*), d(*), trigs(*)

    !  Local scalars:
    REAL(WP) :: c1, c2, c3, c4, c5, qqrt5, qrt5, s1, s2, s3, s4, s5, sin36, sin45, &
         &      sin60, sin72, ssin36, ssin45, ssin60, ssin72
    INTEGER :: i, ia, ib, ibad, ibase, ic, id, ie, if, igo, iink, ijk, j, ja, &
         &      jb, jbase, jc, jd, je, jf, jg, jh, jink, jump, k, kb, kc, kd, ke, kf, &
         &      kstop, l, m

    !  Local arrays:
    REAL(WP) :: a10(nfft), a11(nfft), a20(nfft), a21(nfft), b10(nfft), b11(nfft), b20(nfft), &
         &      b21(nfft)

    !  Intrinsic functions
    INTRINSIC SQRT

    !  Data statements
    DATA sin36/0.587785252292473_wp/, sin72/0.951056516295154_wp/, &
         &      qrt5/0.559016994374947_wp/, sin60/0.866025403784437_wp/


    !  Executable statements

    m = n/ifac
    iink = la*inc1
    jink = la*inc2
    jump = (ifac-1)*jink
    kstop = (n-ifac)/(2*ifac)

    ibad = 1
    IF (lot>nfft) GO TO 180
    ibase = 0
    jbase = 0
    igo = ifac - 1
    IF (igo==7) igo = 6
    ibad = 2
    IF (igo<1 .OR. igo>6) GO TO 180
    GO TO (10,40,70,100,130,160) igo

    ! Coding for factor 2

10  CONTINUE
    ia = 1
    ib = ia + (2*m-la)*inc1
    ja = 1
    jb = ja + jink

    IF (la==m) GO TO 30

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ia+i) + a(ib+i)
          c(jb+j) = a(ia+i) - a(ib+i)
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    ia = ia + iink
    iink = 2*iink
    ib = ib - iink
    ibase = 0
    jbase = jbase + jump
    jump = 2*jump + jink
    IF (ia==ib) GO TO 20
    DO k = la, kstop, la
       kb = k + k
       c1 = trigs(kb+1)
       s1 = trigs(kb+2)
       ibase = 0
       DO l = 1, la
          i = ibase
          j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
          DO ijk = 1, lot
             c(ja+j) = a(ia+i) + a(ib+i)
             d(ja+j) = b(ia+i) - b(ib+i)
             c(jb+j) = c1*(a(ia+i)-a(ib+i)) - s1*(b(ia+i)+b(ib+i))
             d(jb+j) = s1*(a(ia+i)-a(ib+i)) + c1*(b(ia+i)+b(ib+i))
             i = i + inc3
             j = j + inc4
          END DO
          ibase = ibase + inc1
          jbase = jbase + inc2
       END DO
       ia = ia + iink
       ib = ib - iink
       jbase = jbase + jump
    END DO
    IF (ia>ib) GO TO 170
20  CONTINUE
    ibase = 0
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ia+i)
          c(jb+j) = -b(ia+i)
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    GO TO 170
30  CONTINUE
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = 2.0_wp*(a(ia+i)+a(ib+i))
          c(jb+j) = 2.0_wp*(a(ia+i)-a(ib+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    GO TO 170

    ! Coding for factor 3

40  CONTINUE
    ia = 1
    ib = ia + (2*m-la)*inc1
    ic = ib
    ja = 1
    jb = ja + jink
    jc = jb + jink

    IF (la==m) GO TO 60

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ia+i) + a(ib+i)
          c(jb+j) = (a(ia+i)-0.5_wp*a(ib+i)) - (sin60*(b(ib+i)))
          c(jc+j) = (a(ia+i)-0.5_wp*a(ib+i)) + (sin60*(b(ib+i)))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    ia = ia + iink
    iink = 2*iink
    ib = ib + iink
    ic = ic - iink
    jbase = jbase + jump
    jump = 2*jump + jink
    IF (ia==ic) GO TO 50
    DO k = la, kstop, la
       kb = k + k
       kc = kb + kb
       c1 = trigs(kb+1)
       s1 = trigs(kb+2)
       c2 = trigs(kc+1)
       s2 = trigs(kc+2)
       ibase = 0
       DO l = 1, la
          i = ibase
          j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
          DO ijk = 1, lot
             c(ja+j) = a(ia+i) + (a(ib+i)+a(ic+i))
             d(ja+j) = b(ia+i) + (b(ib+i)-b(ic+i))
             c(jb+j) = c1*((a(ia+i)-0.5_wp*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)+ &
                  &            b(ic+i)))) - s1*((b(ia+i)-0.5_wp*(b(ib+i)-b(ic+i)))+(sin60*(a(ib+i)- &
                  &            a(ic+i))))
             d(jb+j) = s1*((a(ia+i)-0.5_wp*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)+ &
                  &            b(ic+i)))) + c1*((b(ia+i)-0.5_wp*(b(ib+i)-b(ic+i)))+(sin60*(a(ib+i)- &
                  &            a(ic+i))))
             c(jc+j) = c2*((a(ia+i)-0.5_wp*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)+ &
                  &            b(ic+i)))) - s2*((b(ia+i)-0.5_wp*(b(ib+i)-b(ic+i)))-(sin60*(a(ib+i)- &
                  &            a(ic+i))))
             d(jc+j) = s2*((a(ia+i)-0.5_wp*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)+ &
                  &            b(ic+i)))) + c2*((b(ia+i)-0.5_wp*(b(ib+i)-b(ic+i)))-(sin60*(a(ib+i)- &
                  &            a(ic+i))))
             i = i + inc3
             j = j + inc4
          END DO
          ibase = ibase + inc1
          jbase = jbase + inc2
       END DO
       ia = ia + iink
       ib = ib + iink
       ic = ic - iink
       jbase = jbase + jump
    END DO
    IF (ia>ic) GO TO 170
50  CONTINUE
    ibase = 0
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ia+i) + a(ib+i)
          c(jb+j) = (0.5_wp*a(ia+i)-a(ib+i)) - (sin60*b(ia+i))
          c(jc+j) = -(0.5_wp*a(ia+i)-a(ib+i)) - (sin60*b(ia+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    GO TO 170
60  CONTINUE
    ssin60 = 2.0_wp*sin60
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = 2.0_wp*(a(ia+i)+a(ib+i))
          c(jb+j) = (2.0_wp*a(ia+i)-a(ib+i)) - (ssin60*b(ib+i))
          c(jc+j) = (2.0_wp*a(ia+i)-a(ib+i)) + (ssin60*b(ib+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    GO TO 170

    ! Coding for factor 4

70  CONTINUE
    ia = 1
    ib = ia + (2*m-la)*inc1
    ic = ib + 2*m*inc1
    id = ib
    ja = 1
    jb = ja + jink
    jc = jb + jink
    jd = jc + jink

    IF (la==m) GO TO 90

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = (a(ia+i)+a(ic+i)) + a(ib+i)
          c(jb+j) = (a(ia+i)-a(ic+i)) - b(ib+i)
          c(jc+j) = (a(ia+i)+a(ic+i)) - a(ib+i)
          c(jd+j) = (a(ia+i)-a(ic+i)) + b(ib+i)
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    ia = ia + iink
    iink = 2*iink
    ib = ib + iink
    ic = ic - iink
    id = id - iink
    jbase = jbase + jump
    jump = 2*jump + jink
    IF (ib==ic) GO TO 80
    DO k = la, kstop, la
       kb = k + k
       kc = kb + kb
       kd = kc + kb
       c1 = trigs(kb+1)
       s1 = trigs(kb+2)
       c2 = trigs(kc+1)
       s2 = trigs(kc+2)
       c3 = trigs(kd+1)
       s3 = trigs(kd+2)
       ibase = 0
       DO l = 1, la
          i = ibase
          j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
          DO ijk = 1, lot
             c(ja+j) = (a(ia+i)+a(ic+i)) + (a(ib+i)+a(id+i))
             d(ja+j) = (b(ia+i)-b(ic+i)) + (b(ib+i)-b(id+i))
             c(jc+j) = c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) - s2*((b(ia+ &
                  &            i)-b(ic+i))-(b(ib+i)-b(id+i)))
             d(jc+j) = s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) + c2*((b(ia+ &
                  &            i)-b(ic+i))-(b(ib+i)-b(id+i)))
             c(jb+j) = c1*((a(ia+i)-a(ic+i))-(b(ib+i)+b(id+i))) - s1*((b(ia+ &
                  &            i)+b(ic+i))+(a(ib+i)-a(id+i)))
             d(jb+j) = s1*((a(ia+i)-a(ic+i))-(b(ib+i)+b(id+i))) + c1*((b(ia+ &
                  &            i)+b(ic+i))+(a(ib+i)-a(id+i)))
             c(jd+j) = c3*((a(ia+i)-a(ic+i))+(b(ib+i)+b(id+i))) - s3*((b(ia+ &
                  &            i)+b(ic+i))-(a(ib+i)-a(id+i)))
             d(jd+j) = s3*((a(ia+i)-a(ic+i))+(b(ib+i)+b(id+i))) + c3*((b(ia+ &
                  &            i)+b(ic+i))-(a(ib+i)-a(id+i)))
             i = i + inc3
             j = j + inc4
          END DO
          ibase = ibase + inc1
          jbase = jbase + inc2
       END DO
       ia = ia + iink
       ib = ib + iink
       ic = ic - iink
       id = id - iink
       jbase = jbase + jump
    END DO
    IF (ib>ic) GO TO 170
80  CONTINUE
    ibase = 0
    sin45 = SQRT(0.5_wp)
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ia+i) + a(ib+i)
          c(jb+j) = sin45*((a(ia+i)-a(ib+i))-(b(ia+i)+b(ib+i)))
          c(jc+j) = b(ib+i) - b(ia+i)
          c(jd+j) = -sin45*((a(ia+i)-a(ib+i))+(b(ia+i)+b(ib+i)))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    GO TO 170
90  CONTINUE
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = 2.0_wp*((a(ia+i)+a(ic+i))+a(ib+i))
          c(jb+j) = 2.0_wp*((a(ia+i)-a(ic+i))-b(ib+i))
          c(jc+j) = 2.0_wp*((a(ia+i)+a(ic+i))-a(ib+i))
          c(jd+j) = 2.0_wp*((a(ia+i)-a(ic+i))+b(ib+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    ! Coding for factor 5

    GO TO 170
100 CONTINUE
    ia = 1
    ib = ia + (2*m-la)*inc1
    ic = ib + 2*m*inc1
    id = ic
    ie = ib
    ja = 1
    jb = ja + jink
    jc = jb + jink
    jd = jc + jink
    je = jd + jink

    IF (la==m) GO TO 120

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ia+i) + (a(ib+i)+a(ic+i))
          c(jb+j) = ((a(ia+i)-0.25_wp*(a(ib+i)+a(ic+i)))+qrt5*(a(ib+i)-a(ic+i))) - &
               &          (sin72*b(ib+i)+sin36*b(ic+i))
          c(jc+j) = ((a(ia+i)-0.25_wp*(a(ib+i)+a(ic+i)))-qrt5*(a(ib+i)-a(ic+i))) - &
               &          (sin36*b(ib+i)-sin72*b(ic+i))
          c(jd+j) = ((a(ia+i)-0.25_wp*(a(ib+i)+a(ic+i)))-qrt5*(a(ib+i)-a(ic+i))) + &
               &          (sin36*b(ib+i)-sin72*b(ic+i))
          c(je+j) = ((a(ia+i)-0.25_wp*(a(ib+i)+a(ic+i)))+qrt5*(a(ib+i)-a(ic+i))) + &
               &          (sin72*b(ib+i)+sin36*b(ic+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    ia = ia + iink
    iink = 2*iink
    ib = ib + iink
    ic = ic + iink
    id = id - iink
    ie = ie - iink
    jbase = jbase + jump
    jump = 2*jump + jink
    IF (ib==id) GO TO 110
    DO k = la, kstop, la
       kb = k + k
       kc = kb + kb
       kd = kc + kb
       ke = kd + kb
       c1 = trigs(kb+1)
       s1 = trigs(kb+2)
       c2 = trigs(kc+1)
       s2 = trigs(kc+2)
       c3 = trigs(kd+1)
       s3 = trigs(kd+2)
       c4 = trigs(ke+1)
       s4 = trigs(ke+2)
       ibase = 0
       DO l = 1, la
          i = ibase
          j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
          DO ijk = 1, lot

             a10(ijk) = (a(ia+i)-0.25_wp*((a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i)))) + &
                  &            qrt5*((a(ib+i)+a(ie+i))-(a(ic+i)+a(id+i)))
             a20(ijk) = (a(ia+i)-0.25_wp*((a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i)))) - &
                  &            qrt5*((a(ib+i)+a(ie+i))-(a(ic+i)+a(id+i)))
             b10(ijk) = (b(ia+i)-0.25_wp*((b(ib+i)-b(ie+i))+(b(ic+i)-b(id+i)))) + &
                  &            qrt5*((b(ib+i)-b(ie+i))-(b(ic+i)-b(id+i)))
             b20(ijk) = (b(ia+i)-0.25_wp*((b(ib+i)-b(ie+i))+(b(ic+i)-b(id+i)))) - &
                  &            qrt5*((b(ib+i)-b(ie+i))-(b(ic+i)-b(id+i)))
             a11(ijk) = sin72*(b(ib+i)+b(ie+i)) + sin36*(b(ic+i)+b(id+i))
             a21(ijk) = sin36*(b(ib+i)+b(ie+i)) - sin72*(b(ic+i)+b(id+i))
             b11(ijk) = sin72*(a(ib+i)-a(ie+i)) + sin36*(a(ic+i)-a(id+i))
             b21(ijk) = sin36*(a(ib+i)-a(ie+i)) - sin72*(a(ic+i)-a(id+i))

             c(ja+j) = a(ia+i) + ((a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i)))
             d(ja+j) = b(ia+i) + ((b(ib+i)-b(ie+i))+(b(ic+i)-b(id+i)))
             c(jb+j) = c1*(a10(ijk)-a11(ijk)) - s1*(b10(ijk)+b11(ijk))
             d(jb+j) = s1*(a10(ijk)-a11(ijk)) + c1*(b10(ijk)+b11(ijk))
             c(je+j) = c4*(a10(ijk)+a11(ijk)) - s4*(b10(ijk)-b11(ijk))
             d(je+j) = s4*(a10(ijk)+a11(ijk)) + c4*(b10(ijk)-b11(ijk))
             c(jc+j) = c2*(a20(ijk)-a21(ijk)) - s2*(b20(ijk)+b21(ijk))
             d(jc+j) = s2*(a20(ijk)-a21(ijk)) + c2*(b20(ijk)+b21(ijk))
             c(jd+j) = c3*(a20(ijk)+a21(ijk)) - s3*(b20(ijk)-b21(ijk))
             d(jd+j) = s3*(a20(ijk)+a21(ijk)) + c3*(b20(ijk)-b21(ijk))

             i = i + inc3
             j = j + inc4
          END DO
          ibase = ibase + inc1
          jbase = jbase + inc2
       END DO
       ia = ia + iink
       ib = ib + iink
       ic = ic + iink
       id = id - iink
       ie = ie - iink
       jbase = jbase + jump
    END DO
    IF (ib>id) GO TO 170
110 CONTINUE
    ibase = 0
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = (a(ia+i)+a(ib+i)) + a(ic+i)
          c(jb+j) = (qrt5*(a(ia+i)-a(ib+i))+(0.25_wp*(a(ia+i)+a(ib+i))-a(ic+i))) - &
               &          (sin36*b(ia+i)+sin72*b(ib+i))
          c(je+j) = -(qrt5*(a(ia+i)-a(ib+i))+(0.25_wp*(a(ia+i)+a(ib+i))-a(ic+i))) - &
               &          (sin36*b(ia+i)+sin72*b(ib+i))
          c(jc+j) = (qrt5*(a(ia+i)-a(ib+i))-(0.25_wp*(a(ia+i)+a(ib+i))-a(ic+i))) - &
               &          (sin72*b(ia+i)-sin36*b(ib+i))
          c(jd+j) = -(qrt5*(a(ia+i)-a(ib+i))-(0.25_wp*(a(ia+i)+a(ib+i))-a(ic+i))) - &
               &          (sin72*b(ia+i)-sin36*b(ib+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    GO TO 170
120 CONTINUE
    qqrt5 = 2.0_wp*qrt5
    ssin36 = 2.0_wp*sin36
    ssin72 = 2.0_wp*sin72
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = 2.0_wp*(a(ia+i)+(a(ib+i)+a(ic+i)))
          c(jb+j) = (2.0_wp*(a(ia+i)-0.25_wp*(a(ib+i)+a(ic+i)))+qqrt5*(a(ib+i)-a(ic+ &
               &          i))) - (ssin72*b(ib+i)+ssin36*b(ic+i))
          c(jc+j) = (2.0_wp*(a(ia+i)-0.25_wp*(a(ib+i)+a(ic+i)))-qqrt5*(a(ib+i)-a(ic+ &
               &          i))) - (ssin36*b(ib+i)-ssin72*b(ic+i))
          c(jd+j) = (2.0_wp*(a(ia+i)-0.25_wp*(a(ib+i)+a(ic+i)))-qqrt5*(a(ib+i)-a(ic+ &
               &          i))) + (ssin36*b(ib+i)-ssin72*b(ic+i))
          c(je+j) = (2.0_wp*(a(ia+i)-0.25_wp*(a(ib+i)+a(ic+i)))+qqrt5*(a(ib+i)-a(ic+ &
               &          i))) + (ssin72*b(ib+i)+ssin36*b(ic+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    GO TO 170

    ! Coding for factor 6

130 CONTINUE
    ia = 1
    ib = ia + (2*m-la)*inc1
    ic = ib + 2*m*inc1
    id = ic + 2*m*inc1
    ie = ic
    if = ib
    ja = 1
    jb = ja + jink
    jc = jb + jink
    jd = jc + jink
    je = jd + jink
    jf = je + jink

    IF (la==m) GO TO 150

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = (a(ia+i)+a(id+i)) + (a(ib+i)+a(ic+i))
          c(jd+j) = (a(ia+i)-a(id+i)) - (a(ib+i)-a(ic+i))
          c(jb+j) = ((a(ia+i)-a(id+i))+0.5_wp*(a(ib+i)-a(ic+i))) - (sin60*(b(ib+ &
               &          i)+b(ic+i)))
          c(jf+j) = ((a(ia+i)-a(id+i))+0.5_wp*(a(ib+i)-a(ic+i))) + (sin60*(b(ib+ &
               &          i)+b(ic+i)))
          c(jc+j) = ((a(ia+i)+a(id+i))-0.5_wp*(a(ib+i)+a(ic+i))) - (sin60*(b(ib+ &
               &          i)-b(ic+i)))
          c(je+j) = ((a(ia+i)+a(id+i))-0.5_wp*(a(ib+i)+a(ic+i))) + (sin60*(b(ib+ &
               &          i)-b(ic+i)))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    ia = ia + iink
    iink = 2*iink
    ib = ib + iink
    ic = ic + iink
    id = id - iink
    ie = ie - iink
    if = if - iink
    jbase = jbase + jump
    jump = 2*jump + jink
    IF (ic==id) GO TO 140
    DO k = la, kstop, la
       kb = k + k
       kc = kb + kb
       kd = kc + kb
       ke = kd + kb
       kf = ke + kb
       c1 = trigs(kb+1)
       s1 = trigs(kb+2)
       c2 = trigs(kc+1)
       s2 = trigs(kc+2)
       c3 = trigs(kd+1)
       s3 = trigs(kd+2)
       c4 = trigs(ke+1)
       s4 = trigs(ke+2)
       c5 = trigs(kf+1)
       s5 = trigs(kf+2)
       ibase = 0
       DO l = 1, la
          i = ibase
          j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
          DO ijk = 1, lot

             a11(ijk) = (a(ie+i)+a(ib+i)) + (a(ic+i)+a(if+i))
             a20(ijk) = (a(ia+i)+a(id+i)) - 0.5_wp*a11(ijk)
             a21(ijk) = sin60*((a(ie+i)+a(ib+i))-(a(ic+i)+a(if+i)))
             b11(ijk) = (b(ib+i)-b(ie+i)) + (b(ic+i)-b(if+i))
             b20(ijk) = (b(ia+i)-b(id+i)) - 0.5_wp*b11(ijk)
             b21(ijk) = sin60*((b(ib+i)-b(ie+i))-(b(ic+i)-b(if+i)))

             c(ja+j) = (a(ia+i)+a(id+i)) + a11(ijk)
             d(ja+j) = (b(ia+i)-b(id+i)) + b11(ijk)
             c(jc+j) = c2*(a20(ijk)-b21(ijk)) - s2*(b20(ijk)+a21(ijk))
             d(jc+j) = s2*(a20(ijk)-b21(ijk)) + c2*(b20(ijk)+a21(ijk))
             c(je+j) = c4*(a20(ijk)+b21(ijk)) - s4*(b20(ijk)-a21(ijk))
             d(je+j) = s4*(a20(ijk)+b21(ijk)) + c4*(b20(ijk)-a21(ijk))

             a11(ijk) = (a(ie+i)-a(ib+i)) + (a(ic+i)-a(if+i))
             b11(ijk) = (b(ie+i)+b(ib+i)) - (b(ic+i)+b(if+i))
             a20(ijk) = (a(ia+i)-a(id+i)) - 0.5_wp*a11(ijk)
             a21(ijk) = sin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(if+i)))
             b20(ijk) = (b(ia+i)+b(id+i)) + 0.5_wp*b11(ijk)
             b21(ijk) = sin60*((b(ie+i)+b(ib+i))+(b(ic+i)+b(if+i)))

             c(jd+j) = c3*((a(ia+i)-a(id+i))+a11(ijk)) - s3*((b(ia+i)+b(id+ &
                  &            i))-b11(ijk))
             d(jd+j) = s3*((a(ia+i)-a(id+i))+a11(ijk)) + c3*((b(ia+i)+b(id+ &
                  &            i))-b11(ijk))
             c(jb+j) = c1*(a20(ijk)-b21(ijk)) - s1*(b20(ijk)-a21(ijk))
             d(jb+j) = s1*(a20(ijk)-b21(ijk)) + c1*(b20(ijk)-a21(ijk))
             c(jf+j) = c5*(a20(ijk)+b21(ijk)) - s5*(b20(ijk)+a21(ijk))
             d(jf+j) = s5*(a20(ijk)+b21(ijk)) + c5*(b20(ijk)+a21(ijk))

             i = i + inc3
             j = j + inc4
          END DO
          ibase = ibase + inc1
          jbase = jbase + inc2
       END DO
       ia = ia + iink
       ib = ib + iink
       ic = ic + iink
       id = id - iink
       ie = ie - iink
       if = if - iink
       jbase = jbase + jump
    END DO
    IF (ic>id) GO TO 170
140 CONTINUE
    ibase = 0
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = a(ib+i) + (a(ia+i)+a(ic+i))
          c(jd+j) = b(ib+i) - (b(ia+i)+b(ic+i))
          c(jb+j) = (sin60*(a(ia+i)-a(ic+i))) - (0.5_wp*(b(ia+i)+b(ic+i))+b(ib+i))
          c(jf+j) = -(sin60*(a(ia+i)-a(ic+i))) - (0.5_wp*(b(ia+i)+b(ic+i))+b(ib+i))
          c(jc+j) = sin60*(b(ic+i)-b(ia+i)) + (0.5_wp*(a(ia+i)+a(ic+i))-a(ib+i))
          c(je+j) = sin60*(b(ic+i)-b(ia+i)) - (0.5_wp*(a(ia+i)+a(ic+i))-a(ib+i))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    GO TO 170
150 CONTINUE
    ssin60 = 2.0_wp*sin60
    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = (2.0_wp*(a(ia+i)+a(id+i))) + (2.0_wp*(a(ib+i)+a(ic+i)))
          c(jd+j) = (2.0_wp*(a(ia+i)-a(id+i))) - (2.0_wp*(a(ib+i)-a(ic+i)))
          c(jb+j) = (2.0_wp*(a(ia+i)-a(id+i))+(a(ib+i)-a(ic+i))) - (ssin60*(b(ib+ &
               &          i)+b(ic+i)))
          c(jf+j) = (2.0_wp*(a(ia+i)-a(id+i))+(a(ib+i)-a(ic+i))) + (ssin60*(b(ib+ &
               &          i)+b(ic+i)))
          c(jc+j) = (2.0_wp*(a(ia+i)+a(id+i))-(a(ib+i)+a(ic+i))) - (ssin60*(b(ib+ &
               &          i)-b(ic+i)))
          c(je+j) = (2.0_wp*(a(ia+i)+a(id+i))-(a(ib+i)+a(ic+i))) + (ssin60*(b(ib+ &
               &          i)-b(ic+i)))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO
    GO TO 170

    ! Coding for factor 8

160 CONTINUE
    ibad = 3
    IF (la/=m) GO TO 180
    ia = 1
    ib = ia + la*inc1
    ic = ib + 2*la*inc1
    id = ic + 2*la*inc1
    ie = id + 2*la*inc1
    ja = 1
    jb = ja + jink
    jc = jb + jink
    jd = jc + jink
    je = jd + jink
    jf = je + jink
    jg = jf + jink
    jh = jg + jink
    ssin45 = SQRT(2.0_wp)

    DO l = 1, la
       i = ibase
       j = jbase
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
       DO ijk = 1, lot
          c(ja+j) = 2.0_wp*(((a(ia+i)+a(ie+i))+a(ic+i))+(a(ib+i)+a(id+i)))
          c(je+j) = 2.0_wp*(((a(ia+i)+a(ie+i))+a(ic+i))-(a(ib+i)+a(id+i)))
          c(jc+j) = 2.0_wp*(((a(ia+i)+a(ie+i))-a(ic+i))-(b(ib+i)-b(id+i)))
          c(jg+j) = 2.0_wp*(((a(ia+i)+a(ie+i))-a(ic+i))+(b(ib+i)-b(id+i)))
          c(jb+j) = 2.0_wp*((a(ia+i)-a(ie+i))-b(ic+i)) + ssin45*((a(ib+i)-a(id+ &
               &          i))-(b(ib+i)+b(id+i)))
          c(jf+j) = 2.0_wp*((a(ia+i)-a(ie+i))-b(ic+i)) - ssin45*((a(ib+i)-a(id+ &
               &          i))-(b(ib+i)+b(id+i)))
          c(jd+j) = 2.0_wp*((a(ia+i)-a(ie+i))+b(ic+i)) - ssin45*((a(ib+i)-a(id+ &
               &          i))+(b(ib+i)+b(id+i)))
          c(jh+j) = 2.0_wp*((a(ia+i)-a(ie+i))+b(ic+i)) + ssin45*((a(ib+i)-a(id+ &
               &          i))+(b(ib+i)+b(id+i)))
          i = i + inc3
          j = j + inc4
       END DO
       ibase = ibase + inc1
       jbase = jbase + inc2
    END DO

    ! Return

170 CONTINUE
    ibad = 0
180 CONTINUE
    ierr = ibad
    RETURN
  END SUBROUTINE rpassm

  SUBROUTINE set99(trigs,ifax,n,ier)

    ! Description:
    !
    ! Computes factors of n & trigonometric functins required by fft99 & fft991cy
    !
    ! Method:
    !
    ! Dimension trigs(n),ifax(1),jfax(10),lfax(7)
    !
    ! subroutine 'set99' - computes factors of n & trigonometric
    ! functins required by fft99 & fft991cy

    IMPLICIT NONE

    !  Scalar arguments
    INTEGER           :: n
    INTEGER ,OPTIONAL :: ier

    !  Array arguments
    REAL(WP) :: trigs(*)
    INTEGER :: ifax(10)

    !  Local scalars:
    REAL(WP) :: angle, del
    INTEGER :: i, ifac, k, l, nfax, nhl, nil, nu

    !  Local arrays:
    INTEGER :: jfax(10), lfax(7)

    !  Intrinsic functions
    INTRINSIC ASIN, COS, MOD, REAL, SIN

    !  Data statements
    DATA lfax/6, 8, 5, 4, 3, 2, 1/


    !  Executable statements
    if (dace% lpio .and. verbose) then
      print *,'set99:',n
      print *,'ifax:',ifax
    endif

    del = 4.0_wp*ASIN(1.0_wp)/REAL(n,wp)
    nil = 0
    nhl = (n/2) - 1
    DO k = nil, nhl
       angle = REAL(k,wp)*del
       trigs(2*k+1) = COS(angle)
       trigs(2*k+2) = SIN(angle)
    END DO

    ! Find factors of n (8,6,5,4,3,2; only one 8 allowed)
    ! Look for sixes first, store factors in descending order
    nu = n
    ifac = 6
    k = 0
    l = 1
10  CONTINUE
    IF (MOD(nu,ifac)/=0) GO TO 30
    k = k + 1
    jfax(k) = ifac
    IF (ifac/=8) GO TO 20
    IF (k==1) GO TO 20
    jfax(1) = 8
    jfax(k) = 6
20  CONTINUE
    nu = nu/ifac
    IF (nu==1) GO TO 40
    IF (ifac/=8) GO TO 10
30  CONTINUE
    l = l + 1
    ifac = lfax(l)
    IF (ifac>1) GO TO 10

    if (dace% lpio) &
      WRITE (nout,'(A,I4,A)') ' n =',n,' - Contains illegal factors'
    IF(PRESENT(ier)) THEN
      ier = 1
    ELSE
      CALL finish('set99','n Contains illegal factors')
    ENDIF

    RETURN

    ! Now reverse order of factors
40  CONTINUE
    nfax = k
    ifax(1) = nfax
    DO i = 1, nfax
       ifax(nfax+2-i) = jfax(i)
    END DO
    ifax(10) = n
    IF (dace% lpio .and. verbose) &
      PRINT *,'ifax:',ifax
    IF(PRESENT(ier)) ier = 0
    RETURN
  END SUBROUTINE set99

END MODULE mo_fft
