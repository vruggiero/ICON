!
!+ Legendre transformation
!
! $Id$
!
MODULE mo_legendre
!
! Description:
!   This module gathers all routines and coefficents from the ECHAM model
!   required for performing Legendre transformations.
!   Usage:
!    1)'call inileg (im,in,ik,igl)' to specify truncation and grid resolution.
!    2)'call gp2sp (gp,sp)' or 'call sp2gp (sp,gp)' to perform transformations.
!    3) goto 1) if different truncation/resolution is required.
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
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_14        2011/11/08 Andreas Rhodin
!  fix comment line
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
!==============================================================================

  use mo_kind, only : wp
  IMPLICIT NONE

  !================
  ! Public entities
  !================
  private
  !------------------------------------------------------------------
  ! Quantities depending on the current truncation/spatial resolution
  !   set by subroutine inileg
  !------------------------------------------------------------------
  public :: gmu     ! sin(gaussian latitudes)
  public :: gw      ! gaussian weights.
  public :: pnm     ! Legendre coefficients (scalar transform)
  public :: anm     ! Legendre coefficients (meridional derivatives)
  public :: nmp     ! displacem.of the first point of columns in spectral space
  public :: nnp     ! number of points on each column
  public :: nm,nn,nk! Truncation
  public :: nsp     ! number of spectral coefficients
  public :: ngl     ! number of gaussian latitudes
  public :: nmp1    ! max zonal wave number              + 1
  public :: nnp1    ! max meridional wave number for m=0 + 1
  !------------------------------------------------
  ! Quantities depending on the truncation/spatial,
  !   specific to the latitude currently processed
  ! ..set by legmod
  !------------------------------------------------
  public :: pnmd    ! Modified coefficients for direct Legendre transform
  public :: anmd    ! ...
  public :: rnmd    !
  !----------------
  ! ..set by leginv
  !----------------
  public :: pnmi    ! Modified coefficients for inverse Legendre transform
  public :: anmi    ! ...
  public :: pnmiuv  !
  public :: anmiuv  !
  !------------------
  ! module procedures
  !------------------
  public :: inileg  ! Set up polynomials needed for the Legendre transforms.
  public :: legmod  ! Calculate modified Legendre polynomials for direct tran.
  public :: leginv  ! Calculate modified Legendre polynomials for inverse tra.
  public :: sp2gp   ! Legendre transform spherical harmonics -> gridpoint
  public :: gp2sp   ! Legendre transform gridpoint -> spherical harmonics
  public :: gauaw   ! Calculate gaussian latitudes and weights
  public :: phcs    ! Computes Legendre polynomials and meridional derivatives
  public :: reord   ! Reorders Legendre polynomials for pentagonal truncation
!==============================================================================
  !======================
  ! Variable declarations
  !======================
  !----------
  ! constants
  !----------
  REAL(WP) ,PARAMETER   :: api = 3.14159265358979323846_wp
  REAL(WP) ,PARAMETER   :: a   = 6371000.0_wp
  !--------------------------------
  ! specification of the truncation
  !--------------------------------
  INTEGER               :: nm = -1  ! max zonal wave number
  INTEGER               :: nn       ! max meridional wave number for m=0
  INTEGER               :: nk       ! max meridional wave number
  !---------------------------------
  ! specification of grid resolution
  !---------------------------------
  INTEGER               :: ngl      ! number of gaussian latitudes
  REAL(WP) ,ALLOCATABLE :: gmu(:)   ! sin(gaussian latitudes)
  REAL(WP) ,ALLOCATABLE :: gw (:)   ! gaussian weights.
  REAL(WP) ,ALLOCATABLE :: cst(:)   ! square of cos(latitude):(1-mu**2)
  !------------
  ! coefficents
  !------------
  REAL(WP) ,ALLOCATABLE :: pnm(:,:) ! Legendre polynominals
  REAL(WP) ,ALLOCATABLE :: anm(:,:)
  INTEGER  ,ALLOCATABLE :: nmp(:)   ! displacement of the first column entry
  INTEGER  ,ALLOCATABLE :: nnp(:)   ! number of points on each column
  INTEGER               :: nhgl     ! (number of gaussian latitudes)/2
  INTEGER               :: nsp      ! number of spectral coefficients
  !-------------------
  ! derived quantities
  !-------------------
  INTEGER               :: nmp1     ! max zonal wave number              + 1
  INTEGER               :: nnp1     ! max meridional wave number for m=0 + 1
  INTEGER               :: nkp1     ! max meridional wave number         + 1
  !--------------------------------------------------
  ! Modified Legendre coefficients for one latitudude
  !   set by subroutine legmod (module mo_legendre)
  !--------------------------------------------------
  REAL(WP), TARGET,  ALLOCATABLE :: pnmd(:)   ! direct Legendre transform
  REAL(WP), TARGET,  ALLOCATABLE :: anmd(:)
  REAL(WP), TARGET,  ALLOCATABLE :: rnmd(:)
  !------------------------------------------------
  !   set by subroutine leginv (module mo_legendre)
  !------------------------------------------------
  REAL(WP), TARGET,  ALLOCATABLE :: pnmi(:)   ! inverse Legendre transform
  REAL(WP), TARGET,  ALLOCATABLE :: anmi(:)
  REAL(WP), TARGET,  ALLOCATABLE :: pnmiuv(:)
  REAL(WP), TARGET,  ALLOCATABLE :: anmiuv(:)
!==============================================================================
CONTAINS
!==============================================================================
  SUBROUTINE gauaw (pa, pw, nlat)

    ! Description:
    !
    ! Compute abscissas and weights for gaussian integration.
    !
    ! Method:
    !

    IMPLICIT NONE

    !  Scalar arguments
    INTEGER,  intent(in)  :: nlat

    !  Array arguments
    REAL(WP), intent(out) :: pa(nlat), pw(nlat)
    ! *pa*  - array, length at least *k,* to receive abscissas.
    ! *pw*  - array, length at least *k,* to receive weights.


    !  Local scalars:
    REAL(WP), PARAMETER :: eps = EPSILON(0.0_wp)
    INTEGER,  PARAMETER :: itemax = 20

    INTEGER :: iter, ins2, isym, jn, jgl
    REAL(WP)    :: za, zw, z, zan
    REAL(WP)    :: zk, zkm1, zkm2, zx, zxn, zldn, zmod

    !  Intrinsic functions
    INTRINSIC ABS, COS, MOD, TAN


    !  Executable statements

    ins2 = nlat/2+MOD(nlat,2)

    ! Find first approximation of the roots of the
    ! Legendre polynomial of degree nlat

    DO jgl = 1, ins2
       z = REAL(4*jgl-1,wp)*api/REAL(4*nlat+2,wp)
       pa(jgl) = COS(z+1._wp/(TAN(z)*REAL(8*nlat**2,wp)))
    END DO

    ! Computes roots and weights
    ! Perform the Newton loop
    ! Find 0 of Legendre polynomial with Newton loop

    DO jgl = 1, ins2

       za = pa(jgl)

       DO iter = 1, itemax+1
          zk = 0.0_wp

          ! Newton iteration step

          zkm2 = 1.0_wp
          zkm1 = za
          zx = za
          DO jn = 2, nlat
             zk = (REAL(2*jn-1,wp)*zx*zkm1-REAL(jn-1,wp)*zkm2)/REAL(jn,wp)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zldn = (REAL(nlat,wp)*(zkm1-zx*zk))/(1._wp-zx*zx)
          zmod = -zk/zldn
          zxn = zx+zmod
          zan = zxn

          ! computes weight

          zkm2 = 1.0_wp
          zkm1 = zxn
          zx = zxn
          DO jn = 2,nlat
             zk = (REAL(2*jn-1,wp)*zx*zkm1-REAL(jn-1,wp)*zkm2)/REAL(jn,wp)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zw = (1.0_wp-zx*zx)/(REAL(nlat*nlat,wp)*zkm1*zkm1)
          za = zan
          IF (ABS(zmod) <= eps) EXIT
       END DO

       pa(jgl) = zan
       pw(jgl) = zw

    ENDDO

!DIR$ IVDEP
!OCL NOVREC

    DO jgl = 1, nlat/2
       isym = nlat-jgl+1
       pa(isym) = -pa(jgl)
       pw(isym) = pw(jgl)
    ENDDO

  END SUBROUTINE gauaw
!==============================================================================
     SUBROUTINE phcs (ppnm, panm, kmp1, kkp1, pmu)

       ! Description:
       !
       ! *phcs* computes the values of the *Legendre polynomials and of their
       ! meridional derivatives at given latitudes for a rhomboidal truncation.
       !
       !
       !
       !        ^
       !   *N*  !
       !        !                            .
       !        !            +             .
       !        !          + +           .
       !        !        +   +         .
       !        !      +     +       .
       !        !    +       +     .
       !        !  +         +   .
       ! *KMAX* !+           + .
       !        !            .
       !        !          . .
       !        !        .   .
       !        !      .     .
       !        !    .       .
       !        !  .         .
       !        !________________________________________>
       !                  *MMAX*                      *M*
       !
       !
       ! Method:
       !
       ! *call* *phcs(ppnm,panm,kmp1,kkp1,pmu)*
       !
       ! *ppnm*   :*Legendre polynomials values.
       ! *panm*   :(mu**2-1)*dpnm/dmu.
       ! *kmp1*   :mmax+1.
       ! *kkp1*   :kmax+1.
       ! *pmu*    :value at which *ppnm* and *panm* are computed.
       !
       ! The *Legendre polynomials are defined as follows:
       !     * p(n,m)(mu)=sqrt((2n+1)*fact(n-m)/fact(n+m))    *
       !     *            /(fact(n)*2**(n+1))                 *
       !     *            *(1-mu**2)**(m/2)                   *
       !     *            *d**(n+m)(mu**2-1)**n/(dmu)**(n+m)  *
       !
       ! with *fact(k)=k*(k-1)**1
       !
       ! They are computed with the following numerically stable
       ! recurrence relation (Belousov,1962):
       !
       !     * p(n,m)(mu)=c(n,m)*p(n-2,m-2)                   *
       !     *           -d(n,m)*mu*p(n-1,m-2)                *
       !     *           +e(n,m)*p(n-1,m)                     *
       !
       ! with
       !     *c(n,m)=sqrt((2n+1)*(m+n-1)*(m+n-3)              *
       !     *           /(2n-3)/(m+n  )/(m+n-2))             *
       !
       !     *d(n,m)=sqrt((2n+1)*(m+n-1)*(n-m+1)              *
       !     *           /(2n-1)/(m+n  )/(m+n-2))             *
       !
       !     *e(n,m)=sqrt((2n+1)*(n-m)                        *
       !     *           /(2n-1)/(n+m))                       *
       !
       ! The derivatives *(panm)* are then computed as follows:
       !
       !     *pa(n,m)=n*f(n+1,m)*p(n+1,m)                     *
       !     *       -(n+1)*f(n,m)*p(n-1,m)                   *
       !
       ! with:
       !
       !     *f(n,m)=sqrt((n**2-m**2)/(4*n**2-1))             *
       !
       ! Results.:
       ! The *Legendre polynomials and their derivatives are stored
       ! column-wise. The following normalisation is used:
       ! Integral over [-1,+1] of *(ppnm**2)* =.5
       !
       ! References:
       ! Belousov,S.L.,1962:Tables of normalised associated Legendre
       ! polynomials.(Mathematical tables series,
       ! Vol 18,Pergamon Press, New York, USA) 379pp
       !

       INTEGER,  intent(in)  :: kmp1, kkp1
       REAL(WP), intent(in)  :: pmu
       REAL(WP), intent(out) :: ppnm(:), panm(:)

       REAL(WP) :: z2mm1, z2q2, zan, zateta, zcnm, zcos2, zcosfak, zcospar, &
            zcostet, zdnm, zenm, zmm1, zn, zn2, zn2m1, znm1, zp, zp2, zq, &
            zq2m1, zsinfak, zsinpar, zsintet, zsqp, zteta, zw, zwm2, zwm2q2, &
            zwnm1, zwq
       INTEGER  :: ik, inmax, inmaxm, ito, iton, jk, jm, jn

       REAL(WP) :: ztemp(3,kmp1+kkp1)

       INTRINSIC ACOS, COS, SIN, SQRT

       !  Executable statements

       ! 1. Initiate recurrence by computing:
       !      *p(0,0)* *p(1,1)* *pa(0,0)* *pa(1,1)*.

       inmax = kkp1+kmp1

       zcos2 = SQRT(1.0_wp-pmu**2)
       zteta = ACOS(pmu)
       zan = 1.0_wp

       ztemp = 0
       ztemp(1,1) = 0.5_wp
       DO jn = 2, inmax
          zsinpar = 0.0_wp
          zcospar = 0.0_wp
          zp = jn-1.0_wp
          zp2 = zp*zp
          zsqp = 1.0_wp/SQRT(zp2+zp)
          zan = zan*SQRT(1.0_wp-1.0_wp/(4*zp2))
          zcosfak = 1.0_wp
          zsinfak = zp*zsqp

          ik = jn
          DO jk = 1, ik, 2
             zq = jk-1.0_wp
             zateta = (zp-zq)*zteta
             zcostet = COS(zateta)
             zsintet = SIN(zateta)
             IF (jn == jk) zcostet = zcostet*0.5_wp
             IF (jk /= 1) THEN
                zcosfak = (zq-1.0_wp)/zq*(2*zp-zq+2.0_wp)/(2*zp-zq+1.0_wp)*zcosfak
                zsinfak = zcosfak*(zp-zq)*zsqp
             END IF
             zcospar = zcospar+zcostet*zcosfak
             zsinpar = zsinpar+zsintet*zsinfak
          END DO

          ztemp(1,jn) = zan*zcospar
          ztemp(2,jn-1) = zan*zsinpar
       END DO

       ppnm(1) = 0.5_wp
       ppnm(1+kkp1) = ztemp(2,1)
       panm(1) = 0.0_wp
       panm(1+kkp1) = pmu*ztemp(2,1)

       ! 2. Complete recurrence

       ! 2.1 First 2 columns.

       DO jn = 2, kkp1
          ppnm(jn) = ztemp(1,jn)
          ppnm(jn+kkp1) = ztemp(2,jn)
          zn2 = 2.0_wp*jn
          zn2m1 = zn2-1.0_wp
          zn = jn
          znm1 = jn-1.0_wp
          panm(jn) = znm1*(pmu*ztemp(1,jn)-SQRT(zn2m1/(zn2-3))*ztemp(1,jn-1))
          panm(jn+kkp1) = zn*pmu*ztemp(2,jn)-SQRT(((zn2+1) &
               *(zn**2-1._wp))/zn2m1)*ztemp(2,jn-1)
       END DO

       ! 2.2 Other columns.

       DO jm = 3, kmp1
          zmm1 = jm-1.0_wp
          z2mm1 = zmm1*2
          ztemp(3,1) = SQRT(1.0_wp+1.0_wp/z2mm1)*zcos2*ztemp(2,1)
          ito = (jm-1)*kkp1
          ppnm(ito+1) = ztemp(3,1)
          panm(ito+1) = zmm1*pmu*ztemp(3,1)
          inmaxm = inmax-jm

          DO jn = 2, inmaxm
             iton = ito+jn
             znm1 = jn-1.0_wp
             zq = z2mm1+znm1-1.0_wp
             zwm2 = zq+znm1
             zw = zwm2+2
             zwq = zw*zq
             zq2m1 = zq*zq-1.0_wp
             zwm2q2 = zwm2*zq2m1
             zwnm1 = zw*znm1
             z2q2 = zq2m1*2
             zcnm = SQRT((zwq*(zq-2.0_wp))/(zwm2q2-z2q2))
             zdnm = SQRT((zwq*(znm1+1.0_wp))/zwm2q2)
             zenm = SQRT(zwnm1/((zq+1.0_wp)*zwm2))
             ztemp(3,jn) = zcnm*ztemp(1,jn)-pmu*(zdnm*ztemp(1,jn+1)-zenm &
                  *ztemp(3,jn-1))
             ppnm(iton) = ztemp(3,jn)
             panm(iton) = (zmm1+znm1)*pmu*ztemp(3,jn) &
                  -SQRT(zwnm1*(zq+1._wp)/zwm2)*ztemp(3,jn-1)
          END DO

          DO jn = 1, inmax
             ztemp(1,jn) = ztemp(2,jn)
             ztemp(2,jn) = ztemp(3,jn)

          END DO

       END DO

     END SUBROUTINE phcs
!==============================================================================
     SUBROUTINE inileg (im, in, ik, igl)
       ! Description:
       !
       ! Set up polynomials needed for the Legendre transforms.
       !
       ! Arguments (intent(in)):
       !
       !   im  = max zonal wave number
       !   in  = max meridional wave number for m=0
       !   ik  = max meridional wave number
       !   ngl = number of gaussian latitudes
       !
       ! Method:
       !
       ! *pointer* dimensions made dynamic, and i/o unit numbers changed
       ! to comply with the i/o scheme.
       !
       ! This subroutine computes, normalises and writes suitably
       ! modified *Legendre polynomials needed for the *Legendre
       ! transforms.
       !
       ! Externals:
       !
       ! *phcs*      called to compute the *Legendre polynomials
       !             and their meridional derivatives.
       ! *reord*     reorders the *Legendre polynomials
       !             and their meridional derivatives.
       !
       ! Authors:
       !
       ! M. Jarraud, ECMWF, March 1982, original source
       ! J. K. Gibson, ECMWF, April 82, changed
       ! L. Kornblueh, MPI, May 1998, f90 rewrite
       ! U. Schulzweida, MPI, May 1998, f90 rewrite
       !
       ! for more details see file AUTHORS
       !

       INTEGER, INTENT(in) :: im, in, ik, igl

       REAL(WP) :: zmu, zrcst
       INTEGER :: jgl, jp

       REAL(WP), ALLOCATABLE :: zanm(:), zanmt(:), zpnm(:), zpnmt(:)

       !  Executable statements

       ! deallocate if repeatedly called with changed parameters
       if (nm >= 0 ) then
         if (nm==im .and. nn==in .and. nk==ik .and. ngl==igl) return
         DEALLOCATE (pnm)
         DEALLOCATE (anm)
         DEALLOCATE (gmu)
         DEALLOCATE (gw )
         DEALLOCATE (cst)
         DEALLOCATE (nnp)
         DEALLOCATE (nmp)
         deallocate (pnmd,anmd,rnmd,pnmi,anmi,pnmiuv,anmiuv)
         print *,'inileg: nm, nn, nk, ngl = ' ,im ,in ,ik ,igl
       endif

       ! 0. setup dimensions

       nm = im
       nn = in
       nk = ik
       ngl = igl

       nmp1   = nm+1
       nnp1   = nn+1
       nkp1   = nk+1

       nhgl   = ngl/2

       nsp = nnp1*nmp1-(nm+nn-nk)*(nm+nn-nk+1)/2

       ALLOCATE (pnm(nsp,nhgl)); pnm=0.0_wp
       ALLOCATE (anm(nsp,nhgl)); anm=0.0_wp

       ALLOCATE (gmu(ngl))
       ALLOCATE (gw (ngl))
       ALLOCATE (cst(nhgl))

       ALLOCATE (nnp(nmp1))
       ALLOCATE (nmp(nmp1))

       ALLOCATE (zanm(nsp), zanmt(nmp1*nnp1+1), zpnm(nsp), zpnmt(nmp1*nnp1+1))

       allocate (pnmd   (nsp))
       allocate (anmd   (nsp))
       allocate (rnmd   (nsp))
       allocate (pnmi   (nsp))
       allocate (anmi   (nsp))
       allocate (pnmiuv (nsp))
       allocate (anmiuv (nsp))

       DO jp = 1, nmp1
          nnp(jp) = MIN(nkp1-jp,nn)+1
       END DO

       nmp(1) = 0
       DO jp = 2, nmp1
          nmp(jp) = nmp(jp-1) + nnp(jp-1)
       END DO

       CALL gauaw(gmu,gw,ngl)
       cst = 1._wp - gmu(1:nhgl)**2

       ! 1. Set up some constants and initiate scan

       DO jgl = 1, nhgl
          zmu = gmu(jgl)
          zrcst = 1._wp/cst(jgl)

          !2. Compute, reorder and normalise Legendre
          !   polynomials and their meridional derivatives.

          ! 2.1 Compute *Legendre polynomials

          CALL phcs (zpnmt, zanmt, nmp1, nnp1, zmu)

          ! *phcs* computes polynomials for a rhomboidal truncation.

          ! 2.2 Reorder *Legendre polynomials

          CALL reord (zpnm,zpnmt,zanm,zanmt)

          ! 2.3 Normalise *Legendre polynomials

          DO jp = 1, nsp
             zpnm(jp) = zpnm(jp)*2.0_wp

             ! *anm* is divided by -(1._wp-mu**2) in order to get *d(pnm)/dmu.*

             zanm(jp) = -zanm(jp)*2.0_wp*zrcst
          END DO

          ! 3. Write Legendre polynomials to module fields

          pnm(:,jgl) = zpnm(:)
          anm(:,jgl) = zanm(:)

       END DO

       DEALLOCATE (zanm, zanmt, zpnm, zpnmt)

     END SUBROUTINE inileg
!==============================================================================
     SUBROUTINE reord (ppnm, ptp, phnm, pth)

       ! Description:
       !
       ! Reorders the *Legendre polynomials and their derivatives generated
       ! on a rhomboidal truncation by *phcs* to get them on the truncation
       ! used by the model.
       !
       ! Method:
       !
       ! *reord* is called from *inicom*. The polynomials are
       ! received *(ptp,pth)* and returned *(ppnm,phnm)* as subroutine
       ! arguments. The other parameters used are obtained from modules.
       !
       ! Results:
       ! *pnm* and *hnm* are returned in the natural order (i.e
       ! column wise ) on the truncation used by the model.
       !
       ! Authors:
       !
       ! M. Jarraud, ECMWF, March 1982, original source
       ! L. Kornblueh, MPI, May 1998, f90 rewrite
       ! U. Schulzweida, MPI, May 1998, f90 rewrite
       !
       ! for more details see file AUTHORS
       !

       REAL(WP), intent(out) :: ppnm(:), phnm(:)
       REAL(WP), intent(in)  :: ptp(:),  pth(:)

       INTEGER :: imp, inp, jm, jn


       !  Executable statements

       ! 1. Reorders polynomials for a pentagonal truncation

       DO jm = 1, nmp1
          inp = nnp(jm)
          imp = nmp(jm)
          DO jn = 1, inp
             ppnm(imp+jn) = ptp((jm-1)*nnp1+jn)
             phnm(imp+jn) = pth((jm-1)*nnp1+jn)
          END DO
       END DO

     END SUBROUTINE reord
!==============================================================================
!+ calculate modified Legendre polynomials for direct tran

SUBROUTINE legmod(kpass)

  ! Description:
  !
  ! Calculate modified Legendre polynomials for direct tran.
  ! (grid-point to spherical harmonics)
  !
  ! Method:
  !
  ! *legmod* computes the modified Legendre polynomials 'pnmd',
  ! 'anmd', 'rnmd' (one latitudinal index) for direc
  ! Legendre transforms using the normalised Legendre polynomials
  ! 'pnm', 'anm' stored in module 'mo_legendre'
  !
  ! *call* legmod(kpass)
  !
  ! *kpass*     The value of the main loop control variable
  !             (half latitudinal index)
  !             when *legmod* is called
  !
  ! Authors:
  !
  ! D. W. Dent, ECMWF, May 1984, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  !
  ! for more details see file AUTHORS
  !

  IMPLICIT NONE

  !  Scalar arguments
  INTEGER, intent(in) :: kpass

  !  Local scalars:
  REAL(WP) :: fnnp, z2w, zsa2
  INTEGER :: imp, innp, inp, jm, jn, jp

  !  Intrinsic functions
  INTRINSIC REAL


  !  Executable statements

!-- 1. Compute modified Legendre polynomials

!      *pnmd*=2*w*pnmd
!      *anmd*=2*w*anmd
!      *rnmd*=2*w*(-n(n+1)/a**2)*pnmd

  z2w = 2._wp*gw(kpass)

  pnmd = pnm(:,kpass)
  anmd = anm(:,kpass)

  DO jp = 1, nsp
    pnmd(jp) = pnmd(jp)*z2w
    anmd(jp) = anmd(jp)*z2w
  END DO

  zsa2 = 1._wp/a**2

  DO jm = 1, nmp1

    inp = nnp(jm)
    imp = nmp(jm)
    innp = jm - 1
    DO jn = 1, inp
      fnnp = REAL(innp+jn,wp)
      rnmd(imp+jn) = -zsa2*pnmd(imp+jn)*(fnnp-1._wp)*(fnnp)

    END DO
  END DO

  RETURN
END SUBROUTINE legmod
!==============================================================================
! calculate modified Legendre polynomials for inverse transform

SUBROUTINE leginv(kpass)

  ! Description:
  !
  ! Calculate modified Legendre polynomials for inverse tra.
  !
  ! Method:
  !
  ! *leginv* computes the modified Legendre polynomials for inver
  ! Legendre transforms using the normalised Legendre polynomials
  ! 'pnm', 'anm' stored in module 'mo_legendre'.
  !
  ! *call* *leginv(kpass)*
  !
  ! *kpass*      The value of the main loop control variable
  !              (half latitudinal index)
  !
  ! Authors:
  !
  ! D. W. Dent, ECMWF, May 1984, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  !
  ! for more details see file AUTHORS
  !

  IMPLICIT NONE

  !  Scalar arguments
  INTEGER, intent(in) :: kpass

  !  Local scalars:
  REAL(WP) :: fnnp, zcst, zsa
  INTEGER  :: imp, innp, inp, is, jm, jn, jp

  !  Intrinsic functions
  INTRINSIC REAL


  !  Executable statements

!-- 1. Compute modified Legendre polynomials
!      suitable for inverse Legendre transforms of u and v

!      *pnmiuv*=pnmi*a*m/(n(n+1)).
!      *anmiuv*=anmi*a*(1.-mu**2)/(n(n+1)).

  pnmi = pnm(:,kpass)
  anmi = anm(:,kpass)

  zcst = cst(kpass)

  DO jm = 1, nmp1

    inp = nnp(jm)
    imp = nmp(jm)
    innp = jm - 1

    IF (jm==1) THEN
      is = 2
    ELSE
      is = 1
    END IF

    DO jn = is, inp
      fnnp = REAL(innp+jn,wp)
      pnmiuv(imp+jn) = pnmi(imp+jn)*a*(jm-1)/((fnnp-1._wp)*(fnnp))
      anmiuv(imp+jn) = anmi(imp+jn)*a*zcst/((fnnp-1._wp)*(fnnp))

    END DO
  END DO

  pnmiuv(nmp(1)+1) = 0._wp
  anmiuv(nmp(1)+1) = 0._wp

!-- 2. Compute inverse Legendre polynomials *anmi*=*anmi/a*

  zsa = 1._wp/a

  DO jp = 1, nsp
    anmi(jp) = anmi(jp)*zsa
  END DO
  RETURN
END SUBROUTINE leginv
!==============================================================================
SUBROUTINE gp2sp (gp,sp)
!-----------------------------------------------------------
! transform scalar grid point fields to spherical harmonics.
! 'inifft' and 'inileg' must be called for initialisation.
! alternating row ordering (adjacent rows are corresponding
! northern, southern rows) assumed.
!-----------------------------------------------------------
  !-----------------
  ! module variables
  !-----------------
  use mo_fft, only: trig, nfax, fft991cy
  !----------
  ! arguments
  !----------
  real(wp), intent(in)  :: gp(:,:,:) ! grid point          (nlp2,nf,ngl)
  real(wp), intent(out) :: sp(:,:,:) ! spherical harmonics (nf,  2, nsp)
  !----------------
  ! local variables
  !----------------
  integer           :: nf    ! number of fields to transform
  integer           :: ilp2  ! number of longitudinal grid points + 2
  integer           :: ilon  ! number of longitudinal grid points
  integer           :: irow  ! row index
  integer           :: ihrow ! half index
  integer           :: jm    !
  integer           :: jnf   ! field index
  integer           :: is
  integer           :: jn
  real(wp), allocatable :: f_sym (:,:), f_asy(:,:), f_wrk(:,:,:), zwork(:)
  !---------------------
  ! allocate work arrays
  !---------------------
  ilp2 = size(gp,1)
  ilon = ilp2 - 2
  nf   = size(gp,2)
  allocate (f_sym(nf,2), f_asy(nf,2), f_wrk(ilp2,nf,2), zwork(ilp2*nf*2))
  !-----------------------
  ! gaussian latitude loop
  !-----------------------
  sp(:,:,:) = 0.0_wp
  DO irow = 1,ngl,2
    ihrow = (irow+1)/2   ! half index
    !------------------------------
    ! prepare Legendre coefficients
    !------------------------------
    CALL legmod(ihrow)   ! modify polynoms
    !---------------------------------------
    ! Fourier transformation irow and irow+1
    !---------------------------------------
    f_wrk = gp(:,:,irow:irow+1)
    CALL fft991cy(f_wrk,zwork,trig,nfax,1,ilp2,ilon,nf*2,-1)
    !-------------------------
    ! Fourier coefficient loop
    !-------------------------
    DO jm = 1,nmp1
      if (2*jm <= ilp2) then
        !-----------------------------------------------
        ! separation in symmetric and antisymmetric part
        !-----------------------------------------------
        DO jnf = 1,nf
          f_sym(jnf,1:2) = 0.5_wp * &
          &     (f_wrk(2*jm-1:2*jm,jnf,1) + f_wrk(2*jm-1:2*jm,jnf,2) )
          f_asy(jnf,1:2) = 0.5_wp * &
          &     ( f_wrk(2*jm-1:2*jm,jnf,1) - f_wrk(2*jm-1:2*jm,jnf,2) )
        END DO
        !------------------------
        ! Legendre transformation
        !------------------------
        DO jn = 1,nnp(jm),2
          !--------------
          ! symmetry loop
          !--------------
          is = nmp(jm)+jn
          sp(:,:,is) = sp(:,:,is) + f_sym(:,:)*pnmd(is)
        END DO
        DO jn = 2,nnp(jm),2
          !------------------
          ! antisymmetry loop
          !------------------
          is = nmp(jm)+jn
          sp(:,:,is) = sp(:,:,is) + f_asy(:,:)*pnmd(is)
        END DO
      endif
    END DO
  END DO
end SUBROUTINE gp2sp
!==============================================================================
SUBROUTINE sp2gp (sp,gp)
!-----------------------------------------------------------
! transform scalar spherical harmonics fields to grid point.
! 'inifft' and 'inileg' must be called for initialisation.
! alternating row ordering (adjacent rows are corresponding
! northern, southern rows) assumed.
!-----------------------------------------------------------
  !-----------------
  ! module variables
  !-----------------
  use mo_fft, only: trig, nfax, fft991cy
  !----------
  ! arguments
  !----------
  real(wp), intent(in)  :: sp(:,:,:) ! spherical harmonics (nf,  2, nsp)
  real(wp), intent(out) :: gp(:,:,:) ! grid point          (nlp2,nf,ngl)
  !----------------
  ! local variables
  !----------------
  integer           :: nf    ! number of fields to transform
  integer           :: ilp2  ! number of longitudinal grid points + 2
  integer           :: ilon  ! number of longitudinal grid points
  integer           :: irow  ! row index
  integer           :: ihrow ! half index
  integer           :: jm    !
  integer           :: jnf   ! field index
  integer           :: is
  integer           :: ins
  real(wp), allocatable :: f_sym (:,:), f_asy(:,:), f_wrk(:,:,:), zwork(:)
  !-------------------------------------------------
  ! single or double precision matrix multiplication
  !-------------------------------------------------
# ifdef CRAY
#   define X_GEMV sgemv
# else
#  define X_GEMV dgemv
# endif
  external X_GEMV
  !--------------------
  ! allocate work space
  !--------------------
  ilp2 = size(gp,1)
  ilon = ilp2 - 2
  nf   = size(gp,2)
  allocate (f_sym(nf,2), f_asy(nf,2), f_wrk(ilp2,nf,2), zwork(ilp2*nf*2))
  !-----------------------
  ! gaussian latitude loop
  !-----------------------
  gp(:,:,:) = 0.0_wp
  DO irow = 1,ngl,2
    ihrow = (irow+1)/2 ! half index
    !------------------------------
    ! prepare Legendre coefficients
    !------------------------------
    CALL leginv(ihrow)
    !-------------------------
    ! Fourier coefficient loop
    !-------------------------
    f_wrk = 0.0_wp
    DO jm = 1, nmp1
      if (2*jm <= ilp2) then
        !------------------------------------------------
        ! inverse Legendre transformation, symmetric part
        !------------------------------------------------
        is  = nmp(jm)+1     ! symmetric part
        ins = (nnp(jm)+1)/2
        IF (ins>0) THEN
          CALL X_GEMV('N',2*nf,ins,1.0_wp,sp(:,:,is:),4*nf,pnmi(is:),2,0.0_wp,f_sym,1)
        ELSE
          f_sym(:,:) = 0.0_wp
        END IF
        !-------------------
        ! antisymmetric part
        !-------------------
        is  = nmp(jm)+2     ! antisymmetric part
        ins = nnp(jm)/2
        IF (ins>0) THEN
          CALL X_GEMV('N',2*nf,ins,1.0_wp,sp(:,:,is:),4*nf,pnmi(is:),2,0.0_wp,f_asy,1)
        ELSE
          f_asy(:,:) = 0.0_wp
        END IF
        !-----------------------------
        ! compose Fourier coefficients
        !-----------------------------
        DO jnf = 1,nf
          !-------------------
          ! northern component
          !-------------------
          f_wrk(2*jm-1:2*jm,jnf,1) = f_sym(jnf,1:2) + f_asy(jnf,1:2)
          !-------------------
          ! southern component
          !-------------------
          f_wrk(2*jm-1:2*jm,jnf,2) = f_sym(jnf,1:2) - f_asy(jnf,1:2)
        ENDDO
      endif
    ENDDO
    !-------------------------------
    ! inverse Fourier transformation
    !-------------------------------
    CALL fft991cy(f_wrk,zwork,trig,nfax,1,ilp2,ilon,nf*2,1)
    gp(:,:,irow:irow+1) = f_wrk
  END DO
END SUBROUTINE sp2gp
!==============================================================================
END MODULE mo_legendre
