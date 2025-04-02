!NEC$ options "-finline-max-depth=4 -finline-max-function-size=50000"

! Source module for the radar forward operator EMVORADO
!
! ---------------------------------------------------------------
! Copyright (C) 2005-2024, DWD, NSSL
! Contact information: ulrich.blahak (at) dwd.de 
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE radar_dualpol_t_matrix2_mod

!=============================================================================
!
! This module includes the necessary code to calculate the complex scattering
! amplitudes of a particle (by default, an oblate spheroid) using a two-layer
! T-matrix method.
!
! History:
!   <unknown>
!   Last change of the original Fortran4 / Fortran77 code: 10 October 2016 by Jeffrey C. Snyder, NSSL
!     (original file names: radar_dualpol_t_matrix2_double_mod.f90, radar_dualpol_t_matrix2_quad_mod.f90
!   Technical Rewrite to Fortran90 December 2020 by Ulrich Blahak, DWD
!     (unified double precision / quad precision code)
!
! Default is double precision.
! Preprocessor flag -DTMATRIX2_QUADPREC to switch from double to quad precision.
!
!=============================================================================

  USE radar_utilities, ONLY : ind2sub3D

  IMPLICIT NONE

  SAVE

  PRIVATE

#ifdef TMATRIX2_QUADPREC
      INTEGER, PARAMETER :: PREC=SELECTED_REAL_KIND(32)    ! quad precision
#else
      INTEGER, PARAMETER :: PREC=SELECTED_REAL_KIND(15)    ! double precision
#endif
      INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(15)

      REAL(PREC) :: dcore,epsinnerreal,epsinnerimag,epsouterreal,epsouterimag
!       common /jeffcommon/ dcore,epsinnerreal,epsinnerimag,epsouterreal,epsouterimag
      REAL(PREC) :: dpart
!      common /vivek2/ dpart


      REAL(PREC) :: dcnr, dcni, ckprr, ckpri, ckr, dckr, conk,aovrb, sigma, ib
!        common /bdycom/ dcnr, dcni, ckprr, ckpri, ckr, dckr, conk, aovrb, sigma, ib

      INTEGER nm, kmv
      REAL(PREC) :: cmi(10), cmv, cm2, twm, prodm
!        common /cmvcom/ nm, kmv, cmi, cmv, cm2, twm, prodm

      REAL(PREC) :: qem, b89
!        common /variam/ qem, b89

      REAL(PREC) :: dcnr2, dcni2, ckr2, dckr2, conk2, aovrb2
!        common /bdycm2/ dcnr2, dcni2, ckr2, dckr2, conk2, aovrb2

      INTEGER nrank, nranki
      REAL(PREC) :: a(40,40,2), b(40,40,2), c(40,40,2), d(40,40,2)
      REAL(PREC) :: x(40,40,2), y(40,40,2), cmxnrm(40)
!        common /mtxcom/ nrank, nranki, a, b, c, d, x, y, cmxnrm
      REAL(PREC) :: tmat(40,40,2)

      REAL(PREC) :: anginc, acans(181,2,2), uang(181), rtsfct, dltang
      INTEGER nuang
!        common /uvccom/ anginc, acans, uang, rtsfct, dltang, nuang

      INTEGER irow, irow1, icol, icol1
      REAL(PREC) :: crow, crowm, crow1, ccol, ccolm, ccol1, crij,crssij
!        common /rowcol/ irow, irow1, crow, crowm, crow1, icol, icol1, ccol, ccolm, ccol1, crij, crssij

      REAL(PREC) :: height
      INTEGER ml, ijk
!        common /brin/ ml, ijk, height

      REAL(PREC) :: epps(4)
      INTEGER nsect
!        common /endpnt/ epps, nsect

      REAL(PREC) :: summ1, summ2
!        common /scatt/ summ1, summ2

      REAL(PREC) :: ptopdia, prhoice, pzhgt, ptempt, palamd, paovrb
!        common /pck_values/ ptopdia, prhoice, pzhgt, ptempt, palamd, paovrb

      INTEGER kgauss
!        common /cgauss/ kgauss

      REAL(PREC) :: theta, sinth, costh
!        common/thtcom/theta,sinth,costh

      REAL(PREC) :: pnmllg(41), bsslsp(41, 31, 3)
      REAL(PREC) :: cneumn(41, 31, 3)
      REAL(PREC) :: bslkpr(41, 31, 3), bslkpi(41, 31, 3)
      REAL(PREC) :: cneumr(41, 31, 3), cneumi(41, 31, 3)
!        common /fnccom/ pnmllg, bsslsp, cneumn, bslkpr, bslkpi,cneumr, cneumi

      REAL(PREC) :: bkpr2(41, 31, 3), bkpi2(41, 31, 3)
      REAL(PREC) :: cnpr2(41, 31, 3), cnpi2(41, 31, 3)
      REAL(PREC) :: bkprr2(41, 31, 3), bkpii2(41,31,3)
!        common /bringi/ bkpr2, bkpi2, cnpr2, cnpi2, bkprr2, bkpii2

      LOGICAL :: debug = .false.
!        common /debugc/ debug

      INTEGER :: number_of_angles = 2

    REAL(PREC), PARAMETER :: &
            cpi = 4.0_PREC * ATAN(1.0_PREC), &
            dtr = cpi / 180.0_PREC, &
            rtd = 180.0_PREC / cpi

  PUBLIC :: tmatrix2

CONTAINS

! JM200526
!  CALL tmatrix2(lambda, D, D-D_shell,ar(?),ar_core(?),m_shell,m_core,fa,fb,fa0,fb0)
!
  SUBROUTINE tmatrix2(wl,idpart,idcore,iaovrb,iaovrb2,mo,mi,f_a,f_b,f_a0,f_b0,ierr)
    
!       #########################################################################
!
!  FIXME: Document the input/output units! Apparently this is now different from JCS 2016 change.
!         And it seems the same between double and quad precision in the original two source files.
!         In the present new F90 source code, it is the same for double and quad by default.
!
!       JM200526:
!       Changed to assume identical units for wl and idpart/core and f*
!
!       JCS - 2016:
!       This subroutine requires the following as inputs:
!          wl - wavelength of the electromagnetic field in cm (single precision float)
!          idpart - equivolume diameter of the entire particle in mm (single prec. float)
!          idcore - equivolume diameter of the inner core in mm (single prec. float)
!          iaovrb - aspect ratio of the entire particle (defined as a/b, where a is
!                   the symmetry axis aligned in the vertical; single prec. float)
!          iaovrb2 - aspect ratio of the core (single prec. float)
!          epso - complex dielectric constant of the outer shell. Both parts of
!                 this should be positive because we assume exp(-jwt) time convention
!          epsi - complex dielectric constant of the inner core. Both parts of
!                 this should be positive because we assume exp(-jwt) time convention
!       The subroutine outputs the following (all in units of mm):
!          f_a - backwards complex scattering amplitude in "a" plane (vertical)
!          f_b - backwards complex scattering amplitude in "b" plane (horizontal)
!          f_a0 - forkwards complex scattering amplitude in "a" plane (vertical)
!          f_b0 - forwards complex scattering amplitude in "b" plane (horizontal)

!       JCS - 6/26/2014:
!       Updated to remove file IO and add passing of relevant parameters
!       via arguments. This improves efficiency when using the
!       subroutine in HUCM.

!       Date: 2/3/2003:
!       PCK mods started this date to implement the Maxwell-Garnett mixing formula
!       when calculating the complex dielectric function for an ice+water combination.


!       Date: 1/28/2003:
!       Mods started on this date to make the program simulate a single layer
!       (one refractive index) particle.
!       New refractive index subroutine added to simulate ice / air mixture.
!       These mods were done with advice of Larry Carey, Walt Petersen, and Bringi.
!       This source code version re-named to tmatrix1_pck.f as of 1/28/2003.

!       Original PCK documentation from NOvember 2002.
!       Parameters of primary interest (axis ratio, temperature, radar wavelength,
!       etc.) now read in from an external initializtion file.
!       NOTE: In subroutine rddata, the refractive index calls have been set
!       so that the inner portion of the scattering particle is ice, and
!       the outer layer is water.
!       #########################################################################

!       Electromagnetic scattering from a lossy dielectric imbedded within
!       another lossy dielectric.
!       exp(-jwt) time convention is used.
!       Perform the numerical integration and fill the a,b,c,d matrices
!       for the outer surface and x,y matrices for the inner surface.

!       Added implicit REAL(PREC) :: statements 11/98 wap

!       Unit 8 outputs just the size information & RCS (not used)
!       Unit 7 outputs final T-matrices
!       Unit 6 outputs intermediate parameters for checking convergence
!       after processing each mode for each size

    IMPLICIT NONE

    REAL(DP), INTENT(in) :: wl
    COMPLEX(DP), INTENT(in) ::  mi,mo
    REAL(DP), INTENT(in) :: idcore,idpart,iaovrb,iaovrb2
    COMPLEX(DP), INTENT(out) ::  f_a,f_a0,f_b,f_b0
    INTEGER, INTENT(out) :: ierr

    REAL(PREC) :: borrp,borip,borrpe,boripe,forrp,forip,forpe,foripe
    COMPLEX(PREC) f_a_q,f_a0_q,f_b_q,f_b0_q,epsi,epso
!
!       function declarations
!
!        REAL(PREC) :: cdmpyr, cdmpyi

!
!       program "global" "structs" if you will . . .
!

    REAL(PREC) :: temp(48, 8, 4)
    REAL(PREC) :: RESULT(48, 8)
    REAL(PREC) :: buff(5)

!
!       local vars
!

    REAL(PREC) :: bdyfct, ckrow, em, fctki, fprod, quanm, rhoice
    REAL(PREC) :: scale, temati, tematr, topdia, zhgt

    INTEGER :: i, ibfct, icom1, icom2, iefct, ifct, iflag, im, ith, j
    INTEGER :: js, k, kk, l, lfct, ll, nn, npts, nsect1


    CHARACTER(len=36) :: runid


!         print *,'kind(epso)=',kind(epso),' kind(f_a)=',kind(f_a),' kind(f_a_q)=',kind(f_a_q)

    ierr = 0

    palamd = REAL(wl,PREC)
    ml     = 1
    dcore  = REAL(idcore,PREC)
    epsi=mi**2.0_PREC
    epso=mo**2.0_PREC
    epsinnerreal = ABS(REAL(epsi))
    epsinnerimag = ABS(AIMAG(epsi))
    epsouterreal = ABS(REAL(epso))
    epsouterimag = ABS(AIMAG(epso))
    dpart        = REAL(idpart,PREC)
    aovrb        = REAL(iaovrb,PREC)
    aovrb2       = REAL(iaovrb2,PREC)


    IF (debug) THEN
      PRINT *,'t_matrix2_double.f: wl=',wl
      PRINT *,'dpart=',dpart,'dcore=',dcore
      PRINT *,'aovrb=', aovrb,'aovrb2=',aovrb2
      PRINT *,'epsinner=',epsinnerreal,epsinnerimag
      PRINT *,'epsouter=',epsouterreal,epsouterimag
    ENDIF

    !-----------------------------------------------
    ! inlined code from rddata():

#ifdef TMATRIX2_QUADPREC
    nrank = 17
    nm = 9
#else
    nrank = 7
    nm = 3
#endif
    nsect = 2
    ib = 8
    anginc = 90._PREC
    nranki = nrank + 1
    
    conk = cpi * dpart / (palamd * (aovrb**(1._PREC / 3._PREC)))
    conk2 = cpi * dcore / (palamd * (aovrb2**(1._PREC / 3._PREC)))

    sigma = 2.0_PREC * palamd / (cpi * 1.0_PREC)

    dcnr2 = epsinnerreal
    dcni2 = epsinnerimag
    dcnr = epsouterreal
    dcni = epsouterimag
        
    IF(nm == 1) THEN
      cmi(1) = 1.0_PREC
    ELSE
      DO i = 1, nm
        cmi(i) = REAL(i - 1,PREC)
      END DO
    ENDIF
    
    !       Preferable to start with kgauss=3 or more
    kgauss = 5
    CALL calenp

    ! end inlined code from rddata():
    !-----------------------------------------------
        
!       Set the number of scattering angles (nuang = 181 maximum and nuang-1
!       should be divisible into 180 by a whole number)
!         nuang= 2
        nuang = number_of_angles
        dltang = 180.0_PREC / (nuang - 1)
        uang(1) = anginc

!        print *,'H B'

        DO i = 2, nuang
          uang(i) = uang(i - 1) + dltang
        END DO

!        print *,'H C'

!       Clear the accumulating answer register (used in addprc).
! JCS - Just do it the simple way - clrtot has equivalence with acans
        acans = 0.0_PREC
        summ1 = 0.0_PREC
        summ2 = 0.0_PREC
        rtsfct = 8.0_PREC / conk

!       Set multiplier b89 dependent on ib value (symmetry indicator).
        IF (ib == 8) THEN
          b89 = 2.0_PREC
        ELSE
          b89 = 1.0_PREC
        END IF
          
        bdyfct = 1.0_PREC

!       Set up a loop for each m value.
        DO im = 1, nm

!          print *,'H E A'

!         set m dependent variables.
          cmv = cmi(im)
          kmv = cmv
          cm2 = cmv**2
          prodm = 1.0_PREC
          IF(kmv > 0) THEN
            em = 2.0_PREC
            quanm = cmv
            DO  ifct = 1, kmv
              quanm = quanm + 1.0_PREC
              prodm = quanm * prodm / 2.0_PREC
            END DO
          ELSE
            em = 1.0_PREC
          END IF

          qem = -2.0_PREC / em
          twm = 2.0_PREC * cmv
            
!         Initialize all matrix areas to zero
! JCS Just do it a simpler way
          a = 0.0_PREC
          b = 0.0_PREC
          c = 0.0_PREC
          d = 0.0_PREC
          x = 0.0_PREC
          y = 0.0_PREC

!         Set up a loop for each row of the matrices.
          crow = 0.0_PREC
          crowm = cmv
          DO irow = 1, nrank
            irow1 = irow + nrank
            crow = crow + 1.0_PREC
            crowm = crowm + 1.0_PREC
            crow1 = crow + 1.0_PREC

!           Set up a loop for each column of the matrices.
            ccol = 0.0_PREC
            ccolm = cmv
            DO icol = 1, nrank
              icol1 = icol + nrank
              ccol = ccol + 1.0_PREC
              ccolm = ccolm + 1.0_PREC
              ccol1 = ccol + 1.0_PREC
!
!             Calculate matrices a,b,c,d associated with the outer surface,
!             following notation by Peterson and Strom, q1(out,re) is stored
!             in a, q1(re,re) is stored in b,q1(re,out) is stored in c and
!             q1(out,out) is stored in d. This completes the matrices required
!             for the outer surface. Note: T-matrix is not calculated for the
!             outer surface. All matrices are transposed in the following code.
!
!             Calculate matrices x,y associated with the embedded or inner
!             surface. q2(out,re) is stored in x and q2(re,re) is stored in y.
!
!             Perform integration using a sequence of 1,3,7,15,31,63,127 and 255
!             point extended gauss-type quadrature formulae.
!             result(48,k) contains the values of the integrals .
!             There are 48 integrations to
!             be performed for each looping through irow and icol. These correspond
!             to 4 sub-matrix elements for each of the 6 matrices (a,b,c,d,x,y)
!             and assiciated real and imaginary parts.
!
              nsect1 = nsect - 1
              DO j=1, nsect1
                js = j
                ith = 0
                CALL quads(epps(j), epps(j + 1), k, RESULT, npts, ith, js, ierr)
                IF (ierr .NE. 0) RETURN
                
                temp(:, k, j) = RESULT(:, k)
              END DO
              
              DO j = 1, nsect1
                a(icol, irow1, 1) = temp(1, k, j) + a(icol, irow1, 1)
                a(icol, irow1, 2) = temp(2, k, j) + a(icol, irow1, 2)
                b(icol, irow1, 1) = temp(3, k, j) + b(icol, irow1, 1)
                b(icol, irow1, 2) = temp(4, k, j) + b(icol, irow1, 2)
                c(icol, irow1, 1) = temp(5, k, j) + c(icol, irow1, 1)
                c(icol, irow1, 2) = temp(6, k, j) + c(icol, irow1, 2)
                d(icol, irow1, 1) = temp(7, k, j) + d(icol, irow1, 1)
                d(icol, irow1, 2) = temp(8, k, j) + d(icol, irow1, 2)
                x(icol, irow1, 1) = temp(9, k, j) + x(icol, irow1, 1)
                x(icol, irow1, 2) = temp(10, k, j) + x(icol, irow1, 2)
                y(icol, irow1, 1) = temp(11, k, j) + y(icol, irow1, 1)
                y(icol, irow1, 2) = temp(12, k, j) + y(icol, irow1, 2)
                a(icol1, irow, 1) = temp(13, k, j) + a(icol1, irow, 1)
                a(icol1, irow, 2) = temp(14, k, j) + a(icol1, irow, 2)
                b(icol1, irow, 1) = temp(15, k, j) + b(icol1, irow, 1)
                b(icol1, irow, 2) = temp(16, k, j) + b(icol1, irow, 2)
                c(icol1, irow, 1) = temp(17, k, j) + c(icol1, irow, 1)
                c(icol1, irow, 2) = temp(18, k, j) + c(icol1, irow, 2)
                d(icol1, irow, 1) = temp(19, k, j) + d(icol1, irow, 1)
                d(icol1, irow, 2) = temp(20, k, j) + d(icol1, irow, 2)
                x(icol1, irow, 1) = temp(21, k, j) + x(icol1, irow, 1)
                x(icol1, irow, 2) = temp(22, k, j) + x(icol1, irow, 2)
                y(icol1, irow, 1) = temp(23, k, j) + y(icol1, irow, 1)
                y(icol1, irow, 2) = temp(24, k, j) + y(icol1, irow, 2)
                a(icol1, irow1, 1) = temp(25, k, j) + a(icol1, irow1, 1)
                a(icol1, irow1, 2) = temp(26, k, j) + a(icol1, irow1, 2)
                b(icol1, irow1, 1) = temp(27, k, j) + b(icol1, irow1, 1)
                b(icol1, irow1, 2) = temp(28, k, j) + b(icol1, irow1, 2)
                c(icol1, irow1, 1) = temp(29, k, j) + c(icol1, irow1, 1)
                c(icol1, irow1, 2) = temp(30, k, j) + c(icol1, irow1, 2)
                d(icol1, irow1, 1) = temp(31, k, j) + d(icol1, irow1, 1)
                d(icol1, irow1, 2) = temp(32, k, j) + d(icol1, irow1, 2)
                x(icol1, irow1, 1) = temp(33, k, j) + x(icol1, irow1, 1)
                x(icol1, irow1, 2) = temp(34, k, j) + x(icol1, irow1, 2)
                y(icol1, irow1, 1) = temp(35, k, j) + y(icol1, irow1, 1)
                y(icol1, irow1, 2) = temp(36, k, j) + y(icol1, irow1, 2)
                a(icol, irow, 1) = temp(37, k, j) + a(icol, irow, 1)
                a(icol, irow, 2) = temp(38, k, j) + a(icol, irow, 2)
                b(icol, irow, 1) = temp(39, k, j) + b(icol, irow, 1)
                b(icol, irow, 2) = temp(40, k, j) + b(icol, irow, 2)
                c(icol, irow, 1) = temp(41, k, j) + c(icol, irow, 1)
                c(icol, irow, 2) = temp(42, k, j) + c(icol, irow, 2)
                d(icol, irow, 1) = temp(43, k, j) + d(icol, irow, 1)
                d(icol, irow, 2) = temp(44, k, j) + d(icol, irow, 2)
                x(icol, irow, 1) = temp(45, k, j) + x(icol, irow, 1)
                x(icol, irow, 2) = temp(46, k, j) + x(icol, irow, 2)
                y(icol, irow, 1) = temp(47, k, j) + y(icol, irow, 1)
                y(icol, irow, 2) = temp(48, k, j) + y(icol, irow, 2)
              END DO

            END DO

!           Calculate the normalization factor (used in addprc)
            ckrow = irow
            IF (kmv <= 0) THEN
              fctki = 1.0_PREC
            ELSE

              IF (irow >= kmv) THEN
                ibfct = irow - kmv + 1
                iefct = irow + kmv
                fprod = ibfct
                fctki = 1.0_PREC
                
                DO lfct = ibfct, iefct
                  fctki = fctki * fprod
                  fprod = fprod + 1.0_PREC
                END DO
              END IF
              
            END IF

            IF (irow < kmv) THEN
              cmxnrm(irow) = 1.0_PREC
            ELSE
              cmxnrm(irow) = 4._PREC * ckrow * (ckrow + 1._PREC) * fctki / (em * (2._PREC * ckrow + 1._PREC))
            END IF
            
          END DO

          nn = 2 * nrank

!         Process computed matrices

!         Calculate t(2) and store in y. x is left free. also print y.
          CALL prcssm(x, y, nrank, nranki,ierr)
          IF (ierr .NE. 0) RETURN

!          print *,'H E G'

!         write(*,*) 'after prcssm'
!         call printm(y,nn,40)
!         normalise t - matrix  of the inner body before coupling
!         with the outer body matrices. this change made on nov 21/84.
          DO l = 1, nrank
            ll=l+nrank

            DO k = 1, nrank
              kk = k + nrank
              scale = cmxnrm(l) / cmxnrm(k)

              DO j = 1, 2
                y(l, k, j) = y(l, k, j) * scale
                y(l, kk, j) = y(l, kk, j) * scale
                y(ll, k, j) = y(ll, k, j) * scale
                y(ll, kk, j) = y(ll, kk, j) * scale
              END DO
            END DO
          END DO

!          print *,'H E H'

!         Calculate b-t(2)*c and store in x.
          DO l = 1, nn
            DO k = 1 , nn
              tematr = 0.0_PREC
              temati = 0.0_PREC

              DO  j = 1, nn
                tematr = cdmpyr(y(l, j, 1), y(l, j, 2), c(j, k, 1), c(j, k, 2)) + tematr
                temati = cdmpyi(y(l, j, 1), y(l, j, 2), c(j, k, 1), c(j, k, 2)) + temati
              END DO

              x(l, k, 1) = b(l, k, 1) - tematr
              x(l, k, 2) = b(l, k, 2) - temati

            END DO
          END DO


!         Calculate a-t(2)*d and store in b.
          DO l = 1, nn
            DO k = 1, nn
              tematr = 0.0_PREC
              temati = 0.0_PREC

              DO j = 1, nn
                tematr = cdmpyr(y(l, j, 1), y(l, j, 2), d(j, k, 1), d(j, k, 2)) + tematr
                temati = cdmpyi(y(l, j, 1), y(l, j, 2), d(j, k, 1), d(j, k, 2)) + temati
              END DO

              b(l, k, 1) = a(l, k, 1) - tematr
              b(l, k, 2) = a(l, k, 2) - temati

            END DO
          END DO

!         Calculate t(1,2)=b(inverse)*x and store in x.
          CALL prcssm(b, x, nrank, nranki, ierr)
          IF (ierr .NE. 0) RETURN

!         Print T matrix for this mode
!         write(7,11) cmi(im)  -smg08
! 11       format(2x,f10.5)
!          do 7  icom1 = 1, 2
!            do 7  icom2 = 1, nn
!          write(7,9) (x(icom2,icom3,icom1),icom3=1,nn)  -smg08
! 9        format(2x,8e15.7) -smg08
! 7        continue

!         call printm(x,nn,40)
!          call addprc(borrp,borip,borrpe,boripe,
!     1          forrp,forip,forpe,foripe)
          CALL addprc(f_a_q,f_a0_q,f_b_q,f_b0_q)
! JCS - Convert to units of mm
! FIXME
! JM190724:
! This is not doing anything here, is it (multipl by 1.0)?
! except for f_a.
! If wl and id* are in the same units (but they
! probably aren't. check that, too), f_* should be in the corresponding
! units, ie there should be no need for conversion. Unless somewhere
! else some unit conversion had been applied (I guess, it's actually a
! leftover from internal refractive index calcs. which would require wl
! to be in a specific unit).
! However, something must be applied somewhere since the output is off
! by 1e-2 (compared to the Mishchenko tmatrix for simple spheroids and
! to coated sphere Mie).
! For now, let's correct for that in the calling routine...
! JM190724:
! Is it correct, that fa is multipl by -1, while fb is not?
! JM190809:
! -> looks like yes, it is. cause test outputs for (quasi-simple) sphere
!    are the same for fa and fb.
! -> however, signs are off compared to Mishchenko tmatrix for f_a/b and
!    for Re(f_a/b0). Fixing this now here, at least for the moment.
! old:
          !f_a_q=conjg(f_a_q)
          !f_a=CMPLX(REAL(f_a_q,DP),AIMAG(f_a_q),DP)*(-1.0d0)
          !f_a0=CMPLX(REAL(f_a0_q,DP),AIMAG(f_a0_q),DP)*1.0d0
          !f_b_q=conjg(f_b_q)
          !f_b=CMPLX(REAL(f_b_q),AIMAG(f_b_q),DP)*1.0d0
          !f_b0=CMPLX(REAL(f_b0_q),AIMAG(f_b0_q),DP)*1.0d0
! new, sign-corrected:
          f_a_q=conjg(f_a_q)
          f_a=CMPLX(REAL(f_a_q,DP),AIMAG(f_a_q),DP)
          f_b_q=conjg(f_b_q)
          f_b=CMPLX(REAL(f_b_q),AIMAG(f_b_q),DP)*(-1.0d0)
          f_a0=CMPLX(-1.0d0*REAL(f_a0_q,DP),AIMAG(f_a0_q),DP)
          f_b0=CMPLX(-1.0d0*REAL(f_b0_q),AIMAG(f_b0_q),DP)

        if (debug) then
!          print *,'f_a_q,f_b_q,f_a0_q,f_b0_q'
!          print *,f_a_q,f_b_q,f_a0_q,f_b0_q
          print *,'In t_matrix2_double.f: f_a,f_b,f_a0,f_b0'
          print *,f_a,f_b,f_a0,f_b0
        endif
!          print *,'H E L'

      END DO

!        print *,'H F'
    
  END SUBROUTINE tmatrix2

!--------------------------------------------------------------
  !  SUBROUTINE rddata
  !
  !    REMOVED AND ACTIVE PARTS INLINED ABOVE!
  !
  !  END SUBROUTINE rddata
     
!-----------------------------------------------------------------
  subroutine gener(xxx, t, ith, js, ierr)

    !      IMPLICIT REAL(PREC) (A-H)
    !      IMPLICIT REAL(PREC) (O-Z)
    !      IMPLICIT INTEGER (I-N)
    IMPLICIT NONE

    ! ARGUMENTS
    REAL(PREC) :: xxx
    INTEGER :: ith, js, ierr

    ! LOCALS
    REAL(PREC) :: xxxr, sqreal, sqimag, qsreal, qsimag, sqrs2, sqis2, rat12r, rat12i
    REAL(PREC) :: sqrea2, sqima2, qsrea2, qsima2
    REAL(PREC) :: srmsin, tempp, bkr, bki, bkr2, bki2, tempnr, tempni, tempdr, tempdi
    REAL(PREC) :: hkr, hki, ckprr2, ckpri2, dkpr2, dkpi2, ckhr2, ckhi2
    REAL(PREC) :: br, hr, hi, hr2, hi2, br2, bi2
    REAL(PREC) :: cmcrco, pnr0c0, pnr0c1, pnr1c0, pnr1c1, a12, b1a, b1b
    REAL(PREC) :: hbkmlr, hbkmli, bbkmlr, bbkmli, hhkmlr, hhkmli
    REAL(PREC) :: bhkmlr, bhkmli, hepsr, hepsi, bepsr, bepsi, hhepsr, hhepsi, bbepsr, bbepsi
    REAL(PREC) :: hb2mlr, hb2mli, bb2mlr, bb2mli, ftheta
    REAL(PREC) :: heps2r,heps2i,beps2r,beps2i, b1, htbkr, htbki, btbkr, btbki
    REAL(PREC) :: sumar, sumr, sumai, sumi, sumbr, sumbi, bthkr, bthki, hthkr, hthki
    REAL(PREC) :: sumcr, sumci, sumdr, sumdi, htbkr2, htbki2, btbkr2, btbki2, tempr, tempi
    REAL(PREC) :: temrz, temiz, sumxr, sumxi, sumgr, sumgi, sumr2, sumi2, srx, six
    REAL(PREC) :: sumyr, sumyi, sumhr, sumhi, temr, temi
    REAL(PREC) :: dd, cr, ci, crh, cih, d2r, d2i, cr2, ci2, sumxr2, sumxi2

    INTEGER :: i, j, k, iswt
    INTEGER :: ierrv(nrank)
    LOGICAL :: lerr(nrank)

    REAL(PREC) :: rhankl(40, 2), rbskpr(40, 2), rbessl(40)
    REAL(PREC) :: hanklr(41), hankli(41), rhskpr(40, 2)
    REAL(PREC) :: rbess2(40, 2), hankr2(41), hanki2(41)
    REAL(PREC) :: rbskp2(40, 2), rhank2(40, 2)
    REAL(PREC) :: t(4, 6, 2)

    !      print *,'START gener'

    lerr = .FALSE.
    
    xxxr = xxx * rtd

    t(:, :, :) = 0.0

    sqreal = crootr(dcnr, dcni)
    sqimag = crooti(dcnr, dcni)
    qsreal = cddvdr(1._PREC,0._PREC, sqreal, sqimag,lerr(1))
    qsimag = cddvdi(1._PREC,0._PREC, sqreal, sqimag)

    !     Calculate sqrt(eps-2)
    !
    sqrs2 = crootr(dcnr2, dcni2)
    sqis2 = crooti(dcnr2, dcni2)

    !     Calculate eps-2/eps-1

    rat12r = cddvdr(dcnr2, dcni2, dcnr, dcni,lerr(1))
    rat12i = cddvdi(dcnr2, dcni2, dcnr, dcni)

    !     calculate sqrt(eps-2/eps-1)
    sqrea2 = crootr(rat12r, rat12i)
    sqima2 = crooti(rat12r, rat12i)

    !     calculate sqrt(eps-1/eps-2)
    qsrea2 = cddvdr(1._PREC, 0._PREC, sqrea2, sqima2,lerr(1))
    qsima2 = cddvdi(1._PREC, 0._PREC, sqrea2, sqima2)
    theta = xxx
    costh = cos(theta)
    sinth = sin(theta)
    srmsin = sinth
    tempp = cpi - theta

    if (lerr(1)) then
      ierr = 90001
      print '(A,I5,A)', 'ERR #',ierr,': div0 in calc of initial eps terms in gener - RETURNING.'
      return
    end if

    IF (ABS(tempp) < REAL(1.0d-8,PREC)) THEN
      sinth = 0.0_PREC
    END IF

    !     Generate the Legendre polynomials.
    CALL genlgp (kmv, twm, theta, sinth, costh, pnmllg)


    !     Evaluate kr and its derivative as a function of theta.
    CALL genkr

    !     Generate all necessary Bessel and Neumann function and their ratios.
    IF((irow==1).AND.(icol==1)) THEN
      CALL genbsl (ith,js,nrank,bsslsp,cneumn)
    END IF

    !     Calculate k*r for outer surface.
    ckprr = sqreal * ckr
    ckpri = sqimag * ckr
    iswt = 1

    if((irow.eq.1).and.(icol.eq.1)) then
      call genbkr(ckprr,ckpri,iswt,ith,js,ierr)
      if (ierr .ne. 0) return
    end if

    DO k = 1, nranki
      hanklr(k) = bslkpr(k, ith, js) - cneumi(k, ith, js)
      hankli(k) = bslkpi(k, ith, js) + cneumr(k, ith, js)
    END DO

    ierrv(:) = 0
    DO k = 1, nrank
      IF (bsslsp(k + 1, ith, js) .EQ. 0d0) THEN
        ierrv(k)=93500+k
      ELSE
        rbessl(k) = bsslsp(k, ith, js) / bsslsp(k + 1, ith, js)
        rhankl(k,1) = cddvdr(bsslsp(k,ith,js),cneumn(k,ith,js),bsslsp(k+1,ith,js),cneumn(k+1,ith,js),lerr(k))
        rhankl(k,2) = cddvdi(bsslsp(k,ith,js),cneumn(k,ith,js),bsslsp(k+1,ith,js),cneumn(k+1,ith,js))
        rbskpr(k,1) = cddvdr(bslkpr(k,ith,js),bslkpi(k,ith,js),bslkpr(k+1,ith,js),bslkpi(k+1,ith,js),lerr(k))
        rbskpr(k,2) = cddvdi(bslkpr(k,ith,js),bslkpi(k,ith,js),bslkpr(k+1,ith,js),bslkpi(k+1,ith,js))

        tempnr=bslkpr(k,ith,js)-cneumi(k,ith,js)
        tempni=bslkpi(k,ith,js)+cneumr(k,ith,js)
        tempdr=bslkpr(k+1,ith,js)-cneumi(k+1,ith,js)
        tempdi=bslkpi(k+1,ith,js)+cneumr(k+1,ith,js)
        rhskpr(k,1)=cddvdr(tempnr,tempni,tempdr,tempdi,lerr(k))
        rhskpr(k,2)=cddvdi(tempnr,tempni,tempdr,tempdi)

        IF (lerr(k)) THEN
          ierrv(k)=93500+k
        ELSE

          bkr = cdmpyr(sqreal,sqimag,rbskpr(k,1),rbskpr(k,2))
          bki = cdmpyi(sqreal,sqimag,rbskpr(k,1),rbskpr(k,2))
          rbskpr(k,1) = bkr
          rbskpr(k,2) = bki
          
          hkr=cdmpyr(sqreal,sqimag,rhskpr(k,1),rhskpr(k,2))
          hki=cdmpyi(sqreal,sqimag,rhskpr(k,1),rhskpr(k,2))
          rhskpr(k,1)=hkr
          rhskpr(k,2)=hki
        END IF
      END IF
    END DO
    
    IF (ANY(ierrv(:) /= 0)) THEN
      ierr=MAXVAL(ierrv)
      PRINT '(A,I5,A,I1,A)', &
           'ERR #',ierr,': div0 at rank=',ierr-93500,' in loop 350 in gener - RETURNING.'
      RETURN
    END IF
          
    call genkr2

    !      print *,'gener 1'

    !     calculate k*r2 for inner surface.
    ckprr2=ckr2*sqreal
    ckpri2=ckr2*sqimag
    !     calculate derivative of k*r2
    dkpr2=dckr2*sqreal
    dkpi2=dckr2*sqimag
    iswt=2
    if((irow.eq.1).and.(icol.eq.1)) then
      call genbkr(ckprr2,ckpri2,iswt,ith,js,ierr)
      if (ierr .ne. 0) return
    end if

    !      print *,'gener 2'

    !     calculate k*r2 for inner surface.
    ckhr2=ckr2*sqrs2
    ckhi2=ckr2*sqis2
    iswt=3
    if((irow.eq.1).and.(icol.eq.1)) then
      call genbkr(ckhr2,ckhi2,iswt,ith,js,ierr)
      if (ierr .ne. 0) return
    end if


    !      print *,'gener 3'

    DO  k=1,nranki
      hankr2(k)=bkpr2(k,ith,js)-cnpi2(k,ith,js)
      hanki2(k)=bkpi2(k,ith,js)+cnpr2(k,ith,js)
    END DO

    !      print *,'gener 4'

    ierrv(:) = 0
    DO k=1,nrank
      rbess2(k,1)=cddvdr(bkpr2(k,ith,js),bkpi2(k,ith,js),bkpr2(k+1,ith,js),bkpi2(k+1,ith,js),lerr(k))
      rbess2(k,2)=cddvdi(bkpr2(k,ith,js),bkpi2(k,ith,js),bkpr2(k+1,ith,js),bkpi2(k+1,ith,js))
      tempnr=bkpr2(k,ith,js)-cnpi2(k,ith,js)
      tempni=bkpi2(k,ith,js)+cnpr2(k,ith,js)
      tempdr=bkpr2(k+1,ith,js)-cnpi2(k+1,ith,js)
      tempdi=bkpi2(k+1,ith,js)+cnpr2(k+1,ith,js)
      rhank2(k,1)=cddvdr(tempnr,tempni,tempdr,tempdi,lerr(k))
      rhank2(k,2)=cddvdi(tempnr,tempni,tempdr,tempdi)

      if (lerr(k)) then
        ierrv(k)=93510+k
      end if

    END DO
    IF (ANY(ierrv(:) /= 0)) THEN
      ierr=MAXVAL(ierrv)
      PRINT '(A,I5,A,I1,A)', &
           'ERR #',ierr,': div0 at rank=',ierr-93510,' in loop 351 in gener - RETURNING.'
      RETURN
    END IF

    !      print *,'gener 5'

    ierrv(:) = 0
    DO k=1,nrank
      rbskp2(k,1)=cddvdr(bkprr2(k,ith,js),bkpii2(k,ith,js),bkprr2(k+1,ith,js),bkpii2(k+1,ith,js),lerr(k))
      IF (lerr(k)) THEN
        ierrv(k)=93520+k
      ELSE
        rbskp2(k,2)=cddvdi(bkprr2(k,ith,js),bkpii2(k,ith,js),bkprr2(k+1,ith,js),bkpii2(k+1,ith,js))

        bkr2=cdmpyr(sqrea2,sqima2,rbskp2(k,1),rbskp2(k,2))
        bki2=cdmpyi(sqrea2,sqima2,rbskp2(k,1),rbskp2(k,2))
        rbskp2(k,1)=bkr2
        rbskp2(k,2)=bki2
      END IF
    END DO
    IF (ANY(ierrv(:) /= 0)) THEN
      ierr=MAXVAL(ierrv)
      PRINT '(A,I5,A,I1,A)', &
           'ERR #',ierr,': div0 at rank=',ierr-93520,' in loop 352 in gener - RETURNING.'
      RETURN
    END IF

    !      print *,'gener 6'

    br = rbessl(irow)
    hr = rhankl(irow,1)
    hi = rhankl(irow,2)
    hr2=rhank2(irow,1)
    hi2=rhank2(irow,2)
    br2=rbess2(irow,1)
    bi2=rbess2(irow,2)

    !      print *,'gener 7'

    !     calculate frequently used variable combinations for use in a,b,c,
    !     d matrices.
    crij = crow+ccol
    crssij = crow*ccol
    cmcrco = cm2-qem*crssij*costh**2
    pnr0c0 = pnmllg(irow)*pnmllg(icol)
    pnr0c1 = pnmllg(irow)*pnmllg(icol+1)
    pnr1c0 = pnmllg(irow+1)*pnmllg(icol)
    pnr1c1 = pnmllg(irow+1)*pnmllg(icol+1)
    b1a = crow*costh*pnr1c1-crowm*pnr0c1
    b1b = ccol*costh*pnr1c1-ccolm*pnr1c0
    bkr = rbskpr(icol,1)
    bki = rbskpr(icol,2)
    hkr=rhskpr(icol,1)
    hki=rhskpr(icol,2)
    hbkmlr=cdmpyr(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js),bslkpr(icol+1,ith,js),bslkpi(icol+1,ith,js))
    hbkmli=cdmpyi(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js),bslkpr(icol+1,ith,js),bslkpi(icol+1,ith,js))
    bbkmlr = bsslsp(irow+1,ith,js)*bslkpr(icol+1,ith,js)
    bbkmli = bsslsp(irow+1,ith,js)*bslkpi(icol+1,ith,js)
    hhkmlr=cdmpyr(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js),hanklr(icol+1),hankli(icol+1))
    hhkmli=cdmpyi(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js),hanklr(icol+1),hankli(icol+1))
    bhkmlr=bsslsp(irow+1,ith,js)*hanklr(icol+1)
    bhkmli=bsslsp(irow+1,ith,js)*hankli(icol+1)
    hepsr = cdmpyr(qsreal,qsimag,hbkmlr,hbkmli)
    hepsi = cdmpyi(qsreal,qsimag,hbkmlr,hbkmli)
    bepsr = cdmpyr(qsreal,qsimag,bbkmlr,bbkmli)
    bepsi = cdmpyi(qsreal,qsimag,bbkmlr,bbkmli)
    hhepsr=cdmpyr(qsreal,qsimag,hhkmlr,hhkmli)
    hhepsi=cdmpyi(qsreal,qsimag,hhkmlr,hhkmli)
    bbepsr=cdmpyr(qsreal,qsimag,bhkmlr,bhkmli)
    bbepsi=cdmpyi(qsreal,qsimag,bhkmlr,bhkmli)

    !      print *,'gener 8'

    !     calculate frequently used variable combinations for use in
    !     x,y matrices.
    bkr2=rbskp2(icol,1)
    bki2=rbskp2(icol,2)
    hb2mlr=cdmpyr(hankr2(irow+1),hanki2(irow+1),bkprr2(icol+1,ith,js),bkpii2(icol+1,ith,js))
    hb2mli=cdmpyi(hankr2(irow+1),hanki2(irow+1),bkprr2(icol+1,ith,js),bkpii2(icol+1,ith,js))
    bb2mlr=cdmpyr(bkpr2(irow+1,ith,js),bkpi2(irow+1,ith,js),bkprr2(icol+1,ith,js),bkpii2(icol+1,ith,js))
    bb2mli=cdmpyi(bkpr2(irow+1,ith,js),bkpi2(irow+1,ith,js),bkprr2(icol+1,ith,js),bkpii2(icol+1,ith,js))
    !! 902   format(2x,4(2x,i3),4(4x,e10.3))
    ftheta=theta*rtd

    !       print *,'gener 9'

    !      write(6,902) js,ith,irow,icol,ftheta,bkpr2(irow+1,ith,js),
    !      bkpi2(irow+1,ith,js),bsslsp(irow+1,ith,js)
    heps2r=cdmpyr(qsrea2,qsima2,hb2mlr,hb2mli)
    heps2i=cdmpyi(qsrea2,qsima2,hb2mlr,hb2mli)
    beps2r=cdmpyr(qsrea2,qsima2,bb2mlr,bb2mli)
    beps2i=cdmpyi(qsrea2,qsima2,bb2mlr,bb2mli)


    !        print *,'gener 11'

    !     fill out elements for equivalent i-submatrix position.
    IF ( (ib==9 .OR. MOD(irow+icol,2) /= 0) .AND. (kmv /= 0) ) THEN

      b1 = b1a+b1b
      htbkr = cdmpyr(hr,hi,bkr,bki)
      htbki = cdmpyi(hr,hi,bkr,bki)
      btbkr = br*bkr
      btbki = br*bki
      tempp=(crow*crow1*bkr+ccol*ccol1*hr-crssij*(crij+2.0_PREC)/ckr)*dckr*sinth
      sumar=tempp*pnr1c1
      sumr=(ckr*(1.0_PREC  +htbkr)-ccol*hr-crow*bkr+crssij/ckr)*b1*ckr+sumar
      sumai = pnr1c1*(crow*crow1*bki+ccol*ccol1*hi)*dckr*sinth
      sumi = (ckr*htbki-ccol*hi-crow*bki)*b1*ckr+sumai
      t(1,1,1)=b89*cmv*srmsin*cdmpyr(sumr,sumi,hbkmlr,hbkmli)
      t(1,1,2)=b89*cmv*srmsin*cdmpyi(sumr,sumi,hbkmlr,hbkmli)
      sumbr = pnr1c1*(crow*crow1*bkr+ccol*ccol1*br-crssij*(crij+2.0_PREC)/ckr)*dckr*sinth
      sumr=(ckr*(1.0_PREC  +btbkr)-ccol*br-crow*bkr+crssij/ckr)*b1*ckr+sumbr
      sumbi = pnr1c1*crow*crow1*bki*dckr*sinth
      sumi = (ckr*btbki-crow*bki)*b1*ckr+sumbi
      t(1,2,1)=b89*cmv*srmsin*cdmpyr(sumr,sumi,bbkmlr,bbkmli)
      t(1,2,2)=b89*cmv*srmsin*cdmpyi(sumr,sumi,bbkmlr,bbkmli)
      bthkr=br*hkr
      bthki=br*hki
      hthkr=cdmpyr(hr,hi,hkr,hki)
      hthki=cdmpyi(hr,hi,hkr,hki)
      sumcr=pnr1c1*(crow*crow1*hkr+ccol*ccol1*br-crssij*(crij+2.0_PREC)/ckr)*dckr*sinth
      sumr=(ckr*(1.0_PREC  +bthkr)-ccol*br-crow*hkr+crssij/ckr)*b1*ckr+sumcr
      sumci=pnr1c1*(crow*crow1*hki)*dckr*sinth
      sumi=(ckr*bthki-crow*hki)*b1*ckr+sumci
      t(1,3,1)=b89*cmv*srmsin*cdmpyr(sumr,sumi,bhkmlr,bhkmli)
      t(1,3,2)=b89*cmv*srmsin*cdmpyi(sumr,sumi,bhkmlr,bhkmli)
      sumdr=pnr1c1*(crow*crow1*hkr+ccol*ccol1*hr-crssij*(crij+2.0_PREC)/ckr)*dckr*sinth
      sumr=(ckr*(1.0_PREC  +hthkr)-ccol*hr-crow*hkr+crssij/ckr)*b1*ckr+sumdr
      sumdi=pnr1c1*(crow*crow1*hki+ccol*ccol1*hi)*dckr*sinth
      sumi=(ckr*hthki-ccol*hi-crow*hki)*b1*ckr+sumdi
      t(1,4,1)=b89*cmv*srmsin*cdmpyr(sumr,sumi,hhkmlr,hhkmli)
      t(1,4,2)=b89*cmv*srmsin*cdmpyi(sumr,sumi,hhkmlr,hhkmli)

      !      print *,'gener 12'

      !     place q2(out,re)in x-matrix and q2(re,re) in y-matrix.
      htbkr2=cdmpyr(hr2,hi2,bkr2,bki2)
      htbki2=cdmpyi(hr2,hi2,bkr2,bki2)
      btbkr2=cdmpyr(br2,bi2,bkr2,bki2)
      btbki2=cdmpyi(br2,bi2,bkr2,bki2)
      tempr=cddvdr(1._PREC,0._PREC,ckprr2,ckpri2,lerr(1))
      IF (lerr(1)) THEN
        ierr=90002
        PRINT '(A,I5,A)', 'ERR #',ierr,': div0 in q2 placement in gener - RETURNING.'
        RETURN
      END IF
      tempi=cddvdi(1._PREC,0._PREC,ckprr2,ckpri2)
      temr=1.0_PREC  +htbkr2
      temrz=cdmpyr(temr,htbki2,ckprr2,ckpri2)
      temiz=cdmpyi(temr,htbki2,ckprr2,ckpri2)
      sumxr=pnr1c1*(crow*crow1*bkr2+ccol*ccol1*hr2-crssij*(crij+2.0_PREC)*tempr)*sinth
      sumxi=pnr1c1*(crow*crow1*bki2+ccol*ccol1*hi2-crssij*(crij+2.0_PREC)*tempi)*sinth
      sumgr=cdmpyr(sumxr,sumxi,dkpr2,dkpi2)
      sumgi=cdmpyi(sumxr,sumxi,dkpr2,dkpi2)
      sumr=(temrz-ccol*hr2-crow*bkr2+crssij*tempr)*b1
      sumi=(temiz-ccol*hi2-crow*bki2+crssij*tempi)*b1
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      srx=sumgr+sumr2
      six=sumgi+sumi2
      t(1,5,1)=b89*cmv*srmsin*cdmpyr(srx,six,hb2mlr,hb2mli)
      t(1,5,2)=b89*cmv*srmsin*cdmpyi(srx,six,hb2mlr,hb2mli)
      sumyr=pnr1c1*(crow*crow1*bkr2+ccol*ccol1*br2-crssij*(crij+2.0_PREC)*tempr)*sinth
      sumyi=pnr1c1*(crow*crow1*bki2+ccol*ccol1*bi2-crssij*(crij+2.0_PREC)*tempi)*sinth
      sumhr=cdmpyr(sumyr,sumyi,dkpr2,dkpi2)
      sumhi=cdmpyi(sumyr,sumyi,dkpr2,dkpi2)
      temr=1.0_PREC+btbkr2
      temrz=cdmpyr(temr,btbki2,ckprr2,ckpri2)
      temiz=cdmpyi(temr,btbki2,ckprr2,ckpri2)
      sumr=(temrz-ccol*br2-crow*bkr2+crssij*tempr)*b1
      sumi=(temiz-ccol*bi2-crow*bki2+crssij*tempi)*b1
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      srx=sumhr+sumr2
      six=sumhi+sumi2
      t(1,6,1)=b89*cmv*srmsin*cdmpyr(srx,six,bb2mlr,bb2mli)
      t(1,6,2)=b89*cmv*srmsin*cdmpyi(srx,six,bb2mlr,bb2mli)

      !      print *,'gener 13'

      !     fill out elements for equivalent l-submatrix position.
      sumr=(ckr*(dcnr+htbkr)-ccol*hr-crow*bkr+crssij/ckr)*b1*ckr+sumar
      sumi=(ckr*(dcni+htbki)-ccol*hi-crow*bki)*b1*ckr+sumai
      t(2,1,1)=-b89*cmv*srmsin*cdmpyr(sumr,sumi,hepsr,hepsi)
      t(2,1,2)=-b89*cmv*srmsin*cdmpyi(sumr,sumi,hepsr,hepsi)
      sumr = (ckr*(dcnr+btbkr)-ccol*br-crow*bkr+crssij/ckr)*b1*ckr+sumbr
      sumi = (ckr*(dcni+btbki)-crow*bki)*b1*ckr+sumbi
      t(2,2,1)=-b89*cmv*srmsin*cdmpyr(sumr,sumi,bepsr,bepsi)
      t(2,2,2)=-b89*cmv*srmsin*cdmpyi(sumr,sumi,bepsr,bepsi)
      sumr=(ckr*(dcnr+bthkr)-ccol*br-crow*hkr+crssij/ckr)*b1*ckr+sumcr
      sumi=(ckr*(dcni+bthki)-crow*hki)*b1*ckr+sumci
      t(2,3,1)=-b89*cmv*srmsin*cdmpyr(sumr,sumi,bbepsr,bbepsi)
      t(2,3,2)=-b89*cmv*srmsin*cdmpyi(sumr,sumi,bbepsr,bbepsi)
      sumr=(ckr*(dcnr+hthkr)-ccol*hr-crow*hkr+crssij/ckr)*b1*ckr+sumdr
      sumi=(ckr*(dcni+hthki)-ccol*hi-crow*hki)*b1*ckr+sumdi
      t(2,4,1)=-b89*cmv*srmsin*cdmpyr(sumr,sumi,hhepsr,hhepsi)
      t(2,4,2)=-b89*cmv*srmsin*cdmpyi(sumr,sumi,hhepsr,hhepsi)
      temr=rat12r+htbkr2
      temi=rat12i+htbki2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz-ccol*hr2-crow*bkr2+crssij*tempr)*b1
      sumi=(temiz-ccol*hi2-crow*bki2+crssij*tempi)*b1
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      srx=sumr2+sumgr
      six=sumi2+sumgi
      t(2,5,1)=-b89*cmv*srmsin*cdmpyr(srx,six,heps2r,heps2i)
      t(2,5,2)=-b89*cmv*srmsin*cdmpyi(srx,six,heps2r,heps2i)
      temr=rat12r+btbkr2
      temi=rat12i+btbki2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz-ccol*br2-crow*bkr2+crssij*tempr)*b1
      sumi=(temiz-ccol*bi2-crow*bki2+crssij*tempi)*b1
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      srx=sumr2+sumhr
      six=sumi2+sumhi
      t(2,6,1)=-b89*cmv*srmsin*cdmpyr(srx,six,beps2r,beps2i)
      t(2,6,2)=-b89*cmv*srmsin*cdmpyi(srx,six,beps2r,beps2i)

    END IF

    !      print *,'gener 14'

    IF (ib /= 8 .OR. MOD(irow+icol,2) == 0) THEN

      !     fill out elements for eqiivalent j-submatrix position.
      a12=cmcrco*pnr1c1+qem*(crow*ccolm*costh*pnr1c0+ccol*crowm*costh*pnr0c1-crowm*ccolm*pnr0c0)
      b1a = ccol*ccol1*b1a
      b1b = crow*crow1*b1b
      b1 = (b1a-b1b)*sinth
      dd=-qem*dckr
      cr = cdmpyr(dcnr,dcni,hr,hi)
      ci = cdmpyi(dcnr,dcni,hr,hi)
      sumr=(ckr*(bkr-cr)+dcnr*crow-ccol)*a12*ckr+(b1a-dcnr*b1b)*sinth*dd
      sumi=(ckr*(bki-ci)+dcni*crow)*a12*ckr-(dcni*b1b)*sinth*dd
      t(3,1,1)=b89*srmsin*cdmpyr(sumr,sumi,hepsr,hepsi)
      t(3,1,2)=b89*srmsin*cdmpyi(sumr,sumi,hepsr,hepsi)
      cr = br*dcnr
      ci = br*dcni
      sumr=(ckr*(bkr-cr)+dcnr*crow-ccol)*a12*ckr+(b1a-dcnr*b1b)*sinth*dd
      sumi=(ckr*(bki-ci)+dcni*crow)*a12*ckr-(dcni*b1b)*sinth*dd
      t(3,2,1)=b89*srmsin*cdmpyr(sumr,sumi,bepsr,bepsi)
      t(3,2,2)=b89*srmsin*cdmpyi(sumr,sumi,bepsr,bepsi)
      sumr=(ckr*(hkr-cr)+crow*dcnr-ccol)*a12*ckr+(b1a-dcnr*b1b)*sinth*dd
      sumi=(ckr*(hki-ci)+dcni*crow)*a12*ckr-(dcni*b1b)*sinth*dd
      t(3,3,1)=b89*srmsin*cdmpyr(sumr,sumi,bbepsr,bbepsi)
      t(3,3,2)=b89*srmsin*cdmpyi(sumr,sumi,bbepsr,bbepsi)
      crh=cdmpyr(dcnr,dcni,hr,hi)
      cih=cdmpyi(dcnr,dcni,hr,hi)
      sumr=(ckr*(hkr-crh)+crow*dcnr-ccol)*a12*ckr+(b1a-dcnr*b1b)*sinth*dd
      sumi=(ckr*(hki-cih)+crow*dcni)*a12*ckr-(dcni*b1b)*sinth*dd
      t(3,4,1)=b89*srmsin*cdmpyr(sumr,sumi,hhepsr,hhepsi)
      t(3,4,2)=b89*srmsin*cdmpyi(sumr,sumi,hhepsr,hhepsi)
      d2r=-qem*dkpr2
      d2i=-qem*dkpi2
      cr2=cdmpyr(rat12r,rat12i,hr2,hi2)
      ci2=cdmpyi(rat12r,rat12i,hr2,hi2)
      temr=bkr2-cr2
      temi=bki2-ci2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz+crow*rat12r-ccol)*a12
      sumi=(temiz+crow*rat12i)*a12
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      sumxr=(b1a-rat12r*b1b)*sinth
      sumxi=-rat12i*b1b*sinth
      sumxr2=cdmpyr(sumxr,sumxi,d2r,d2i)
      sumxi2=cdmpyi(sumxr,sumxi,d2r,d2i)
      srx=sumr2+sumxr2
      six=sumi2+sumxi2
      t(3,5,1)=b89*srmsin*cdmpyr(srx,six,heps2r,heps2i)
      t(3,5,2)=b89*srmsin*cdmpyi(srx,six,heps2r,heps2i)
      cr2=cdmpyr(rat12r,rat12i,br2,bi2)
      ci2=cdmpyi(rat12r,rat12i,br2,bi2)
      temr=bkr2-cr2
      temi=bki2-ci2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz+rat12r*crow-ccol)*a12
      sumi=(temiz+rat12i*crow)*a12
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      sumxr=(b1a-rat12r*b1b)*sinth
      sumxi=-rat12i*b1b*sinth
      sumxr2=cdmpyr(sumxr,sumxi,d2r,d2i)
      sumxi2=cdmpyi(sumxr,sumxi,d2r,d2i)
      srx=sumr2+sumxr2
      six=sumi2+sumxi2
      t(3,6,1)=b89*srmsin*cdmpyr(srx,six,beps2r,beps2i)
      t(3,6,2)=b89*srmsin*cdmpyi(srx,six,beps2r,beps2i)

      !      print *,'gener 15'

      !     fill out elements for equivalent k-submatrix.
      sumr = (ckr*(bkr-hr)+crow-ccol)*a12*ckr+b1*dd
      sumi = (ckr*(bki-hi))*a12*ckr
      t(4,1,1)=b89*srmsin*cdmpyr(sumr,sumi,hbkmlr,hbkmli)
      t(4,1,2)=b89*srmsin*cdmpyi(sumr,sumi,hbkmlr,hbkmli)
      sumr = (ckr*(bkr-br)+crow-ccol)*a12*ckr+b1*dd
      sumi = bki*a12*ckr**2
      t(4,2,1)=b89*srmsin*cdmpyr(sumr,sumi,bbkmlr,bbkmli)
      t(4,2,2)=b89*srmsin*cdmpyi(sumr,sumi,bbkmlr,bbkmli)
      sumr=(ckr*(hkr-br)+crow-ccol)*a12*ckr+b1*dd
      sumi=(ckr*hki)*a12*ckr
      t(4,3,1)=b89*srmsin*cdmpyr(sumr,sumi,bhkmlr,bhkmli)
      t(4,3,2)=b89*srmsin*cdmpyi(sumr,sumi,bhkmlr,bhkmli)
      sumr=(ckr*(hkr-hr)+crow-ccol)*a12*ckr+b1*dd
      sumi=(ckr*(hki-hi))*a12*ckr
      t(4,4,1)=b89*srmsin*cdmpyr(sumr,sumi,hhkmlr,hhkmli)
      t(4,4,2)=b89*srmsin*cdmpyi(sumr,sumi,hhkmlr,hhkmli)
      temr=bkr2-hr2
      temi=bki2-hi2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz+crow-ccol)*a12
      sumi=temiz*a12
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      sumxr=b1*d2r
      sumxi=b1*d2i
      srx=sumr2+sumxr
      six=sumi2+sumxi
      t(4,5,1)=b89*srmsin*cdmpyr(srx,six,hb2mlr,hb2mli)
      t(4,5,2)=b89*srmsin*cdmpyi(srx,six,hb2mlr,hb2mli)
      temr=bkr2-br2
      temi=bki2-bi2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz+crow-ccol)*a12
      sumi=temiz*a12
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      srx=sumr2+sumxr
      six=sumi2+sumxi
      t(4,6,1)=b89*srmsin*cdmpyr(srx,six,bb2mlr,bb2mli)
      t(4,6,2)=b89*srmsin*cdmpyi(srx,six,bb2mlr,bb2mli)

    END IF

    !      print *,'END gener'

  end subroutine gener
!-----------------------------------------------------------------
      subroutine quads(af,bf,k,result,npts,ith,js,ierr)
!     Numerical integration using Gauss-Legendre method.

      IMPLICIT NONE

      INTEGER   ,INTENT(inout)  :: npts, ith
      INTEGER   ,INTENT(in)     :: js
      REAL(PREC),INTENT(in)     :: af, bf
      INTEGER   ,INTENT(out)    :: k
      REAL(PREC),INTENT(out)    :: RESULT(48,8)
      INTEGER   ,INTENT(out)    :: ierr
      
      INTEGER    :: jj, kk, ii, l, i, iold, inew, m, j
      REAL(PREC) :: temp1
      REAL(PREC) :: xf
      REAL(PREC) :: sum, diff, temp2
      
      REAL(PREC)  :: p(381), funct(48,127), fzero(4,6,2), acum(48), test(48), t(4,6,2)
      REAL(PREC), PARAMETER, DIMENSION(54)  :: &
           d1 = [ &
                  7.74596669241483E-01_PREC, 5.55555555555557E-01_PREC, 8.88888888888889E-01_PREC, &
                  2.68488089868333E-01_PREC, 9.60491268708019E-01_PREC, 1.04656226026467E-01_PREC, &
                  4.34243749346802E-01_PREC, 4.01397414775962E-01_PREC, 4.50916538658474E-01_PREC, &
                  1.34415255243784E-01_PREC, 5.16032829970798E-02_PREC, 2.00628529376989E-01_PREC, &
                  9.93831963212756E-01_PREC, 1.70017196299402E-02_PREC, 8.88459232872258E-01_PREC, &
                  9.29271953151245E-02_PREC, 6.21102946737228E-01_PREC, 1.71511909136392E-01_PREC, &
                  2.23386686428967E-01_PREC, 2.19156858401588E-01_PREC, 2.25510499798206E-01_PREC, &
                  6.72077542959908E-02_PREC, 2.58075980961766E-02_PREC, 1.00314278611795E-01_PREC, &
                  8.43456573932111E-03_PREC, 4.64628932617579E-02_PREC, 8.57559200499902E-02_PREC, &
                  1.09578421055925E-01_PREC, 9.99098124967666E-01_PREC, 2.54478079156187E-03_PREC, &
                  9.81531149553739E-01_PREC, 1.64460498543878E-02_PREC, 9.29654857429739E-01_PREC, &
                  3.59571033071293E-02_PREC, 8.36725938168868E-01_PREC, 5.69795094941234E-02_PREC, &
                  7.02496206491528E-01_PREC, 7.68796204990037E-02_PREC, 5.31319743644374E-01_PREC, &
                  9.36271099812647E-02_PREC, 3.31135393257977E-01_PREC, 1.05669893580235E-01_PREC, &
                  1.12488943133187E-01_PREC, 1.11956873020953E-01_PREC, 1.12755256720769E-01_PREC, &
                  3.36038771482077E-02_PREC, 1.29038001003512E-02_PREC, 5.01571393058995E-02_PREC, &
                  4.21763044155885E-03_PREC, 2.32314466399103E-02_PREC, 4.28779600250078E-02_PREC, &
                  5.47892105279628E-02_PREC, 1.26515655623007E-03_PREC, 8.22300795723591E-03_PREC  ]
      REAL(PREC), PARAMETER, DIMENSION(54)  :: &
           d2 = [ &
                  1.79785515681282E-02_PREC, 2.84897547458336E-02_PREC, 3.84398102494556E-02_PREC, &
                  4.68135549906281E-02_PREC, 5.28349467901166E-02_PREC, 5.59784365104763E-02_PREC, &
                  9.99872888120358E-01_PREC, 3.63221481845531E-04_PREC, 9.97206259372224E-01_PREC, &
                  2.57904979468569E-03_PREC, 9.88684757547428E-01_PREC, 6.11550682211726E-03_PREC, &
                  9.72182874748583E-01_PREC, 1.04982469096213E-02_PREC, 9.46342858373402E-01_PREC, &
                  1.54067504665595E-02_PREC, 9.10371156957005E-01_PREC, 2.05942339159128E-02_PREC, &
                  8.63907938193691E-01_PREC, 2.58696793272147E-02_PREC, 8.06940531950218E-01_PREC, &
                  3.10735511116880E-02_PREC, 7.39756044352696E-01_PREC, 3.60644327807826E-02_PREC, &
                  6.62909660024781E-01_PREC, 4.07155101169443E-02_PREC, 5.77195710052045E-01_PREC, &
                  4.49145316536321E-02_PREC, 4.83618026945841E-01_PREC, 4.85643304066732E-02_PREC, &
                  3.83359324198731E-01_PREC, 5.15832539520484E-02_PREC, 2.77749822021825E-01_PREC, &
                  5.39054993352661E-02_PREC, 1.68235251552208E-01_PREC, 5.54814043565595E-02_PREC, &
                  5.63443130465928E-02_PREC, 5.62776998312542E-02_PREC, 5.63776283603847E-02_PREC, &
                  1.68019385741038E-02_PREC, 6.45190005017574E-03_PREC, 2.50785696529497E-02_PREC, &
                  2.10881524572663E-03_PREC, 1.16157233199551E-02_PREC, 2.14389800125039E-02_PREC, &
                  2.73946052639814E-02_PREC, 6.32607319362634E-04_PREC, 4.11150397865470E-03_PREC, &
                  8.98927578406411E-03_PREC, 1.42448773729168E-02_PREC, 1.92199051247278E-02_PREC, &
                  2.34067774953141E-02_PREC, 2.64174733950583E-02_PREC, 2.79892182552381E-02_PREC  ]
      REAL(PREC), PARAMETER, DIMENSION(54)  :: &
           d3 = [ &
                  1.80739564445388E-04_PREC, 1.28952408261042E-03_PREC, 3.05775341017553E-03_PREC, &
                  5.24912345480885E-03_PREC, 7.70337523327974E-03_PREC, 1.02971169579564E-02_PREC, &
                  1.29348396636074E-02_PREC, 1.55367755558440E-02_PREC, 1.80322163903913E-02_PREC, &
                  2.03577550584721E-02_PREC, 2.24572658268161E-02_PREC, 2.42821652033366E-02_PREC, &
                  2.57916269760242E-02_PREC, 2.69527496676331E-02_PREC, 2.77407021782797E-02_PREC, &
                  2.81388499156271E-02_PREC, 9.99982430354891E-01_PREC, 5.05360952078625E-05_PREC, &
                  9.99598799671912E-01_PREC, 3.77746646326985E-04_PREC, 9.98316635318407E-01_PREC, &
                  9.38369848542380E-04_PREC, 9.95724104698407E-01_PREC, 1.68114286542147E-03_PREC, &
                  9.91495721178104E-01_PREC, 2.56876494379402E-03_PREC, 9.85371499598521E-01_PREC, &
                  3.57289278351730E-03_PREC, 9.77141514639705E-01_PREC, 4.67105037211432E-03_PREC, &
                  9.66637851558417E-01_PREC, 5.84344987583563E-03_PREC, 9.53730006425761E-01_PREC, &
                  7.07248999543356E-03_PREC, 9.38320397779592E-01_PREC, 8.34283875396818E-03_PREC, &
                  9.20340025470011E-01_PREC, 9.64117772970252E-03_PREC, 8.99744899776941E-01_PREC, &
                  1.09557333878379E-02_PREC, 8.76513414484705E-01_PREC, 1.22758305600827E-02_PREC, &
                  8.50644494768350E-01_PREC, 1.35915710097655E-02_PREC, 8.22156254364980E-01_PREC, &
                  1.48936416648152E-02_PREC, 7.91084933799848E-01_PREC, 1.61732187295777E-02_PREC, &
                  7.57483966380512E-01_PREC, 1.74219301594641E-02_PREC, 7.21423085370098E-01_PREC, &
                  1.86318482561388E-02_PREC, 6.82987431091078E-01_PREC, 1.97954950480975E-02_PREC  ]
      REAL(PREC), PARAMETER, DIMENSION(54)  :: &
           d4 = [ &
                  6.42276642509760E-01_PREC, 2.09058514458120E-02_PREC, 5.99403930242243E-01_PREC, &
                  2.19563663053178E-02_PREC, 5.54495132631931E-01_PREC, 2.29409642293877E-02_PREC, &
                  5.07687757533716E-01_PREC, 2.38540521060385E-02_PREC, 4.59130011989833E-01_PREC, &
                  2.46905247444876E-02_PREC, 4.08979821229888E-01_PREC, 2.54457699654648E-02_PREC, &
                  3.57403837831532E-01_PREC, 2.61156733767061E-02_PREC, 3.04576441556714E-01_PREC, &
                  2.66966229274503E-02_PREC, 2.50678730303482E-01_PREC, 2.71855132296248E-02_PREC, &
                  1.95897502711100E-01_PREC, 2.75797495664819E-02_PREC, 1.40424233152560E-01_PREC, &
                  2.78772514766137E-02_PREC, 8.44540400837110E-02_PREC, 2.80764557938172E-02_PREC, &
                  2.81846489497457E-02_PREC, 2.81763190330167E-02_PREC, 2.81888141801924E-02_PREC, &
                  8.40096928705192E-03_PREC, 3.22595002508787E-03_PREC, 1.25392848264749E-02_PREC, &
                  1.05440762286332E-03_PREC, 5.80786165997757E-03_PREC, 1.07194900062519E-02_PREC, &
                  1.36973026319907E-02_PREC, 3.16303660822264E-04_PREC, 2.05575198932735E-03_PREC, &
                  4.49463789203206E-03_PREC, 7.12243868645840E-03_PREC, 9.60995256236391E-03_PREC, &
                  1.17033887476570E-02_PREC, 1.32087366975291E-02_PREC, 1.39946091276191E-02_PREC, &
                  9.03727346587510E-05_PREC, 6.44762041305726E-04_PREC, 1.52887670508776E-03_PREC, &
                  2.62456172740443E-03_PREC, 3.85168761663987E-03_PREC, 5.14855847897819E-03_PREC, &
                  6.46741983180368E-03_PREC, 7.76838777792199E-03_PREC, 9.01610819519566E-03_PREC, &
                  1.01788775292361E-02_PREC, 1.12286329134080E-02_PREC, 1.21410826016683E-02_PREC  ]
      REAL(PREC), PARAMETER, DIMENSION(54)  :: &
           d5 = [ &     
                  1.28958134880121E-02_PREC, 1.34763748338165E-02_PREC, 1.38703510891399E-02_PREC, &
                  1.40694249578135E-02_PREC, 2.51578703842806E-05_PREC, 1.88873264506505E-04_PREC, &
                  4.69184924247851E-04_PREC, 8.40571432710723E-04_PREC, 1.28438247189701E-03_PREC, &
                  1.78644639175865E-03_PREC, 2.33552518605716E-03_PREC, 2.92172493791781E-03_PREC, &
                  3.53624499771678E-03_PREC, 4.17141937698409E-03_PREC, 4.82058886485126E-03_PREC, &
                  5.47786669391895E-03_PREC, 6.13791528004137E-03_PREC, 6.79578550488277E-03_PREC, &
                  7.44682083240758E-03_PREC, 8.08660936478883E-03_PREC, 8.71096507973207E-03_PREC, &
                  9.31592412806942E-03_PREC, 9.89774752404876E-03_PREC, 1.04529257229060E-02_PREC, &
                  1.09781831526589E-02_PREC, 1.14704821146939E-02_PREC, 1.19270260530193E-02_PREC, &
                  1.23452623722438E-02_PREC, 1.27228849827324E-02_PREC, 1.30578366883530E-02_PREC, &
                  1.33483114637252E-02_PREC, 1.35927566148124E-02_PREC, 1.37898747832410E-02_PREC, &
                  1.39386257383068E-02_PREC, 1.40382278969086E-02_PREC, 1.40881595165083E-02_PREC, &
                  9.99997596379750E-01_PREC, 6.93793643241083E-06_PREC, 9.99943996207055E-01_PREC, &
                  5.32752936697805E-05_PREC, 9.99760490924434E-01_PREC, 1.35754910949228E-04_PREC, &
                  9.99380338025023E-01_PREC, 2.49212400482998E-04_PREC, 9.98745614468096E-01_PREC, &
                  3.89745284473282E-04_PREC, 9.97805354495956E-01_PREC, 5.54295314930373E-04_PREC, &
                  9.96514145914890E-01_PREC, 7.40282804244503E-04_PREC, 9.94831502800622E-01_PREC, &
                  9.45361516858527E-04_PREC, 9.92721344282788E-01_PREC, 1.16748411742996E-03_PREC  ]
      REAL(PREC), PARAMETER, DIMENSION(54)  :: &
           d6 = [ &
                  9.90151370400771E-01_PREC, 1.40490799565515E-03_PREC, 9.87092527954033E-01_PREC, &
                  1.65611272815445E-03_PREC, 9.83518657578632E-01_PREC, 1.91971297101387E-03_PREC, &
                  9.79406281670862E-01_PREC, 2.19440692536384E-03_PREC, 9.74734459752401E-01_PREC, &
                  2.47895822665757E-03_PREC, 9.69484659502459E-01_PREC, 2.77219576459345E-03_PREC, &
                  9.63640621569812E-01_PREC, 3.07301843470258E-03_PREC, 9.57188216109859E-01_PREC, &
                  3.38039799108691E-03_PREC, 9.50115297521293E-01_PREC, 3.69337791702565E-03_PREC, &
                  9.42411565191083E-01_PREC, 4.01106872407503E-03_PREC, 9.34068436157727E-01_PREC, &
                  4.33264096809299E-03_PREC, 9.25078932907077E-01_PREC, 4.65731729975685E-03_PREC, &
                  9.15437587155765E-01_PREC, 4.98436456476553E-03_PREC, 9.05140358813263E-01_PREC, &
                  5.31308660518706E-03_PREC, 8.94184568335557E-01_PREC, 5.64281810138445E-03_PREC, &
                  8.82568840247341E-01_PREC, 5.97291956550816E-03_PREC, 8.70293055548114E-01_PREC, &
                  6.30277344908575E-03_PREC, 8.57358310886234E-01_PREC, 6.63178124290190E-03_PREC, &
                  8.43766882672707E-01_PREC, 6.95936140939044E-03_PREC, 8.29522194637402E-01_PREC, &
                  7.28494798055382E-03_PREC, 8.14628787655138E-01_PREC, 7.60798966571904E-03_PREC, &
                  7.99092290960843E-01_PREC, 7.92794933429486E-03_PREC, 7.82919394118284E-01_PREC, &
                  8.24430376303287E-03_PREC, 7.66117819303759E-01_PREC, 8.55654356130769E-03_PREC, &
                  7.48696293616938E-01_PREC, 8.86417320948252E-03_PREC, 7.30664521242183E-01_PREC, &
                  9.16671116356077E-03_PREC, 7.12033155362253E-01_PREC, 9.46368999383007E-03_PREC  ]
      REAL(PREC), PARAMETER, DIMENSION(54)  :: &
           d7 = [ &
                  6.92813769779114E-01_PREC, 9.75465653631741E-03_PREC, 6.73018830230419E-01_PREC, &
                  1.00391720440569E-02_PREC, 6.52661665410019E-01_PREC, 1.03168123309476E-02_PREC, &
                  6.31756437711193E-01_PREC, 1.05871679048852E-02_PREC, 6.10318113715188E-01_PREC, &
                  1.08498440893373E-02_PREC, 5.88362434447664E-01_PREC, 1.11044611340069E-02_PREC, &
                  5.65905885423653E-01_PREC, 1.13506543159806E-02_PREC, 5.42965666498311E-01_PREC, &
                  1.15880740330440E-02_PREC, 5.19559661537457E-01_PREC, 1.18163858908302E-02_PREC, &
                  4.95706407918762E-01_PREC, 1.20352707852796E-02_PREC, 4.71425065871658E-01_PREC, &
                  1.22444249816120E-02_PREC, 4.46735387662029E-01_PREC, 1.24435601907140E-02_PREC, &
                  4.21657686626164E-01_PREC, 1.26324036435421E-02_PREC, 3.96212806057616E-01_PREC, &
                  1.28106981638774E-02_PREC, 3.70422087950079E-01_PREC, 1.29782022395374E-02_PREC, &
                  3.44307341599437E-01_PREC, 1.31346900919602E-02_PREC, 3.17890812068477E-01_PREC, &
                  1.32799517439305E-02_PREC, 2.91195148518247E-01_PREC, 1.34137930851101E-02_PREC, &
                  2.64243372410927E-01_PREC, 1.35360359349562E-02_PREC, 2.37058845589829E-01_PREC, &
                  1.36465181025713E-02_PREC, 2.09665238243181E-01_PREC, 1.37450934430019E-02_PREC, &
                  1.82086496759252E-01_PREC, 1.38316319095064E-02_PREC, 1.54346811481378E-01_PREC, &
                  1.39060196013255E-02_PREC, 1.26470584372302E-01_PREC, 1.39681588065169E-02_PREC, &
                  9.84823965981194E-02_PREC, 1.40179680394566E-02_PREC, 7.04069760428552E-02_PREC, &
                  1.40553820726499E-02_PREC, 4.22691647653637E-02_PREC, 1.40803519625536E-02_PREC  ]
      REAL(PREC), PARAMETER, DIMENSION(3)  :: &
           d8 = [ 1.40938864107825E-02_PREC, 1.40928450691604E-02_PREC, 1.40944070900962E-02_PREC  ]


      p(1:54)    = d1
      p(55:108)  = d2
      p(109:162) = d3
      p(163:216) = d4
      p(217:270) = d5
      p(271:324) = d6
      p(325:378) = d7
      p(379:381) = d8
     
!        print *,'START quad'

      IF(af /= bf) THEN
        
        sum=(bf+af)/2.0_PREC
        diff=(bf-af)/2.0_PREC
!     one point formula.
!     set up variable combinations for use in evaluation of integrands
        ith=ith+1
        CALL gener(sum,t,ith,js,ierr)
        IF (ierr /= 0) RETURN

!!$ UB: replaced the next nested 3 loops by vectorized syntax:
!!$        DO ii=1,4
!!$          DO jj=1,6
!!$            DO kk=1,2
!!$              !     jj=1,6 corresponds to matrices a thro d and x thro y. ii=1,4
!!$              !     corresponds to filling out eq. i,l,j,k positoons or sub-matrices
!!$              !     kk=1,2 corresponds to real and imaginary parts.
!!$              fzero(ii,jj,kk)=t(ii,jj,kk)
!!$              l=(ii-1)*12+(jj-1)*2+kk
!!$              RESULT(l,1)=2.0_PREC * t(ii,jj,kk) * diff
!!$            END DO
!!$          END DO
!!$        END DO
        fzero(:,:,:)=t(:,:,:)
        DO l=1,48
          !! recover ii,jj,kk from 1D-address l:
          CALL ind2sub3d(l, 2, 6, kk, jj, ii)
          RESULT(l,1)=2.0_PREC * t(ii,jj,kk) * diff
        END DO

        i=0
        iold=0
        inew=1
        acum(:)=0.0_PREC

        DO k=2, kgauss
        
          IF (k > 2) THEN
            
            acum(:)=0.0_PREC
!     contribution from function values already computed.
            DO j=1,iold
              i=i+1
              DO l=1,48
                acum(l) = acum(l) + p(i)*funct(l,j)
              END DO
            END DO
!     contribution from new values.

          END IF
      
          iold = iold + inew
          DO j = inew, iold
            
            i=i+1
            xf=p(i)*diff
            
            temp1 = sum+xf
            ith=ith+1
            CALL gener(temp1,t,ith,js,ierr)
            IF (ierr .NE. 0) RETURN

            DO l=1,48
              CALL ind2sub3d(l, 2, 6, kk, jj, ii)
              funct(l,j)=t(ii,jj,kk)
            END DO
            
            temp2=sum-xf
            ith=ith+1
            CALL gener(temp2,t,ith,js,ierr)
            IF (ierr .NE. 0) RETURN

            DO l=1,48
              CALL ind2sub3d(l, 2, 6, kk, jj, ii)
              funct(l,j)=funct(l,j)+t(ii,jj,kk)
            END DO
            i=i+1
            DO l=1,48
              acum(l)=acum(l)+p(i)*funct(l,j)
            END DO

          END DO
          
          inew=iold+1
          i=i+1
          DO l=1,48
            CALL ind2sub3d(l, 2, 6, kk, jj, ii)
            RESULT(l,k)=(acum(l)+p(i)*fzero(ii,jj,kk))*diff
          END DO

        END DO

!     normal termination.
        k = kgauss
        npts = inew + iold
      
    ELSE
      
      k=2
      RESULT(:,1:2) = 0.0_PREC
      npts=0

    END IF
!      print *,'END quads'


  END SUBROUTINE quads
!-----------------------------------------------------------------
  SUBROUTINE genlgp (kmv, twm, theta, sinth, costh, pnmllg)

!     A routine to generate Legendre polynomials.
!     The index on the function is incremented by one.

    IMPLICIT NONE

    INTEGER,    INTENT(in)  :: kmv
    REAL(PREC), INTENT(in)  :: twm, theta, sinth, costh
    REAL(PREC), INTENT(out) :: pnmllg(41)

    REAL(PREC) :: dtwm, pla, plb, cnmul, cnm, cnmm, plc
    INTEGER :: ilg, ibeg, ilgr

    dtwm = twm + 1.0_PREC

    !     this is special case when theta equals cpiand m=0.
    !     when theta equals cpi all integrands are 0 and any values can be
    !     put in pnmllg(41).here we have put them equal to 0.
    IF ((sinth == 0.0_PREC  ).AND.(kmv == 0)) THEN
        
      DO ilg = 1,nranki
        pnmllg(ilg)=0.0_PREC
      END DO

    ELSE IF (theta == 0.0_PREC .AND. kmv /= 1) THEN

      DO ilg = 1,nranki
        pnmllg(ilg)=0.0_PREC
      END DO

    ELSE
        
      !     at this point theta lies strictly between 0 and cpi.
      IF (theta == 0.0_PREC) THEN

        pnmllg(1)=0.0_PREC
        pnmllg(2)=1.0_PREC
        pla=1.0_PREC
          
        plb = dtwm*costh*pla
        pnmllg(kmv+2) = plb
        ibeg = kmv+3
          
      ELSE
          
        IF (kmv <= 0)  THEN
          !     the special case when m = 0.
          pla=1.0_PREC/sinth
          plb = costh*pla
          pnmllg(1) = pla
          pnmllg(2) = plb
          ibeg = 3
        ELSE
          !     general case for m not equal to 0.
          DO ilg = 1,kmv
            pnmllg(ilg)=0.0_PREC
          END DO
          IF ((sinth == 0.0_PREC  ).AND.(kmv == 1)) THEN
            pla=0.0_PREC
          ELSE
            pla = prodm*sinth**(kmv-1)
          END IF
          pnmllg(kmv+1) = pla
          
          plb = dtwm*costh*pla
          pnmllg(kmv+2) = plb
          ibeg = kmv+3
        END IF
          
        !     do recursion formula for all remaining legendre polynomials.
        cnmul = ibeg+ibeg-3.0_PREC
        cnm=2.0_PREC
        cnmm = dtwm
        DO ilgr = ibeg,nranki
          plc = (cnmul*costh*plb-cnmm*pla)/cnm
          pnmllg(ilgr) = plc
          pla   = plb
          plb   = plc
          cnmul = cnmul + 2.0_PREC
          cnm   = cnm   + 1.0_PREC
          cnmm  = cnmm  + 1.0_PREC
        END DO

      END IF
    END IF
    
  END SUBROUTINE genlgp
!-----------------------------------------------------------------
  REAL(PREC) FUNCTION cdabx(a, b)

    REAL(PREC) :: a, b, e, f, g

    IF (a /= 0.0_PREC) THEN ! 4, 22, 4

      IF (b /= 0.0_PREC) THEN
        e = MAX(a, b)
        f = MIN(a, b)
        g = f / e
        cdabx = ABS(e) * SQRT(1.0_PREC + g * g)
      ELSE
        cdabx = ABS(a)
      END IF
      
    ELSE

      IF(b /= 0.0_PREC) THEN ! 28, 26, 28
        cdabx = ABS(b)
      ELSE
        cdabx = 0.0_PREC
      END IF

    END IF
  
  END FUNCTION cdabx
!-----------------------------------------------------------------
  REAL(PREC) FUNCTION cdmpyr(a, b, c, d)

    REAL(PREC) :: a, b, c, d

    cdmpyr = a * c - b * d
    
  END FUNCTION cdmpyr
!-----------------------------------------------------------------
  REAL(PREC) FUNCTION cdmpyi(a, b, c, d)

    REAL(PREC) :: a, b, c, d

    cdmpyi = b * c + a * d
    
  END FUNCTION cdmpyi
!-----------------------------------------------------------------
  REAL(PREC) FUNCTION cddvdr(a, b, c, d, lerr)

    REAL(PREC) :: a, b, c, d, e, f
    LOGICAL    :: lerr

    e = c * c + d * d
    IF (e == 0.0_PREC) THEN
      lerr = .TRUE.
      cddvdr = 0.0_PREC
    ELSE    
      f = a * c + b * d
      cddvdr = f / e
    END IF

  END FUNCTION cddvdr
!-----------------------------------------------------------------
  REAL(PREC) FUNCTION cddvdi(a, b, c, d)

    REAL(PREC) :: a, b, c, d, e, f

    e = c * c + d * d

    ! JM200929:
    ! we don't need error catchment here, since cddvdi is always
    ! preceeded by cddvdr called with exactly the same c and d's
    IF (e == 0.0_PREC) THEN
      cddvdi = 0.0_PREC
    ELSE
      f = b * c - a * d
      cddvdi = f / e
    END IF
    
  END FUNCTION cddvdi
!-----------------------------------------------------------------
  SUBROUTINE genbsl (ith,js,nrank,bsslsp,cneumn)

    !     Generate Bessel and Neumann functions for real arguments.
    !     The index on the function is incremented by one.

    IMPLICIT NONE

    INTEGER,    INTENT(in)  :: ith, js, nrank
    REAL(PREC), INTENT(out) :: bsslsp(41, 31, 3), cneumn(41, 31, 3)
    
    INTEGER    :: nval, i, ip, ierror1, ierror2, iend, ie, je, jb
    REAL(PREC) :: pckr, answr, ansa, ansb, conn, ansc, cskrx, snkrx
    REAL(PREC) :: ckr2, cmuln, snsa, snsb, snsc, quanbt, quannt, thtprt, pckrr

!     set up a loop to get 2 successive bessel functions
    nval = nrank-1
    pckr = ckr
    DO i=1,4
      ! try the first one:
      CALL bessel(nval,pckr,answr,ierror1)
      IF(ierror1 <= 0) THEN
        ansa = answr
      ELSE
        ansa = -HUGE(1.0_PREC)
      END IF
      ! try the second one:
      CALL bessel(nval+1,pckr,answr,ierror2)
      IF(ierror2 <= 0) THEN
        ansb = answr
      ELSE
        ansb = -HUGE(1.0_PREC)       
      END IF
      ! If one of them failed, try again to find a pair with higher order:
      IF (ierror1 > 0 .OR. ierror2 > 0) THEN
        nval = nval + nrank
      ELSE
        EXIT
      END IF
    END DO
    !     program unable to generate bessel functions
    IF (ierror1 > 0 .OR. ierror2 > 0) THEN
      WRITE(6,'(///,5x,"unable to generate bessel functions",//)')

    ELSE
      
      ! set up for proper recursion of the bessel functons
      IF(nval-nrank > 0) THEN
        iend=nval-nrank
        conn=2*(nval-1)+1.0_PREC
        DO ip=1,iend
          ansc = conn*ansa/pckr-ansb
          conn = conn-2.0_PREC
          ansb = ansa
          ansa = ansc
        END DO
      END IF
      
!     program is ready to recurse downward into bessel function
      bsslsp(nranki,ith,js)=ansb
      bsslsp(nranki-1,ith,js)=ansa
      conn=(REAL(nrank+nrank-1,PREC))
      ie=nranki-2
      je=ie
      DO jb=1,je
        ansc=conn*ansa/pckr-ansb
        bsslsp(ie,ith,js)=ansc
        ansb=ansa
        ansa=ansc
        ie=ie-1
        conn=conn-2.0_PREC
      END DO
!     generate neumann functions
      cskrx= cos(pckr)/pckr
      snkrx= sin(pckr)/pckr
      ckr2=pckr**2
      cmuln=3.0_PREC
      snsa=-cskrx
      snsb=-cskrx/pckr-snkrx
      cneumn(1,ith,js)=snsa
      cneumn(2,ith,js)=snsb
      DO i=3,nranki
        snsc=cmuln*snsb/pckr-snsa
        cneumn(i,ith,js)=snsc
        snsa=snsb
        snsb=snsc
        cmuln=cmuln+2.0_PREC
      END DO
!     perform wronskian test on orders 0 and 1 and orders nrank-1 and nrank.
      quanbt=abs(ckr2*(bsslsp(2,ith,js)*cneumn(1,ith,js)-bsslsp(1,ith,js)*cneumn(2,ith,js))-1.0_PREC)
      quannt=abs(ckr2*(bsslsp(nranki,ith,js)*cneumn(nrank,ith,js)-bsslsp(nrank,ith,js)*cneumn(nranki,ith,js))-1.0_PREC)
      IF (quanbt-REAL(1.0e-10,PREC) >= 0.0_PREC .OR. quannt-REAL(1.0e-10,PREC) >= 0.0_PREC) THEN
        thtprt = rtd*theta
        IF (debug) THEN
          WRITE(6,'(/,10x,"theta=",f9.4,"kr=",f10.4,"bessel test=",e12.5,"neumann test=",e12.5)') &
               thtprt,pckrr,quanbt,quannt
        ENDIF
      END IF

    END IF
    
  END SUBROUTINE genbsl
!-----------------------------------------------------------------
  SUBROUTINE bessel(norder, argmnt, answr, ierror)
    
    IMPLICIT NONE
    REAL(PREC) :: argmnt, answr, acr, apr, ci, cn, cni, fact
    REAL(PREC) :: prod, sum, topr, x
    INTEGER    :: norder, ierror, ifct, n, i
    LOGICAL    :: lconverged

    ierror = 0
    n = norder
    x = argmnt
    cn = n
    sum = 1.0_PREC
    apr = 1.0_PREC
    topr = -0.5_PREC * x * x
    ci = 1.0_PREC
    cni = 2 * n + 3.0_PREC
    cni = 2 * n + 3.0_PREC

    lconverged = .FALSE.
    DO i = 1, 100
      acr = topr * apr / (ci * cni)
      sum = sum + acr
      IF (ABS(acr / sum) - REAL(1.0e-20,PREC) <= 0.0_PREC) THEN
        lconverged = .TRUE.
        EXIT
      END IF
      apr = acr
      ci = ci + 1.0_PREC
      cni = cni + 2.0_PREC
    END DO

    IF (.NOT. lconverged) THEN
      ierror = 1
    ELSE

!       the series has converged
      prod = 2 * n + 1.0_PREC
      fact = 1.0_PREC

      IF (n > 0) THEN ! 160, 160, 120
        DO ifct = 1, n
          fact = fact * x / prod
          prod = prod - 2.0_PREC
        END DO
      END IF
      answr = fact * sum

    END IF
    
  END SUBROUTINE bessel
!-----------------------------------------------------------------
      SUBROUTINE genbkr(xxr,xxi,iswt,ith,js,ierr)

!     Generate Bessel functions for complex arguments.
!     The index on the function is incremented by one.

      IMPLICIT NONE

      REAL(PREC) :: xxr, xxi
      INTEGER :: iswt, ith, js, ierr

      REAL(PREC) :: pckrr, pckri, rm, eposx, enegx, ex1, ex2, csr, csi, ar, ai
      REAL(PREC) :: f, cfr, cfi, temr, temi, tr, ti, pxr, pxi, ccr, cci
      REAL(PREC) :: quabt, quant, quanbt, quannt, pckr2r, pckr2i, thtprt, alarge
      INTEGER :: nval, ie, k, ij, ii, itest, i
      REAL(PREC) :: t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i, t5r, t5i

      REAL(PREC) :: bsr(41),bsi(41),cnr(41),cni(41),rjr(301),rji(301),pr(20),pi(20)
      INTEGER :: ipntr(20)
      LOGICAL :: lerr

      lerr = .false.

      alarge=REAL(1.0e30,PREC)
      pckrr=xxr
      pckri=xxi
      rm = cdabx(pckrr,pckri)
      nval=50
      IF ((rm > 25.0_PREC  ).AND.(rm <= 150.0_PREC  )) THEN
        nval=2*rm
      ELSE IF (rm > 150.0_PREC  ) THEN
        nval=300
      END IF

!     generate bessel functions.
!
      rjr(nval+1)=0.0_PREC
      rji(nval+1)=0.0_PREC
      rji(nval)=0.0_PREC
      rjr(nval)=1.0_PREC
      IF (rm > 25.0_PREC  ) THEN
        rjr(nval) = REAL(1.0e-30,PREC)
      ELSE IF (rm > 10.0_PREC  ) THEN
        rjr(nval) = REAL(1.0e-20,PREC)
      ELSE IF (rm > 2.0_PREC  ) THEN
        rjr(nval) = REAL(1.0e-10,PREC)
      END IF
      
      ie=nval+2
      eposx=exp(pckri)
      enegx=exp(-pckri)
      ex1 = (enegx+eposx)/2.0_PREC
      ex2 = (enegx-eposx)/2.0_PREC
      csr =  sin(pckrr)*ex1
      csi = - cos(pckrr)*ex2
      ar = cddvdr(csr,csi,pckrr,pckri,lerr)
      if (lerr) then
        ierr=80001
        print '(A,I5,A)', 'ERR #',ierr,': div0 in ar/i calc in genbkr - RETURNING.'
        return
      end if
      ai = cddvdi(csr,csi,pckrr,pckri)
      k=0

!      print *,'genbkr 2'

      DO i=2,nval
        ij=ie-i
        f=REAL(2*ij-1,PREC)
        cfr=cddvdr(f,0._PREC,pckrr,pckri,lerr)
        if (lerr) then
          ierr=10100+i
          print '(A,I5,A,I1,A)', &
            'ERR #',ierr,': div0 at val=',i,' cfr/i calc in loop 10 of genbkr - RETURNING.'
          return
        end if
        cfi=cddvdi(f,0._PREC,pckrr,pckri)
        rjr(ij-1) = cdmpyr(rjr(ij),rji(ij),cfr,cfi)-rjr(ij+1)
        rji(ij-1) = cdmpyi(rjr(ij),rji(ij),cfr,cfi)-rji(ij+1)
        temr= abs(rjr(ij-1))
        temi= abs(rji(ij-1))
! JM190724:
! Seemed very fishy to loop and do nothing calc-wise until values
! exceeded a very large number. Which also led to k remaining 0, hence
! code crashing when accessing ipntr(k) below.
! FIXED: turned the check around, which made the code run through. and
!        provide correct results (apart from a 1d2 factor).
!        if (temr.le.alarge.and.temi.le.alarge)  go to 10
        IF (temr <= alarge .OR. temi <= alarge)  THEN
!
!******* bump up the pointer k and store the current value
!******* of rj at p(k), and its position at ipntr(k)
!******* (taking care of real and im. parts)
!
          ii=ij+1
          IF (ii > nranki) THEN
!
!******* (ij+1).gt.nranki  so normalize without storing
!******* rj and its position
!
            tr=cddvdr(rjr(ij),rji(ij),rjr(ij-1),rji(ij-1),lerr)
            IF (lerr) THEN
              ierr=10100+i
              PRINT '(A,I5,A,I1,A)', &
                   'ERR #',ierr,': div0 at val=',i,' tr/i calc in loop 10 of genbkr - RETURNING.'
              RETURN
            END IF
            ti=cddvdi(rjr(ij),rji(ij),rjr(ij-1),rji(ij-1))
            rjr(ij)=tr
            rji(ij)=ti
            rjr(ij-1)=1.0_PREC
            rji(ij-1)=0.0_PREC
            
          ELSE
!
!******* (ij+1).le.nranki  so normalize and store rj and its position
!
            k=k+1
            IF (k > 20) THEN
              IF (debug) THEN
                WRITE (6,*) ' pr,pi,ipntr arrays are too small'
              ENDIF
              STOP
            ELSE
              pr(k)=rjr(ij-1)
              pi(k)=rji(ij-1)
              ipntr(k)=ij+1
              rjr(ij-1)=1.0_PREC
              rji(ij-1)=0.0_PREC
              tr=cddvdr(rjr(ij),rji(ij),pr(k),pi(k),lerr)
              IF (lerr) THEN
                ierr=17100+i
                PRINT '(A,I5,A,I1,A)', &
                     'ERR #',ierr,': div0 at val=',i,' tr/i calc in loop 7000/10 of genbkr - RETURNING.'
                RETURN
              END IF
              ti=cddvdi(rjr(ij),rji(ij),pr(k),pi(k))
              rjr(ij)=tr
              rji(ij)=ti
              ii=ipntr(k)-2
            END IF
          END IF

        END IF
        
      END DO

!
!******* backsolution
!
      pxr = cddvdr(ar,ai,rjr(1),rji(1),lerr)
      if (lerr) then
        ierr=80002
        print '(A,I5,A)', 'ERR #',ierr,': div0 in pxr/i calc in genbkr - RETURNING.'
        return
      end if
      pxi = cddvdi(ar,ai,rjr(1),rji(1))
      DO i=1,nranki
! JCS - This is going to try ipntr(0) when k=0, which is out of bounds. Haven't been able to track down cause/fix
!      print *,'k=',k,' i=',i,' ipntr(k)=',ipntr(k)
        IF (ipntr(k) /= i) THEN
          bsr(i)=cdmpyr(rjr(i),rji(i),pxr,pxi)
          bsi(i)=cdmpyi(rjr(i),rji(i),pxr,pxi)

        ELSE
!
!******* make correction for px  ...  divide   by p(k) (real and im)
!
          tr=cddvdr(pxr,pxi,pr(k),pi(k),lerr)
          IF (lerr) THEN
            ierr=15100+i
            PRINT '(A,I5,A,I1,A)', &
                 'ERR #',ierr,': div0 at val=',i,' tr/i calc in loop 5001/10 of genbkr - RETURNING.'
            RETURN
          END IF
          ti=cddvdi(pxr,pxi,pr(k),pi(k))
          pxr=tr
          pxi=ti
          k=k-1
          bsr(i)=cdmpyr(rjr(i),rji(i),pxr,pxi)
          bsi(i)=cdmpyi(rjr(i),rji(i),pxr,pxi)
        END IF
      END DO

!     generate neumann functions for test.
!
      ccr=cos(pckrr)*ex1
      cci=sin(pckrr)*ex2
      cnr(1)=-cddvdr(ccr,cci,pckrr,pckri,lerr)
      if (lerr) then
        ierr=80003
        print '(A,I5,A)', 'ERR #',ierr,': div0 in cnr/i calc in genbkr - RETURNING.'
        return
      end if
      cni(1)=-cddvdi(ccr,cci,pckrr,pckri)
      cnr(2)=cddvdr(cnr(1),cni(1),pckrr,pckri,lerr)-ar
      cni(2)=cddvdi(cnr(1),cni(1),pckrr,pckri)-ai

      DO i=3,nranki
        f=REAL(2*i-3,PREC)
        cfr=cddvdr(f,0._PREC,pckrr,pckri,lerr)
        cfi=cddvdi(f,0._PREC,pckrr,pckri)
        cnr(i)=cdmpyr(cnr(i-1),cni(i-1),cfr,cfi)-cnr(i-2)
        cni(i)=cdmpyi(cnr(i-1),cni(i-1),cfr,cfi)-cni(i-2)
      END DO

      IF(iswt == 1) THEN ! go to 101

        DO i=1,nranki
          bslkpr(i,ith,js)=bsr(i)
          bslkpi(i,ith,js)=bsi(i)
          cneumr(i,ith,js)=cnr(i)
          cneumi(i,ith,js)=cni(i)
        END DO

      ELSE

        IF (iswt == 2) THEN
          DO i=1,nranki
            bkpr2(i,ith,js)=bsr(i)
            bkpi2(i,ith,js)=bsi(i)
            cnpr2(i,ith,js)=cnr(i)
            cnpi2(i,ith,js)=cni(i)
          END DO
        ELSE
          DO i=1,nranki
            bkprr2(i,ith,js)=bsr(i)
            bkpii2(i,ith,js)=bsi(i)
          END DO
        END IF
      END IF

!     Perform two tests on Bessel and Neumann functions.  First test is
!     most accurate for large arguments and the second is most accurate
!     for smaller arguments.  If either test is passed, functions are good.

!     For large arguments abs(bessel) should equal abs(neumann).

      
      IF ( (MAX(cnr(1),cni(1)) .EQ. 0d0) .OR. &
           (max(cnr(nranki),cni(nranki)) .eq. 0d0) ) then
        ierr=80004
        PRINT '(A,I5,A)', 'ERR #',ierr,': div0 in qua[b/n]t calc in genbkr - RETURNING.'
        PRINT *, cnr(1),cni(1),cnr(nranki),cni(nranki)
        return
      END IF
      
      quabt=cdabx(bsr(1),bsi(1) )/cdabx(cnr(1),cni(1))-1.0_PREC
      quant=cdabx(bsr(nranki),bsi(nranki))/cdabx(cnr(nranki),cni(nranki))-1.0_PREC
      IF(( ABS(quabt) <= REAL(1.0d-8,PREC)) .OR. ( ABS(quant) <= REAL(1.0d-8,PREC))) THEN

        !     Perform Wronskian test if large argument test fails.
        pckr2r = pckrr**2-pckri**2
        pckr2i = 2.0_PREC  *pckrr*pckri

        !     bessel test
        t1r=cdmpyr(bsr(2),bsi(2),cnr(1),cni(1))
        t1i=cdmpyi(bsr(2),bsi(2),cnr(1),cni(1))
        t2r=cdmpyr(bsr(1),bsi(1),cnr(2),cni(2))
        t2i=cdmpyi(bsr(1),bsi(1),cnr(2),cni(2))
        t3r = cddvdr(t1r,t1i,t2r,t2i,lerr)-1.0_PREC
        IF (lerr) THEN
          ierr=80004
          PRINT '(A,I5,A)', 'ERR #',ierr,': div0 in t3r/i calc in genbkr - RETURNING.'
          RETURN
        END IF
        t3i = cddvdi(t1r,t1i,t2r,t2i)
        t4r = cdmpyr(t2r,t2i,t3r,t3i)
        t4i = cdmpyi(t2r,t2i,t3r,t3i)
        t5r = cdmpyr(pckr2r,pckr2i,t4r,t4i)-1.0_PREC
        t5i = cdmpyi(pckr2r,pckr2i,t4r,t4i)
        
        quanbt = cdabx(t5r,t5i)

        !     neumann test
        t1r=cdmpyr(bsr(nranki),bsi(nranki),cnr(nrank),cni(nrank))
        t1i=cdmpyi(bsr(nranki),bsi(nranki),cnr(nrank),cni(nrank))
        t2r=cdmpyr(bsr(nrank),bsi(nrank),cnr(nranki),cni(nranki))
        t2i=cdmpyi(bsr(nrank),bsi(nrank),cnr(nranki),cni(nranki))
        t3r = cddvdr(t1r,t1i,t2r,t2i,lerr)-1.0_PREC
        IF (lerr) THEN
          ierr=80004
          PRINT '(A,I5,A)', 'ERR #',ierr,': div0 in t3r/i calc in genbkr - RETURNING.'
          RETURN
        END IF
        t3i = cddvdi(t1r,t1i,t2r,t2i)
        t4r = cdmpyr(t2r,t2i,t3r,t3i)
        t4i = cdmpyi(t2r,t2i,t3r,t3i)
        t5r = cdmpyr(pckr2r,pckr2i,t4r,t4i)-1.0_PREC
        t5i = cdmpyi(pckr2r,pckr2i,t4r,t4i)
 
        quannt = cdabx(t5r,t5i)
        IF((quanbt > REAL(1.0e-08,PREC)) .OR. (quannt > REAL(1.0e-08,PREC))) THEN

          thtprt = rtd*theta
          !      if (debug) then
          !        write(6,50) thtprt,pckrr,pckri,quabt,quant,quanbt,quannt
          !   50 format(//,2x,8h theta =,f9.4,6h, kr =,2f10.4,7h,abt =,f12.5,7h,a
          !     1nt =,e12.5,7h, nbt =,e12.5,7h, nnt =,e12.5)
          !      endif
          
        END IF
      
      END IF

    END SUBROUTINE genbkr
!-----------------------------------------------------------------
    REAL(PREC) FUNCTION crootr(a, b)

      REAL(PREC) :: a, b, dmag, angle

      dmag = (a * a + b * b)**0.25_PREC
      angle = 0.5_PREC * ATAN2(b, a)
      crootr = dmag * COS(angle)
      
    END FUNCTION crootr
!-----------------------------------------------------------------
    REAL(PREC) FUNCTION crooti(a, b)

      REAL(PREC) :: a, b, dmag, angle

      dmag = (a * a + b * b)**0.25_PREC
      angle = 0.5_PREC * ATAN2(b, a)
      crooti = dmag * SIN(angle)
    END FUNCTION crooti
!-----------------------------------------------------------------
    SUBROUTINE prcssm(af,bf,nr,nri,ierr)

!     A routine to solve the equation t = (a-inverse)*b  ( all matrices
!     are transposed) using Gauss-Jordan elimination.

      IMPLICIT NONE

      REAL(PREC) :: af, bf, aijmax, t1, t2, arat
      INTEGER :: nr, nri, ierr, n, i, j, k, jmax, fl(40), jj
      LOGICAL :: lerr
      DIMENSION af(40,40,2),bf(40,40,2)
      DIMENSION aijmax(2),arat(2)

      lerr = .false.

      n = 2*nr

      !     start reduction of the a matrix.
      DO i = 1,n
        !     search for the maximum element in the ith row of the a-matrix.
        aijmax(1) = af(i,1,1)
        aijmax(2) = af(i,1,2)
        jmax = 1
        DO j = 2,n
          IF(cdabx(af(i,j,1),af(i,j,2)) > cdabx(aijmax(1),aijmax(2))) THEN
            aijmax(1) = af(i,j,1)
            aijmax(2) = af(i,j,2)
            jmax = j
          END IF
        END DO
          !     if aijmax is zero ( as it will be for any row (or column) where the
!     index m is .gt. the index n, i.e., the legendre functions force those
!     matrix elements to zero),then the matrix is singular so solve the
!     reduced matrix (order = 2*(nrank-m)).
        IF(cdabx(aijmax(1),aijmax(2)) > 0.0_PREC  ) THEN
!     normalize the ith row by aijmax (jmax element of the ith row).
          DO j = 1,n
            t1 = af(i,j,1)
            t2 = af(i,j,2)
            af(i,j,1) = cddvdr(t1,t2,aijmax(1),aijmax(2),lerr)
            IF (lerr) THEN
              ierr=70001
              PRINT '(A,I5,A)', 'ERR #',ierr,': div0 in af calc in prcssm - RETURNING.'
              RETURN
            END IF
            af(i,j,2) = cddvdi(t1,t2,aijmax(1),aijmax(2))
            !     normalize the ith row of b.
            t1 = bf(i,j,1)
            t2 = bf(i,j,2)
            bf(i,j,1) = cddvdr(t1,t2,aijmax(1),aijmax(2),lerr)
            bf(i,j,2) = cddvdi(t1,t2,aijmax(1),aijmax(2))
          END DO
          !     use row transformations to get zeros above and below the jmax
          !     element of the ith row of a.  apply same row transformations
          !     to the b matrix.
          loop_70: DO k = 1,n
            IF (k == i) CYCLE loop_70
            arat(1) = -af(k,jmax,1)
            arat(2) = -af(k,jmax,2)
            loop_50: DO j = 1,n
              IF(cdabx(af(i,j,1),af(i,j,2)) > 0.0_PREC  ) THEN
                af(k,j,1) = cdmpyr(arat(1),arat(2),af(i,j,1),af(i,j,2))+af(k,j,1)
                af(k,j,2) = cdmpyi(arat(1),arat(2),af(i,j,1),af(i,j,2))+af(k,j,2)
              END IF
            END DO loop_50
            af(k,jmax,1)=0.0_PREC
            af(k,jmax,2)=0.0_PREC
            loop_60: DO j=1,n
              IF(cdabx(bf(i,j,1),bf(i,j,2)) > 0.0_PREC  ) THEN
                bf(k,j,1) = cdmpyr(arat(1),arat(2),bf(i,j,1),bf(i,j,2))+bf(k,j,1)
                bf(k,j,2) = cdmpyi(arat(1),arat(2),bf(i,j,1),bf(i,j,2))+bf(k,j,2)
              END IF
            END DO loop_60
          END DO loop_70
        ELSE
          jmax = i
        END IF

        !     store row counter (i) in top element of jmax column.  thus,
        !     the vector will contain the location of the pivot
        !     (unity) element of each column (after reduction).
        fl(jmax) = i
      END DO
      !      write(*,*) 'after line 80 prcssm'
      !     the reduction of a is complete.  perform row interchanges
      !     as indicated in the first row of a.
      loop_120: DO i = 1,n
        k=i
        !     put the integer value in fl into k.
        DO 
          k = fl(k)
          IF (k-i >= 0) EXIT
        END DO
        IF (k > i) THEN
          !     if k(1,i) is less than i, then that row has already been
          !     involved in an interchange, and we use k(1,k) until we get
          !     a value of k greater than i (corresponding to a row stored
          !     below the ith row).
          DO j=1,n
            arat(1) = bf(i,j,1)
            arat(2) = bf(i,j,2)
            bf(i,j,1) = bf(k,j,1)
            bf(i,j,2) = bf(k,j,2)
            bf(k,j,1) = arat(1)
            bf(k,j,2) = arat(2)
          END DO
        END IF
      END DO loop_120

    END SUBROUTINE prcssm

!-----------------------------------------------------------------

    SUBROUTINE addprc(f_a_q,f_a0_q,f_b_q,f_b0_q)

!     A routine to obtain the scattered field coefficients and calculate
!     the differential scattering cross section in the azimuthal plane.

      IMPLICIT NONE

      COMPLEX(PREC), INTENT(out) :: f_a_q,f_a0_q,f_b_q,f_b0_q

      COMPLEX(PREC), PARAMETER :: ci = (0.0_PREC,1.0_PREC)
      COMPLEX(PREC) :: cim
      REAL(PREC) :: cn, ci1r, ci1i, ci2r, ci2i, p1, p2
      REAL(PREC) :: s1r, s1i, s2r, s2i, sum1, sum2, temp1, temp2
      REAL(PREC) :: cnorm, cir, cii, f1r, f1i, g1r, g1i, f2r, f2i, g2r, g2i
      REAL(PREC) :: extpp, extper
      REAL(PREC) :: forrp, forip, forpe, foripe, borrp, borip, borrpe, boripe
      REAL(PREC) :: xhor, yver
      INTEGER :: n, np, n1, nr2, i, ii, iu, iup, j

      REAL(PREC)::  zxold(181),zyold(181),ab1(40,2),ab2(40,2),fg1(40,2),fg2(40,2),fgans(181,2,2)

      tmat=x

!     Generate the legendre functions for the incident angle.
      IF (anginc == 0.0_PREC) THEN
        costh=1.0_PREC
        sinth=0.0_PREC
        theta=0.0_PREC
      ELSE
        IF (anginc == 180.0_PREC) THEN
          costh=-1.0_PREC
          sinth=0.0_PREC
          theta=0.0_PREC
        ELSE
          theta = dtr*anginc
          sinth=SIN(theta)
          costh=COS(theta)
        END IF
      END IF
      
      CALL genlgp (kmv, twm, theta, sinth, costh, pnmllg)


!     Generate the incident field coefficients -- ab1 = theta polarization
!     and ab2 = phi polarization.
      DO n=1, nrank
        np = n + nrank
        cn=REAL(n,PREC)
        n1 = n+1
        ci1r=REAL(ci**n)
        ci1i=AIMAG(ci**n)
        ci2r=REAL(ci**n1)
        ci2i=AIMAG(ci**n1)
        p1 = cn*costh*pnmllg(n1)-(cn+cmv)*pnmllg(n)
        p2 = cmv*pnmllg(n1)
        ab1(n,1) = -ci1r*p2
        ab1(n,2) = -ci1i*p2
        ab1(np,1) = ci2r*p1
        ab1(np,2) = ci2i*p1
        ab2(n,1) = ci1r*p1
        ab2(n,2) = ci1i*p1
        ab2(np,1) = -ci2r*p2
        ab2(np,2) = -ci2i*p2
      END DO

!     The scattered field coefficients = the row vector of incident field
!     coefficients times the t-transposed matrix.

      nr2 = 2*nrank
      DO j = 1,nr2
        s1r=0.0_PREC
        s1i=0.0_PREC
        s2r=0.0_PREC
        s2i=0.0_PREC
        DO i = 1,nr2
          s1r = s1r+cdmpyr(ab1(i,1),ab1(i,2),tmat(i,j,1),tmat(i,j,2))
          s1i = s1i+cdmpyi(ab1(i,1),ab1(i,2),tmat(i,j,1),tmat(i,j,2))
          s2r = s2r+cdmpyr(ab2(i,1),ab2(i,2),tmat(i,j,1),tmat(i,j,2))
          s2i = s2i+cdmpyi(ab2(i,1),ab2(i,2),tmat(i,j,1),tmat(i,j,2))
        END DO
        fg1(j,1) = s1r
        fg1(j,2) = s1i
        fg2(j,1) = s2r
        fg2(j,2) = s2i
      END DO

!     Calculate scattering cossections normalized for parallel and
!     perpendicular polz
      sum1=0.0_PREC
      sum2=0.0_PREC
      DO i=1,nrank
        ii=i+nrank
        temp1=fg1(i,1)**2+fg1(i,2)**2+fg1(ii,1)**2+fg1(ii,2)**2
        temp2=fg2(i,1)**2+fg2(i,2)**2+fg2(ii,1)**2+fg2(ii,2)**2
        sum1=sum1+temp1/cmxnrm(i)
        sum2=sum2+temp2/cmxnrm(i)
      END DO
!     normalize scattering crossections.
      sum1=(rtsfct*2.0_PREC/conk)*sum1
      sum2=(rtsfct*2.0_PREC/conk)*sum2
!         normalize w.r.t. eq. spherical dia.
          cnorm=aovrb**(-2._PREC/3._PREC)
          sum1=sum1*cnorm
          sum2=sum2*cnorm
!     accumulate results for each m value
      summ1=sum1+summ1
      summ2=sum2+summ2

!     Evaluate the scattered field at each scattering angle.

      DO iu = 1,nuang
!     Generate the Legendre multipliers.

        IF (uang(iu) == 0.0_PREC) THEN
          costh=1.0_PREC
          sinth=0.0_PREC
          theta=0.0_PREC
        ELSE
          IF (uang(iu) == 180.0_PREC) THEN
            costh=-1.0_PREC
            sinth=0.0_PREC
            theta=0.0_PREC
          ELSE
            theta = dtr*uang(iu)
            sinth=SIN(theta)
            costh=COS(theta)
          END IF
        END IF
        
        CALL genlgp (kmv, twm, theta, sinth, costh, pnmllg)

        fgans(iu,1,1)=0.0_PREC
        fgans(iu,1,2)=0.0_PREC
        fgans(iu,2,1)=0.0_PREC
        fgans(iu,2,2)=0.0_PREC
        DO n = 1,nrank
          np = n+nrank
          n1 = n+1
          cn=REAL(n,PREC)
          p1 = cn*costh*pnmllg(n1)-(cn+cmv)*pnmllg(n)
          p2 = cmv*pnmllg(n1)
          cim = (-ci)**n1
          cir=REAL(cim)
          cii=AIMAG(cim)
          f1r = fg1(n,1)*p2
          f1i = fg1(n,2)*p2
          g1r = -fg1(np,2)*p1
          g1i = fg1(np,1)*p1
          fgans(iu,1,1) = fgans(iu,1,1)+cdmpyr(cir,cii,f1r+g1r,f1i+g1i)/cmxnrm(n)
          fgans(iu,1,2) = fgans(iu,1,2)+cdmpyi(cir,cii,f1r+g1r,f1i+g1i)/cmxnrm(n)
          f2r = fg2(n,1)*p1
          f2i = fg2(n,2)*p1
          g2r = -fg2(np,2)*p2
          g2i = fg2(np,1)*p2
          fgans(iu,2,1) = fgans(iu,2,1)-cdmpyr(cir,cii,f2r+g2r,f2i+g2i)/cmxnrm(n)
          fgans(iu,2,2) = fgans(iu,2,2)-cdmpyi(cir,cii,f2r+g2r,f2i+g2i)/cmxnrm(n)
        END DO

        !     The normalized diff.scat.cross sect. is given by ((8/ka)*fgans)**2
        !     scale fgans to calculate diff. scat. cross sect. (rtsfct = 8/ka)

        fgans(iu,1,1) = rtsfct*fgans(iu,1,1)
        fgans(iu,1,2) = rtsfct*fgans(iu,1,2)
        fgans(iu,2,1) = rtsfct*fgans(iu,2,1)
        fgans(iu,2,2) = rtsfct*fgans(iu,2,2)
      END DO

!     Accumulate the results for each m value.

      DO iup = 1,nuang
        acans(iup,1,1) = acans(iup,1,1)+fgans(iup,1,1)
        acans(iup,1,2) = acans(iup,1,2)+fgans(iup,1,2)
        acans(iup,2,1) = acans(iup,2,1)+fgans(iup,2,1)
        acans(iup,2,2) = acans(iup,2,2)+fgans(iup,2,2)
      END DO

!     Calculate the extinction crossections.
      extpp=acans(1,1,2)*rtsfct/4.0_PREC
      extper=acans(1,2,2)*rtsfct/4.0_PREC

!     Normalize wrt equivalent spherical diameter
      extpp=extpp*cnorm
      extper=extper*cnorm

!     Calculate forward and backward amplitude in far zone.
!     sigma equals 4.0/k
      forrp=sigma*acans(1,1,1)/rtsfct
      forip=sigma*acans(1,1,2)/rtsfct
      forpe=sigma*acans(1,2,1)/rtsfct
      foripe=sigma*acans(1,2,2)/rtsfct
      borrp=sigma*acans(nuang,1,1)/rtsfct
      borip=sigma*acans(nuang,1,2)/rtsfct
      borrpe=sigma*acans(nuang,2,1)/rtsfct
      boripe=sigma*acans(nuang,2,2)/rtsfct

!     Calculate normalized radar crossections for both polarizations
      xhor=acans(nuang,1,1)**2+acans(nuang,1,2)**2
      yver=acans(nuang,2,1)**2+acans(nuang,2,2)**2

!     Normalize wrt equivalent spherical diameter
      xhor= xhor*cnorm
      yver= yver*cnorm

!     Store the scattering results

      f_a_q=cmplx(borrp,borip,PREC)
      f_a0_q=cmplx(forrp,forip,PREC)
      f_b_q=cmplx(borrpe,boripe,PREC)
      f_b0_q=cmplx(forpe,foripe,PREC)

      end subroutine addprc

!-----------------------------------------------------------------
      subroutine genkr

! JCS - Try quad precision

!       calculate ckr and dckr as a function of theta for a oblate spheroid
        IMPLICIT NONE

        REAL(PREC) :: bovra, qb

        bovra = 1.0_PREC / aovrb
        qb = 1.000_PREC / sqrt((bovra * costh)**2 + sinth**2)
        ckr = conk * qb
        dckr = conk * costh * sinth * (bovra**2 - 1.000_PREC) * qb**3
!       thet1=theta*rtd

      end subroutine genkr
!-----------------------------------------------------------------
      subroutine genkr2

        IMPLICIT NONE

!       calculate ckr2 and dckr2 for sph-cone-oblate

        REAL(PREC) :: bovra2, qb2

        bovra2 = 1.0_PREC / aovrb2

!       inner particle is an oblate spheroid
        qb2 = 1._PREC / sqrt((bovra2 * costh)**2 + sinth**2)
        ckr2 = conk2 * qb2
        dckr2 = conk2 * costh * sinth * (bovra2**2 - 1.0_PREC) * qb2**3

!       inner particle is a sphere
!       ckr2=conk2
!       dckr2=0.0
!       thet1=theta*rtd
! 10     format(2x,3(e15.7,4x))

      end subroutine genkr2
!----------------------------------------------------------------
      subroutine calenp
        IMPLICIT NONE

        epps(1) = 0._PREC
        epps(2) = cpi / 2._PREC

      end subroutine calenp

END MODULE radar_dualpol_t_matrix2_mod
