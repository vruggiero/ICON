!
!+ icosahedral grid, routines adapted from GME
!
MODULE mo_ico_grid
!
! Description:
!   Adapted routines from GME for the icosahedral grid
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
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
! V1_28        2014/02/26 Harald Anlauf
!  Protect Intel-specific compiler directives by #ifdef
! V1_45        2015-12-15 Harald Anlauf
!  OpenMP parallelize loop nest
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
!------------------------------------------------------------------------------

  !===========================
  ! icosahedral-hexagonal grid
  !===========================
  !-------------
  ! Modules used
  !-------------
  USE mo_kind,       ONLY: wp
  USE mo_exception,  ONLY: finish
  USE mo_atm_decomp, ONLY: t_atm_dec
#ifndef NOMPI
  USE mo_mpi_dace,   ONLY: dace,         &
                           p_min, p_max, &
                           p_real,       &
                           MPI_STATUS_SIZE
#endif

  IMPLICIT NONE

  !----------------
  ! public entities
  !----------------
  private
  public :: factorize_nir
  public :: global_coordinates
  public :: nspoke2

  !===========================================================================
  !
  ! Valid for all grids

  !--------------------------------------------------
  ! Normal spokes in two fashions (direct neighbours)
  !--------------------------------------------------
  INTEGER, PARAMETER :: &
       nsp11 =  1, nsp12 =  0, nsp13 = -1,  &
       nsp14 = -1, nsp15 =  0, nsp16 =  1,  &
       nsp21 =  0, nsp22 =  1, nsp23 =  1,  &
       nsp24 =  0, nsp25 = -1, nsp26 = -1

  INTEGER, PARAMETER :: nspoke (12) = &
       (/ 1,  0, -1, -1,  0,  1,  0,  1,  1,  0, -1, -1 /)

  !------------------------------
  ! same with 2-dimensional shape
  !------------------------------
  Integer, Parameter :: &
    nspoke2(2,6) =&    ! Offsets of 6(5) neighbours
                       ! [1:6,j1/j2]
    Reshape(Source = (/ 1,  0, -1, -1,  0,  1,     &
                        0,  1,  1,  0, -1, -1 /),  &
            Shape  = (/ 2, 6 /), Order = (/ 2, 1 /) )

  !------------------------------------------------
  ! Spokes for first and second grid distances away
  !------------------------------------------------
  INTEGER, PARAMETER :: nspokes (12,4) = RESHAPE ( &
       SOURCE = (/  1,  0, -1, -2,  0,  1,  0,  1,  2,  2, -1, -1,    &
                    2,  0, -1, -1,  0,  1, -2,  1,  1,  0, -1, -2,    &
                    1,  0, -1, -2, -2, -1,  0,  1,  1,  1,  0, -1,    &
                    1,  1, -1, -1,  0,  1,  0,  1,  1,  0, -1, -1 /), &
                SHAPE = (/ 12, 4 /) )

#ifndef NOMPI
  INTEGER :: mpi_err
#endif

  !===========================================================================
  !
  ! Specific for each grid level

  TYPE triangular_grid
    type (t_atm_dec) :: dc
    !================================
    ! Base grid description parameter
    !================================
    INTEGER :: ni
    !---------------------------------------------------------------
    ! Factorize ni into 3**ni3 * 2**ni2 with ni3 0 or 1 and ni2 > 1,
    !---------------------------------------------------------------
    INTEGER :: ni2, nir
    !------------------------------------------------------------
    ! Number of grids allowed for multigrid grid level generation
    !------------------------------------------------------------
    INTEGER :: ni2_mg
    !------------------------------------------------------------
    ! Grid generation level for later use of multigrid grid level
    ! generation
    !------------------------------------------------------------
    INTEGER, POINTER :: nlev_mg (:,:,:)
    !-------------------
    ! Dimensions of grid
    !-------------------
    INTEGER :: ng1s, ng1e, ng1sm1, ng1sm2, ng1ep1, ng1ep2
    INTEGER :: ng2s, ng2e, ng2sm1, ng2sm2, ng2ep1, ng2ep2
    INTEGER :: ng3s, ng3e
    INTEGER :: ngg1s, ngg1e, ngg2s, ngg2e
    !-------------------
    ! Number of diamonds
    !-------------------
    INTEGER :: nd
    !----------------------------------------------------------
    ! First and second dimension indices of the mirrored points
    !----------------------------------------------------------
    INTEGER :: ni1mrp(8), ni2mrp(8)
    !-----------------------------------------------------------------
    ! Cartesian (x,y,z) coordinates of the location vector of the grid
    ! points (nodes) on the unit-sphere; phys. dim. ( - )
    ! Full domain of the diamonds
    !
    ! rank 1: i index on subdomain (ith diamond) [0:ni]
    ! rank 2: j index on subdomain (ith diamond) [1:ni+1]
    ! rank 3: x, y, z coordinates                [1:3]
    ! rank 4: number of diamond                  [1:10]
    !-----------------------------------------------------------------
    REAL(wp), POINTER :: xnglob(:,:,:,:)
    REAL(wp), POINTER :: xn(:,:,:,:)
    !---------------------------------------------------------------
    ! geographical longitude of the gridpoints; phys. dim. (radians)
    !---------------------------------------------------------------
    REAL(wp), POINTER :: rlon(:,:,:)
    !--------------------------------------------------------------
    ! geographical latitude of the gridpoints; phys. dim. (radians)
    !--------------------------------------------------------------
    REAL(wp), POINTER :: rlat(:,:,:)
    !-------------------
    ! Coriolis Parameter
    !-------------------
    REAL(wp), POINTER :: corio(:,:,:)
    !------------------------------------------------------------------------
    ! area (excess) of triangles (l'Hulier - multiply by Re**2 for real area)
    !------------------------------------------------------------------------
    REAL(wp), POINTER :: area(:,:,:)
    !------------------------------------------------------------------
    ! reciprocal of the areas of the hexagon related to the gridpoints;
    ! phys. dim. ( - )
    !------------------------------------------------------------------
    REAL(wp), POINTER :: rarn(:,:)
    !--------------------------------------------------------------------
    ! weights for the calculation of diagnostic quantities (mean values);
    ! phys. dim. ( - )
    !--------------------------------------------------------------------
    REAL(wp), POINTER :: hexwgt(:,:)
    !-----------------------------------------------------------------------
    ! Cartesian (x,y) coordinates of the unit vector into the local easterly
    ! (northerly) direction on the unit-sphere (the z component is zero);
    !-----------------------------------------------------------------------
    REAL(wp), POINTER :: erlon(:, :, :, :)
    REAL(wp), POINTER :: erlat(:, :, :, :)
    !---------------------------------------------------------------------
    ! local longitude of the 6 (5) neighbouring gridpoints relative to the
    ! local (xn, erlon, erlat) system; phys. dim. (radians)
    !---------------------------------------------------------------------
    REAL(wp), POINTER :: eta(:, :, :)
    !--------------------------------------------------------------------------
    ! local latitude of the 6 (5) neighbouring gridpoints relative to the local
    ! (xn, erlon, erlat) system; phys. dim. (radians)
    !--------------------------------------------------------------------------
    REAL(wp), POINTER :: chi(:, :, :)
    !------------------------------------------------------------------------
    ! cosine of the rotation angle between the local coordinate system of the
    ! home node and the 6 (5) neighbouring gridpoints; phys. dim. ( - )
    !------------------------------------------------------------------------
    REAL(wp), POINTER :: cpsi(:, :, :)
    !----------------------------------------------------------------------
    ! sine of the rotation angle between the local coordinate system of the
    ! home node and the 6 (5) neighbouring gridpoints; phys. dim. ( - )
    !----------------------------------------------------------------------
    REAL(wp), POINTER :: spsi(:, :, :)
    !------------------------------------------------------------------
    ! factor 1/det needed to compute the barycentric coordinates of the
    ! departure or midpoint of the trajectory; phys. dim. (radians**2)
    !------------------------------------------------------------------
    REAL(wp), POINTER :: bary(:, :, :)
    !---------------------------------------------------------------------
    ! coefficients needed for the calculation of the gradients in eta- and
    ! chi-direction; phys. dim. (1/m)
    !---------------------------------------------------------------------
    REAL(wp), POINTER :: grd(:, :, :, :)
    !---------------------------------------------------------------------
    ! coefficients needed for the calculation of the Laplacian in eta- and
    ! chi-direction; phys. dim. ( - )
    !---------------------------------------------------------------------
    REAL(wp), POINTER :: rlap(:, :, :)
    !-----------------------------------------------
    ! minimum/maximum mesh width and triangle height
    !-----------------------------------------------
    REAL :: dxmin, dxmax, dhmin, dhmax

  END TYPE triangular_grid

#ifdef DEBUG
  CHARACTER(len=10) :: field_name
#endif

!=============================================================================

CONTAINS

  !===========================================================================
  SUBROUTINE generate_grid (grid, dc, kni)

!    USE mo_atm_decomp,     ONLY: setup_parallel,       &
!                                 domain_decomposition, &
!                                 setup_communication

    TYPE (triangular_grid) ,INTENT(inout) :: grid
    TYPE (t_atm_dec)       ,INTENT(in)    :: dc
    integer                ,INTENT(in)    :: kni

    INTEGER :: ierror = 0    ! Error flag; 0 if no error occured
    !------------------------------------------------------------
    ! 1. Set grid dimensions and allocate grid information fields
    !------------------------------------------------------------
    ! Number of diamonds (either 5 or 10, only 10 tested ...)
    !--------------------------------------------------------
    grid% nd = 10
    !-----------------------------------------------------------
    ! Factorize ni into 3**i3 * 2**i2 with i3 0 or 1 and i2 > 1,
    ! if possible, otherwise abort the program
    !-----------------------------------------------------------
    grid% ni = kni

    if (grid% ni < 4) call finish ('construct_atm_grid','ni must be >= 4')
    call factorize_nir (grid% ni, grid% ni2, grid% nir, ierror)
    if (ierror/=0) call finish ('construct_atm_grid','cannot factorize')

    !---------------------
    ! Set boundary indices
    !---------------------
    grid% ng1s = 0
    grid% ng1e = grid% ni

    grid% ng2s = 1
    grid% ng2e = grid% ni + 1

    grid% ng3s = 1
    grid% ng3e = 1

    grid% ng1sm1 = grid% ng1s - 1
    grid% ng1sm2 = grid% ng1s - 2
    grid% ng1ep1 = grid% ng1e + 1
    grid% ng1ep2 = grid% ng1e + 2

    grid% ng2sm1 = grid% ng2s - 1
    grid% ng2sm2 = grid% ng2s - 2
    grid% ng2ep1 = grid% ng2e + 1
    grid% ng2ep2 = grid% ng2e + 2

    grid% ngg1s = 0
    grid% ngg1e = grid% ni

    grid% ngg2s = 1
    grid% ngg2e = grid% ni + 1
    !-------------------------------------------
    ! Define number of possible multigrid levels
    !-------------------------------------------
    IF (grid% nir /= 1) THEN
      grid% ni2_mg = grid% ni2
    ELSE
      grid% ni2_mg = grid% ni2 - 1
    END IF
    !----------------
    ! Allocate fields
    !----------------
    ALLOCATE ( grid% xnglob (0           :grid% ni    , 1           :grid% ni+1  , 3, grid% nd) )
    ALLOCATE ( grid% xn     (grid% ng1sm2:grid% ng1ep2, grid% ng2sm2:grid% ng2ep2, 3, grid% nd) )
    ALLOCATE ( grid% rlon   (grid% ng1s  :grid% ng1e  , grid% ng2s  :grid% ng2e  ,    grid% nd) )
    ALLOCATE ( grid% rlat   (grid% ng1s  :grid% ng1e  , grid% ng2s  :grid% ng2e  ,    grid% nd) )
    ALLOCATE ( grid% corio  (grid% ng1s  :grid% ng1e  , grid% ng2s  :grid% ng2e  ,    grid% nd) )
    ALLOCATE ( grid% area   (grid% ng1s  :grid% ng1ep1, grid% ng2sm1:grid% ng2e  , 2    ) )

    ALLOCATE ( grid% erlon  (grid% ng1sm2:grid% ng1ep2, grid% ng2sm2:grid% ng2ep2, 2, grid% nd) )
    ALLOCATE ( grid% erlat  (grid% ng1sm2:grid% ng1ep2, grid% ng2sm2:grid% ng2ep2, 3, grid% nd) )

    ALLOCATE ( grid% eta    (grid% ng1sm1:grid% ng1ep1, grid% ng2sm1:grid% ng2ep1, 7) )
    ALLOCATE ( grid% chi    (grid% ng1sm1:grid% ng1ep1, grid% ng2sm1:grid% ng2ep1, 7) )
    ALLOCATE ( grid% cpsi   (grid% ng1sm1:grid% ng1ep1, grid% ng2sm1:grid% ng2ep1, 6) )
    ALLOCATE ( grid% spsi   (grid% ng1sm1:grid% ng1ep1, grid% ng2sm1:grid% ng2ep1, 6) )

    ALLOCATE ( grid% bary   (grid% ng1s  :grid% ng1e  , grid% ng2s  :grid% ng2e  , 6) )

    ALLOCATE ( grid% rarn   (grid% ng1sm1:grid% ng1ep1, grid% ng2sm1:grid% ng2ep1   ) )
    ALLOCATE ( grid% hexwgt (grid% ng1s  :grid% ng1e  , grid% ng2s  :grid% ng2e     ) )

    ALLOCATE ( grid% grd    (grid% ng1sm1:grid% ng1ep1, grid% ng2sm1:grid% ng2ep1, 7, 2) )
    ALLOCATE ( grid% rlap   (grid% ng1sm1:grid% ng1ep1, grid% ng2sm1:grid% ng2ep1, 7) )
    !---------------------------------------------------
    ! setup MPI, domain-decomposition, and communication
    !---------------------------------------------------
!    CALL setup_parallel
!
!    CALL domain_decomposition (grid% ni, grid% ng1s, grid% ng1e, grid% ng2s, grid% ng2e)
!
!    CALL setup_communication (grid% ngg1s, grid% ngg1e, grid% ngg2s, grid% ngg2e, grid% nd, &
!                              grid% ng1s,  grid% ng1e,  grid% ng2s,  grid% ng2e,  grid% ni)
!
    !-------------------------------------------------------------
    ! 2. Calculate the arrays xn, rlon, rlat, and corio for all 10
    ! diamonds
    !-------------------------------------------------------------
!    CALL global_coordinates (grid% xn, grid% rlon, grid% rlat, &
!      grid% ni, grid% ni2, grid% ni3, grid% nd, &
!      (/grid% ng1s, grid% ng2s, grid% ng3s, grid% ng3s/), &
!      (/grid% ng1e, grid% ng2e, grid% ng3e, grid% ng3e/), ierror)
!
!    IF (ierror /= 0) THEN
!      WRITE (0,*) '  Error in SUBROUTINE global_coordinates!'
!      STOP
!    ENDIF
    !----------------------------------------------------------
    ! 3. Calculate the minimum mesh width and minimum height of
    ! the grid
    !----------------------------------------------------------
    CALL dxdhmin (grid, ierror)

    IF (ierror /= 0) THEN
      WRITE (0,*) '  Error in SUBROUTINE dxdhmin!'
      STOP
    ENDIF
    !-----------------------------------------------------------------
    ! 4. Calculate the arrays erlon, erlat, eta, chi, cpsi, spsi, bary
    ! which describe the local coordinate system
    !-----------------------------------------------------------------
    CALL local_coordinates (grid, ierror)

    IF (ierror /= 0) THEN
      WRITE (0,*) '  Error in SUBROUTINE local_coordinates!'
      STOP
    ENDIF
    !-------------------------------------------------------
    ! 5. Calculate the areas of the triangles for diamond #1
    !-------------------------------------------------------
    CALL triangle_area (grid, ierror)

    IF (ierror /= 0) THEN
      WRITE (0,*) '  Error in SUBROUTINE triangle_area!'
      STOP
    ENDIF
    !-------------------------------------------------------------------
    ! 6. Calculate the areas of the hexagons for extended diamond #1 and
    ! the weights of the gridpoints for diagnostic evaluations
    !-------------------------------------------------------------------
    CALL hexagon_area (grid, dc, ierror)

    IF (ierror /= 0) THEN
      WRITE (0,*) '  Error in SUBROUTINE hexagon_area!'
      STOP
    ENDIF
    !-----------------------------------------------------------------
    ! 7. Calculate the gradient  operator for the extended diamond #1,
    ! calculate the Laplacian operator for the extended diamond #1
    !-----------------------------------------------------------------
    CALL differentiation_operators (grid, ierror)

    IF (ierror /= 0) THEN
      WRITE (0,*) '  Error in SUBROUTINE  differentiation_operators!'
      STOP
    ENDIF

  END SUBROUTINE generate_grid
  !===========================================================================
  SUBROUTINE factorize_nir (kni, kni2, knir, kerror)

    INTEGER, INTENT(in)  :: kni
    INTEGER, INTENT(out) :: kni2    ! number n of factors of 2: 2^n
    INTEGER, INTENT(out) :: knir    ! remaining factor        : 1,3,5,7,9,..
    INTEGER, INTENT(out) :: kerror  ! always 0

    INTEGER  :: mx

    kerror = 0

    mx    = kni
    kni2  = 0
    knir  = 1

    DO WHILE (mx > 1)
      IF (MOD(mx,2) == 0) THEN
        kni2  = kni2+1
        mx    = mx/2
      ELSE
        knir = mx
        mx   = 1
      ENDIF
    ENDDO

  END SUBROUTINE factorize_nir
  !---------------------------------------------------------------------------
  SUBROUTINE global_coordinates (xn, rlon, rlat, &
                                 ni, ni2, nir, nd, lb, ub, kerror)

    USE mo_constants, ONLY: pid5  !,omcor

    integer  ,intent(in)  :: lb(4), ub(4)
    real(wp) ,intent(out) :: xn   (lb(1):,lb(2):,lb(3):,lb(4):)
    real(wp) ,intent(out) :: rlon (lb(1):,lb(2):,lb(3):,lb(4):)
    real(wp) ,intent(out) :: rlat (lb(1):,lb(2):,lb(3):,lb(4):)
    integer  ,intent(in)  :: ni, ni2, nir, nd
    integer  ,intent(out) :: kerror

    REAL(wp) :: zw   ,   & ! the spherical angle in an icosahedron
                           ! subtended by two vertices.
                zcosw,   & ! cosine(zw)
                zsinw,   & ! sine  (zw)
                zsgn ,   & ! zsgn is a hemisphere factor.
                           ! zsgn =  1.0 is north  (diamonds 1- 5)
                           ! zsgn = -1.0 is south  (diamonds 6-10)
                zrlon,   & ! longitude of diamond vertices
                zgamma,  & ! fraction of great circle angle
                zchord,  & ! Cartesian distance between two points
                ztheta,  & ! Great circle angle between two points
                zalpha,  & ! Weighting factor
                zbeta      ! Weighting factor

    INTEGER :: mcosv( nd) ! meridian angle locations of the 10
                               ! non-polar vertices in units of Pi/5
    INTEGER :: ml,           & ! recursive index interval
               ml2,          & ! recursive bisected index interval
               ml3,          & ! trisected index interval
               mi1,          & ! recursive row index of new node
               mi2,          & ! recursive column index of new node
               mm              ! recursive number of subdivisions

    INTEGER :: j1, j2, jd, jb  ! Loop indices

    kerror = 0

    ! 1. Calculate the Cartesian coordinates of the gridpoints of the
    ! icosahedral grid on the unit sphere.  The grid resolution
    ! corresponds to a subdivision of the edges of the original
    ! icosahedral triangles into ni equal parts.

    ! Compute angles associated with the icosahedron.

    zw      = 2*ACOS(1.0_wp/(2*SIN(Pid5)))
    zcosw   = COS(zw)
    zsinw   = SIN(zw)

    ! Compute the local array mcosv, i.e. the meridian angle locations

    DO jd = 1, nd
      IF (MOD(jd,2) == 1) THEN
        mcosv((jd+1)/2) = -1 + (jd - 1) - nd*((jd - 1)/7)
      ELSE
        mcosv(jd/2+5)   = -1 + (jd - 1) - nd*((jd - 1)/7)
      ENDIF
    ENDDO

    ! Loop over the ten diamonds computing diamond vertices (x,y,z)
    ! coordinates and then iteratively bisecting them ni2 times.
    ! First a trisection is performed, if required (ni3=1).

    DO jd = 1, nd     ! Loop over the diamonds

      ! Toggle the hemisphere

      IF (jd >= 6) THEN
        zsgn = -1.0_wp    ! southern
      ELSE
        zsgn =  1.0_wp    ! northern
      ENDIF

      ! Compute the meridian angle for each diamond "home" vertex.

      zrlon = mcosv(jd)*Pid5

      ! Every diamond has one vertex at a pole (N or S).
      ! Label this point (0,1,,) in each diamond, and
      ! initialize it to have the (x,y,z) coordinates of
      ! the pole point on the unit sphere.

      xn (  0,    1, 1, jd) =  0.0_wp
      xn (  0,    1, 2, jd) =  0.0_wp
      xn (  0,    1, 3, jd) =  zsgn

      ! Now initialize the (x,y,z) coordinates of the "home" vertex,
      ! which defines which diamond we are talking about, and label
      ! this point (ni,1,,).

      xn ( ni,    1, 1, jd) =  zsinw*COS(zrlon)
      xn ( ni,    1, 2, jd) =  zsinw*SIN(zrlon)
      xn ( ni,    1, 3, jd) =  zcosw*zsgn

      ! Next initialize the (x,y,z) coordinates for the corner of the
      ! diamond on the same latitude as the (ni,1,,) vertex, which
      ! is (0,ni+1,,)

      xn (  0, ni+1, 1, jd) =  zsinw*COS(zrlon + 2*Pid5)
      xn (  0, ni+1, 2, jd) =  zsinw*SIN(zrlon + 2*Pid5)
      xn (  0, ni+1, 3, jd) =  zcosw*zsgn

      ! Initialize the last diamond vertex, which is located
      ! in the opposite hemisphere as (ni,ni+1,,)

      xn ( ni, ni+1, 1, jd) =  zsinw*COS(zrlon + Pid5)
      xn ( ni, ni+1, 2, jd) =  zsinw*SIN(zrlon + Pid5)
      xn ( ni, ni+1, 3, jd) = -zcosw*zsgn

      ! First a trisection is performed, if required (ni3=1).

      select case (nir)
      case default
        call finish ('global_coordinates','grid root is not 1 or 3')
      case (1)
      case (3)
        !---------------------------------
        ! Trisect the rows of the diamond.
        !---------------------------------
        ml3 =  ni/3
        DO j1 = 1,2
          DO j2 = 1,2
            mi1    = (j1-1)* ni
            mi2    = j2*ml3 + 1
            zgamma = REAL(j2,wp)/3.0_wp

            ! Calculate "zchord", the Cartesian distance between x1 and x2

            zchord =  SQRT ( &
                 (xn (mi1, ni+1,1,jd) - xn (mi1,1,1,jd))**2 + &
                 (xn (mi1, ni+1,2,jd) - xn (mi1,1,2,jd))**2 + &
                 (xn (mi1, ni+1,3,jd) - xn (mi1,1,3,jd))**2 )

            ! Calculate "ztheta", the great circle angle between x1 and x2

            ztheta = 2*ASIN (0.5_wp*zchord)

            ! Calculate the weighting factors which follow from the condition
            ! that x is a point on the unit-sphere, too.

            zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
            zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)

            ! Store the (x,y,z) coordinates of the point x

            xn (mi1,mi2,1,jd) = &
                 zalpha*xn (mi1,1   ,1,jd) + &
                 zbeta *xn (mi1, ni+1,1,jd)
            xn (mi1,mi2,2,jd) = &
                 zalpha*xn (mi1,1   ,2,jd) + &
                 zbeta *xn (mi1, ni+1,2,jd)
            xn (mi1,mi2,3,jd) = &
                 zalpha*xn (mi1,1   ,3,jd) + &
                 zbeta *xn (mi1, ni+1,3,jd)

          ENDDO
        ENDDO

        ! Trisect the columns of the diamond.

        DO j1 = 1,2
          DO j2 = 1,2
            mi1    = j2*ml3
            mi2    = (j1-1)* ni + 1
            zgamma = REAL(j2,wp)/3.0_wp

            ! Calculate "zchord", the Cartesian distance between x1 and x2

            zchord =  SQRT ( &
                  (xn ( ni,mi2,1,jd) &
                 - xn (0,mi2,1,jd))**2  &
                 +(xn ( ni,mi2,2,jd) &
                 - xn (0,mi2,2,jd))**2  &
                 +(xn ( ni,mi2,3,jd) &
                 - xn (0,mi2,3,jd))**2 )

            ! Calculate "ztheta", the great circle angle between x1 and x2

            ztheta = 2*ASIN (0.5_wp*zchord)

            ! Calculate the weighting factors which follow from the condition
            ! that x is a point on the unit-sphere, too.

            zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
            zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)

            ! Store the (x,y,z) coordinates of the point x

            xn (mi1,mi2,1,jd) = &
                 zalpha*xn (0  ,mi2,1,jd) + &
                 zbeta *xn ( ni,mi2,1,jd)
            xn (mi1,mi2,2,jd) = &
                 zalpha*xn (0  ,mi2,2,jd) + &
                 zbeta *xn ( ni,mi2,2,jd)
            xn (mi1,mi2,3,jd) = &
                 zalpha*xn (0  ,mi2,3,jd) + &
                 zbeta *xn ( ni,mi2,3,jd)

          ENDDO
        ENDDO

        ! Trisect the diagonal of the diamond.

        DO j2 = 1,2
          mi1 =  ni - j2*ml3
          mi2 =   1 + j2*ml3
          zgamma = REAL(j2,wp)/3.0_wp

          ! Calculate "zchord", the Cartesian distance between x1 and x2

          zchord = SQRT ( &
                (xn (0, ni+1,1,jd) &
               - xn ( ni,1,1,jd))**2 &
               +(xn (0, ni+1,2,jd) &
               - xn ( ni,1,2,jd))**2 &
               +(xn (0, ni+1,3,jd) &
               - xn ( ni,1,3,jd))**2 )

          ! Calculate "ztheta", the great circle angle between x1 and x2

          ztheta = 2*ASIN (0.5_wp*zchord)

          ! Calculate the weighting factors which follow from the condition
          ! that x is a point on the unit-sphere, too.

          zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
          zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)

          ! Store the (x,y,z) coordinates of the point x

          xn (mi1,mi2,1,jd) = &
               zalpha*xn ( ni,1    ,1,jd) +  &
               zbeta *xn (0  , ni+1,1,jd)
          xn (mi1,mi2,2,jd) = &
               zalpha*xn ( ni,1    ,2,jd) +  &
               zbeta *xn (0  , ni+1,2,jd)
          xn (mi1,mi2,3,jd) = &
               zalpha*xn ( ni,1    ,3,jd) +  &
               zbeta *xn (0  , ni+1,3,jd)

        ENDDO

        ! Compute coordinates of icosahedral triangle centers.

        CALL triangle_center (xn, ni, jd)

      END SELECT     ! End of trisection

      ! Find the coordinates of the triangle nodes by iteratively
      ! bisecting the diamond intervals.

      DO jb = 0,  ni2-1
        mm  = (nir)*(2**jb)
        ml  =  ni/mm
        ml2 = ml/2

        ! Compute the rows of the diamond.

        DO j1 = 1,mm+1
          DO j2 = 1,mm
            mi1 = (j1-1)*ml
            mi2 = (j2-1)*ml + ml2 + 1
            zgamma = 0.5_wp

            ! Calculate "zchord", the Cartesian distance between x1 and x2

            zchord = SQRT (  &
                 (xn (mi1,mi2+ml2,1,jd) - xn (mi1,mi2-ml2,1,jd))**2 +  &
                 (xn (mi1,mi2+ml2,2,jd) - xn (mi1,mi2-ml2,2,jd))**2 +  &
                 (xn (mi1,mi2+ml2,3,jd) - xn (mi1,mi2-ml2,3,jd))**2 )

            ! Calculate "ztheta", the great circle angle between x1 and x2

            ztheta = 2*ASIN (0.5_wp*zchord)

            ! Calculate the weighting factors which follow from the condition
            ! that x is a point on the unit-sphere, too.

            zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
            zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)

            ! Store the (x,y,z) coordinates of the point x

            xn (mi1,mi2,1,jd) = zalpha*xn (mi1,mi2-ml2,1,jd) +  &
                 zbeta *xn (mi1,mi2+ml2,1,jd)
            xn (mi1,mi2,2,jd) = zalpha*xn (mi1,mi2-ml2,2,jd) +  &
                 zbeta *xn (mi1,mi2+ml2,2,jd)
            xn (mi1,mi2,3,jd) = zalpha*xn (mi1,mi2-ml2,3,jd) +  &
                 zbeta *xn (mi1,mi2+ml2,3,jd)

          ENDDO
        ENDDO

        ! Compute the columns of diamond.

        DO j1 = 1,mm+1
          DO j2 = 1,mm
            mi1 = (j2-1)*ml + ml2
            mi2 = (j1-1)*ml + 1
            zgamma = 0.5_wp

            ! Calculate "zchord", the Cartesian distance between x1 and x2

            zchord = SQRT (  &
                 (xn (mi1+ml2,mi2,1,jd) - xn (mi1-ml2,mi2,1,jd))**2 +  &
                 (xn (mi1+ml2,mi2,2,jd) - xn (mi1-ml2,mi2,2,jd))**2 +  &
                 (xn (mi1+ml2,mi2,3,jd) - xn (mi1-ml2,mi2,3,jd))**2 )

            ! Calculate "ztheta", the great circle angle between x1 and x2

            ztheta = 2*ASIN (0.5_wp*zchord)

            ! Calculate the weighting factors which follow from the condition
            ! that x is a point on the unit-sphere, too.

            zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
            zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)

            ! Store the (x,y,z) coordinates of the point x

            xn (mi1,mi2,1,jd) = zalpha*xn (mi1-ml2,mi2,1,jd) +  &
                 zbeta *xn (mi1+ml2,mi2,1,jd)
            xn (mi1,mi2,2,jd) = zalpha*xn (mi1-ml2,mi2,2,jd) +  &
                 zbeta *xn (mi1+ml2,mi2,2,jd)
            xn (mi1,mi2,3,jd) = zalpha*xn (mi1-ml2,mi2,3,jd) +  &
                 zbeta *xn (mi1+ml2,mi2,3,jd)
          ENDDO
        ENDDO

        ! Compute the diagonals of the diamond.

        DO j1 = 1,mm
          DO j2 = 1,mm
            mi1 = (j1-1)*ml + ml2
            mi2 = (j2-1)*ml + ml2 + 1
            zgamma = 0.5_wp

            ! Calculate "zchord", the Cartesian distance between x1 and x2

            zchord = SQRT (  &
                 (xn (mi1+ml2,mi2-ml2,1,jd) -       &
                 xn (mi1-ml2,mi2+ml2,1,jd))**2 +  &
                 (xn (mi1+ml2,mi2-ml2,2,jd) -      &
                 xn (mi1-ml2,mi2+ml2,2,jd))**2 +  &
                 (xn (mi1+ml2,mi2-ml2,3,jd) -      &
                 xn (mi1-ml2,mi2+ml2,3,jd))**2 )

            ! Calculate "ztheta", the great circle angle between x1 and x2

            ztheta = 2*ASIN (0.5_wp*zchord)

            ! Calculate the weighting factors which follow from the condition
            ! that x is a point on the unit-sphere, too.

            zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
            zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)

            ! Store the (x,y,z) coordinates of the point x

            xn (mi1,mi2,1,jd) =  &
                 zalpha*xn (mi1-ml2,mi2+ml2,1,jd) +  &
                 zbeta *xn (mi1+ml2,mi2-ml2,1,jd)
            xn (mi1,mi2,2,jd) =  &
                 zalpha*xn (mi1-ml2,mi2+ml2,2,jd) +  &
                 zbeta *xn (mi1+ml2,mi2-ml2,2,jd)
            xn (mi1,mi2,3,jd) =  &
                 zalpha*xn (mi1-ml2,mi2+ml2,3,jd) +  &
                 zbeta *xn (mi1+ml2,mi2-ml2,3,jd)

          ENDDO
        ENDDO

      ENDDO       ! end loop over bisections

      ! Set xnglob to 0 if it is less than 2.5 e-14 to avoid round-off
      ! errors


      WHERE (ABS(xn (:,:,:,jd)) < 2.5e-14_wp)
        xn (:,:,:,jd) = 0.0_wp
      END WHERE

#ifdef DEBUG
      field_name = 'zxn (*)'
      DO j1 = 1,3
        WRITE (field_name(5:5),'(i1)') j1
        CALL print_array ( &
             xn (lb(1):ub(1),lb(2):ub(2),j1,jd), &
             lb(1), ub(1), lb(2), ub(2), &
             lb(1), ub(1), lb(2), ub(2), jd, 1, 0, &
             field_name(1:6), '-', 'glo_coor', 'h', &
             1000._wp, 0._wp,  kerror)
      ENDDO
#endif

    ENDDO         ! end loop over diamonds

    ! Extend the nodal array by two rows/colums around
    ! This is done by the subroutine *xd_p* for all diamonds and the
    ! three Cartesian (x,y,z) coordinates simultaneously.

!    grid% xn (lb(1):ub(1),lb(2):ub(2),:,:) = &
!         xn (lb(1):ub(1),lb(2):ub(2),:,:)
!
!    CALL set_boundaries ( xn, &
!         lb(1)-2, ub(1)+2, lb(2)-2, ub(2)+2, 1, 3, nd, 2 )

#ifdef DEBUG
    DO jd = 1, nd     ! Loop over the diamonds
      field_name = 'pxn (*)'
      DO j1 = 1,3
        WRITE (field_name(5:5),'(i1)') j1
        CALL print_array ( &
              xn (lb(1)-2:ub(1)+2,lb(2)-2:ub(2)+2,j1,jd), &
             lb(1)-2, ub(1)+2, lb(2)-2, ub(2)+2, &
             lb(1)-2, ub(1)+2, lb(2)-2, ub(2)+2, jd, 1, 0, &
             field_name(1:6), '-', 'glo_coor', 'h', &
             1000._wp, 0._wp,  kerror)
      ENDDO
    ENDDO         ! end loop over diamonds
#endif

    ! Calculate the longitude "prlon" and the latitude "rlat";
    ! only for the core of the diamonds, not the extended ones.

!$omp parallel do private(j1,j2,jd) schedule(static)
    DO jd = 1, nd     ! Loop over the diamonds

      DO j2 = lb(2), ub(2)
        DO j1 = lb(1), ub(1)
          rlon (j1,j2,1,jd) = ATAN2(xn (j1,j2,2,jd), xn (j1,j2,1,jd) + 1.e-20_wp)
          rlat (j1,j2,1,jd) = ASIN (xn (j1,j2,3,jd))
!         corio(j1,j2,jd) = 2*Omcor*xn (j1,j2,3,jd)
        ENDDO
      ENDDO

    ENDDO         ! end loop over diamonds
!$omp end parallel do

#ifdef DEBUG
    DO jd = 1, nd     ! Loop over the diamonds
      CALL print_array ( &
           rlon (lb(1):ub(1),lb(2):ub(2),1,jd), &
           lb(1)  , ub(1)  , lb(2)  , ub(2)  , &
           lb(1)  , ub(1)  , lb(2)  , ub(2)  , jd, 1, 0, &
           'rlon  ', 'Degree','glo_coor', 'h', &
           572.96_wp , 0._wp,  kerror)
      CALL print_array ( &
           rlat (lb(1):ub(1),lb(2):ub(2),1,jd), &
           lb(1)  , ub(1)  , lb(2)  , ub(2)  , &
           lb(1)  , ub(1)  , lb(2)  , ub(2)  , jd, 1, 0, &
           'rlat  ', 'Degree', 'glo_coor', 'h', &
           572.96_wp , 0._wp,  kerror)
!     CALL print_array ( &
!          corio(lb(1):ub(1),lb(2):ub(2),jd), &
!          lb(1)  , ub(1)  , lb(2)  , ub(2)  , &
!          lb(1)  , ub(1)  , lb(2)  , ub(2)  , jd, 1, 0, &
!          'corio ', '1/s'  , 'glo_coor', 'h', &
!          1.e5_wp   , 0._wp,  kerror)
    ENDDO         ! end loop over diamonds
#endif

  END SUBROUTINE global_coordinates
  !---------------------------------------------------------------------------
  PURE SUBROUTINE triangle_center (xn, ni, kjd)

    real(wp), INTENT(inout) :: xn (0:,:,:,:)
    INTEGER,  INTENT(in)    :: ni
    INTEGER,  INTENT(in)    :: kjd

    REAL(wp) :: zxnorm
    REAL(wp) :: zlon, zlat

    INTEGER :: j,           & ! loop index
               mi1,         & ! index of center point
               mi2            ! index of top or bottom diamond corner

    DO j = 1, 2 ! Loop over the two triangles

      mi1  = ni/3 * j
      mi2  = 1 + (j - 1)*ni

      xn (mi1,mi1+1,1,kjd) = xn (mi2-1   ,mi2      ,1,kjd) +   &
                                      xn (ni,1     ,1,kjd) +   &
                                      xn (0   ,ni+1,1,kjd)
      xn (mi1,mi1+1,2,kjd) = xn (mi2-1   ,mi2      ,2,kjd) +   &
                                      xn (ni,1     ,2,kjd) +   &
                                      xn (0   ,ni+1,2,kjd)
      xn (mi1,mi1+1,3,kjd) = xn (mi2-1   ,mi2      ,3,kjd) +   &
                                      xn (ni,1     ,3,kjd) +   &
                                      xn (0   ,ni+1,3,kjd)

      ! Normalize to unit-sphere

      zxnorm = 1./SQRT (xn (mi1,mi1+1,1,kjd)**2 +  &
                        xn (mi1,mi1+1,2,kjd)**2 +  &
                        xn (mi1,mi1+1,3,kjd)**2)

      xn (mi1,mi1+1,1,kjd) = zxnorm*xn (mi1,mi1+1,1,kjd)
      xn (mi1,mi1+1,2,kjd) = zxnorm*xn (mi1,mi1+1,2,kjd)
      xn (mi1,mi1+1,3,kjd) = zxnorm*xn (mi1,mi1+1,3,kjd)

      call xyz2ll (xn (mi1,mi1+1,1,kjd), xn (mi1,mi1+1,2,kjd),&
                   xn (mi1,mi1+1,3,kjd), zlon, zlat)

    ENDDO    ! End of loop over triangles

  END SUBROUTINE triangle_center
  !-----------------------------------------------------------------------------
  SUBROUTINE set_boundaries (p, dc, kip1s, kip1e, kip2s, kip2e, kip3s, kip3e, knd, kr)

    ! Comment on kip3s and kip3e: This dimension is used either for levels
    ! or sometimes as counter of independent coordinates (like in xn)

    type (t_atm_dec) ,intent(in) :: dc

    INTEGER, INTENT(in)     :: kip1s, kip1e, kip2s, kip2e, kip3s, kip3e, knd, kr

    REAL(wp), INTENT(inout) :: p(kip1s:kip1e, kip2s:kip2e, kip3s:kip3e, knd)

    REAL(wp) :: zbuff_send(kip3s:kip3e,dc% np_send_max)
    REAL(wp) :: zbuff_recv(kip3s:kip3e,dc% np_recv_max)

    INTEGER :: j1, j2, j3, jd, jp, jp1, j
    INTEGER :: mpr, mpl, mreq

#ifndef NOMPI
    INTEGER :: mpi_req                       (dace% npe)
    INTEGER :: mpi_statuses (MPI_STATUS_SIZE, dace% npe)
#endif
    mpl = kip3e - kip3s + 1 ! number of layers

    ! Put all points we have to send into the send buffer

    DO j3 = kip3s, kip3e
      DO jp = 0, dc% nproc1*dc% nproc2-1

        IF (kr == 1) mpr = dc% np_send_1(jp)
        IF (kr == 2) mpr = dc% np_send_2(jp)

!dir$ ivdep
#if 0
!dir$ cache_bypass zbuff_send, p
#endif

        DO j = dc% np_send_s(jp), dc% np_send_s(jp) + mpr - 1

          j1 = dc% idx_send(1,j)
          j2 = dc% idx_send(2,j)
          jd = dc% idx_send(3,j)

          zbuff_send(j3,j) = p(j1,j2,j3,jd)

        ENDDO
      ENDDO
    ENDDO

    ! Setup non blocking receives on all processors from which
    ! we expect data

    mreq = 0

    DO jp = 0, dc% nproc1*dc% nproc2-1

      IF (jp == dc% myproc) CYCLE ! Don't receive from ourself

      IF (kr == 1) mpr =  dc% np_recv_1(jp)
      IF (kr == 2) mpr =  dc% np_recv_2(jp)

      IF (mpr > 0) THEN
        mreq = mreq+1
#ifndef NOMPI
        CALL MPI_Irecv( zbuff_recv(kip3s,dc% np_recv_s(jp)), mpr*mpl,  &
             p_real, jp, 1, dc% comm,           &
             mpi_req(mreq), mpi_err)
#endif
      ENDIF

    ENDDO

    ! Send data to all processors who want to have data from us

    DO jp1 = 0, dc% nproc1*dc% nproc2-1

      jp = MOD( jp1+dc% myproc, dc% nproc1*dc% nproc2 )

      IF (kr == 1) mpr =  dc% np_send_1(jp)
      IF (kr == 2) mpr =  dc% np_send_2(jp)

      IF (jp == dc% myproc) THEN

        ! Don't send to ourself, just copy

!dir$ ivdep
#if 0
!dir$ cache_bypass zbuff_send, zbuff_recv
#endif

        DO j3 = kip3s, kip3e
          DO j = 1, mpr
            zbuff_recv(j3,dc% np_recv_s(jp)+j-1) =       &
                 zbuff_send(j3,dc% np_send_s(jp)+j-1)
          ENDDO
        ENDDO

      ELSE

        IF (mpr > 0) THEN
#ifndef NOMPI
          CALL MPI_Send( zbuff_send(kip3s,dc% np_send_s(jp)), mpr*mpl,  &
               p_real, jp, 1, dc% comm,          &
               mpi_err)
#endif
        ENDIF

      ENDIF

    ENDDO

    ! Wait until all receives have finished

    IF (mreq > 0) THEN
#ifndef NOMPI
      CALL MPI_Waitall(mreq, mpi_req, mpi_statuses, mpi_err)
#endif
    ENDIF

    ! Fill received data into it's location

    DO j3 = kip3s, kip3e
      DO jp = 0, dc% nproc1*dc% nproc2-1

        IF (kr == 1) mpr =  dc% np_recv_1(jp)
        IF (kr == 2) mpr =  dc% np_recv_2(jp)

!dir$ ivdep
#if 0
!dir$ cache_bypass zbuff_recv, p
#endif

        DO j = dc% np_recv_s(jp), dc% np_recv_s(jp) + mpr - 1

          j1 = dc% idx_recv(1,j)
          j2 = dc% idx_recv(2,j)
          jd = dc% idx_recv(3,j)

          p(j1,j2,j3,jd) = zbuff_recv(j3,j)

        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE set_boundaries
  !---------------------------------------------------------------------------
  SUBROUTINE dxdhmin (grid, kerror)

    USE mo_constants, ONLY: Pid5, Re
#ifndef NOMPI
    USE mo_atm_decomp,  ONLY: lmpi
#endif

    TYPE (triangular_grid), INTENT(inout) ,target :: grid
    INTEGER,                INTENT(out)           :: kerror


    ! scalar product of the location vectors of the
    ! home node and neighbouring node # 1 - 6

    REAL(wp) :: zsp1, zsp2, zsp3, zsp4, zsp5, zsp6

    ! scalar product of the location vectors of
    ! neigbouring nodes 1 and 2, 2 and 3, ...

    REAL(wp) :: zsp12, zsp23, zsp34, zsp45, zsp56, zsp61, zlon, zlat

    REAL(wp) :: zdxmax, zdxmin, & ! maximum/minimum mesh width of grid
                zdhmax, zdhmin    ! maximum/minimum height of triangle

#ifndef NOMPI
    REAL(wp) :: zdxdhs(2), zdxdhr(2) ! for MPI_Allreduce
#endif

    INTEGER  :: j1, j2 ! loop indices

    type (t_atm_dec) ,pointer :: dc

    dc => grid% dc

    zdxmin = 1.e20_wp
    zdhmin = 1.e20_wp
    zdxmax = 0.0_wp
    zdhmax = 0.0_wp

    kerror = 0

    ! The scalar product between the nodal vectors pxn of the home
    ! node and the 6 (5) surrounding nodes is equal to the cosine of
    ! the angle epsilon between the nodes and the home node. Epsilon
    ! is a measure of the great circle arc distance of the gridpoints.
    ! The calculation is performed only for diamond 1 since the relative
    ! distances are the same in all diamonds.

    DO j2 = grid% ng2s, grid% ng2e
      DO j1  = grid% ng1s, grid% ng1e
        zsp1   = grid% xn (j1,j2,1,1)*grid% xn (j1+nsp11,j2+nsp21,1,1) +  &
                 grid% xn (j1,j2,2,1)*grid% xn (j1+nsp11,j2+nsp21,2,1) +  &
                 grid% xn (j1,j2,3,1)*grid% xn (j1+nsp11,j2+nsp21,3,1)

        call xyz2ll(grid% xn (j1+nsp11,j2+nsp21,1,1),grid% xn (j1+nsp11,j2+nsp21,2,1), &
                    grid% xn (j1+nsp11,j2+nsp21,3,1), zlon, zlat)

        zsp1   = ACOS (zsp1)
        zdxmin = MIN (zdxmin, zsp1)
        zdxmax = MAX (zdxmax, zsp1)
        zsp2   = grid% xn (j1,j2,1,1)*grid% xn (j1+nsp12,j2+nsp22,1,1) +  &
                 grid% xn (j1,j2,2,1)*grid% xn (j1+nsp12,j2+nsp22,2,1) +  &
                 grid% xn (j1,j2,3,1)*grid% xn (j1+nsp12,j2+nsp22,3,1)

        call xyz2ll(grid% xn (j1+nsp12,j2+nsp22,1,1),grid% xn (j1+nsp12,j2+nsp22,2,1), &
                    grid% xn (j1+nsp12,j2+nsp22,3,1), zlon, zlat)

        zsp2   = ACOS (zsp2)
        zdxmin = MIN (zdxmin, zsp2)
        zdxmax = MAX (zdxmax, zsp2)
        zsp3   = grid% xn (j1,j2,1,1)*grid% xn (j1+nsp13,j2+nsp23,1,1) +  &
                 grid% xn (j1,j2,2,1)*grid% xn (j1+nsp13,j2+nsp23,2,1) +  &
                 grid% xn (j1,j2,3,1)*grid% xn (j1+nsp13,j2+nsp23,3,1)

        call xyz2ll(grid% xn (j1+nsp13,j2+nsp23,1,1),grid% xn (j1+nsp13,j2+nsp23,2,1), &
             grid% xn (j1+nsp13,j2+nsp23,3,1), zlon, zlat)

        zsp3   = ACOS (zsp3)
        zdxmin = MIN (zdxmin, zsp3)
        zdxmax = MAX (zdxmax, zsp3)
        zsp4   = grid% xn (j1,j2,1,1)*grid% xn (j1+nsp14,j2+nsp24,1,1) +  &
                 grid% xn (j1,j2,2,1)*grid% xn (j1+nsp14,j2+nsp24,2,1) +  &
                 grid% xn (j1,j2,3,1)*grid% xn (j1+nsp14,j2+nsp24,3,1)

        call xyz2ll(grid% xn (j1+nsp14,j2+nsp24,1,1),grid% xn (j1+nsp14,j2+nsp24,2,1), &
                    grid% xn (j1+nsp14,j2+nsp24,3,1), zlon, zlat)

        zsp4   = ACOS (zsp4)
        zdxmin = MIN (zdxmin, zsp4)
        zdxmax = MAX (zdxmax, zsp4)
        zsp5   = grid% xn (j1,j2,1,1)*grid% xn (j1+nsp15,j2+nsp25,1,1) +  &
                 grid% xn (j1,j2,2,1)*grid% xn (j1+nsp15,j2+nsp25,2,1) +  &
                 grid% xn (j1,j2,3,1)*grid% xn (j1+nsp15,j2+nsp25,3,1)

        call xyz2ll(grid% xn (j1+nsp15,j2+nsp25,1,1),grid% xn (j1+nsp15,j2+nsp25,2,1), &
                    grid% xn (j1+nsp15,j2+nsp25,3,1), zlon, zlat)

        zsp5   = ACOS (zsp5)
        zdxmin = MIN (zdxmin, zsp5)
        zdxmax = MAX (zdxmax, zsp5)
        zsp6   = grid% xn (j1,j2,1,1)*grid% xn (j1+nsp16,j2+nsp26,1,1) +  &
                 grid% xn (j1,j2,2,1)*grid% xn (j1+nsp16,j2+nsp26,2,1) +  &
                 grid% xn (j1,j2,3,1)*grid% xn (j1+nsp16,j2+nsp26,3,1)

        call xyz2ll(grid% xn (j1+nsp16,j2+nsp26,1,1),grid% xn (j1+nsp16,j2+nsp26,2,1), &
                    grid% xn (j1+nsp16,j2+nsp26,3,1), zlon, zlat)

        zsp6   = ACOS (zsp6)
        zdxmin = MIN (zdxmin, zsp6)
        zdxmax = MAX (zdxmax, zsp6)
      ENDDO
    ENDDO

    ! Now the same procedure for the nodes 1 and 2, 2 and 3 and so on.
    ! To avoid problems at the 4 special points (corners of the
    ! diamond with only 5 neighbours) the calculation omits the
    ! outer row/column.

    DO j2 = grid% ng2s+1, grid% ng2e-1
      DO j1  = grid% ng1s+1, grid% ng1e-1
        zsp12  = grid% xn (j1+nsp11,j2+nsp21,1,1)*grid% xn (j1+nsp12,j2+nsp22,1,1) +  &
                 grid% xn (j1+nsp11,j2+nsp21,2,1)*grid% xn (j1+nsp12,j2+nsp22,2,1) +  &
                 grid% xn (j1+nsp11,j2+nsp21,3,1)*grid% xn (j1+nsp12,j2+nsp22,3,1)
        zsp12  = ACOS (zsp12)
        zdxmin = MIN (zdxmin, zsp12)
        zdxmax = MAX (zdxmax, zsp12)
        zsp23  = grid% xn (j1+nsp12,j2+nsp22,1,1)*grid% xn (j1+nsp13,j2+nsp23,1,1) +  &
                 grid% xn (j1+nsp12,j2+nsp22,2,1)*grid% xn (j1+nsp13,j2+nsp23,2,1) +  &
                 grid% xn (j1+nsp12,j2+nsp22,3,1)*grid% xn (j1+nsp13,j2+nsp23,3,1)
        zsp23  = ACOS (zsp23)
        zdxmin = MIN (zdxmin, zsp23)
        zdxmax = MAX (zdxmax, zsp23)
        zsp34  = grid% xn (j1+nsp13,j2+nsp23,1,1)*grid% xn (j1+nsp14,j2+nsp24,1,1) +  &
                 grid% xn (j1+nsp13,j2+nsp23,2,1)*grid% xn (j1+nsp14,j2+nsp24,2,1) +  &
                 grid% xn (j1+nsp13,j2+nsp23,3,1)*grid% xn (j1+nsp14,j2+nsp24,3,1)
        zsp34  = ACOS (zsp34)
        zdxmin = MIN (zdxmin, zsp34)
        zdxmax = MAX (zdxmax, zsp34)
        zsp45  = grid% xn (j1+nsp14,j2+nsp24,1,1)*grid% xn (j1+nsp15,j2+nsp25,1,1) +  &
                 grid% xn (j1+nsp14,j2+nsp24,2,1)*grid% xn (j1+nsp15,j2+nsp25,2,1) +  &
                 grid% xn (j1+nsp14,j2+nsp24,3,1)*grid% xn (j1+nsp15,j2+nsp25,3,1)
        zsp45  = ACOS (zsp45)
        zdxmin = MIN (zdxmin, zsp45)
        zdxmax = MAX (zdxmax, zsp45)
        zsp56  = grid% xn (j1+nsp15,j2+nsp25,1,1)*grid% xn (j1+nsp16,j2+nsp26,1,1) +  &
                 grid% xn (j1+nsp15,j2+nsp25,2,1)*grid% xn (j1+nsp16,j2+nsp26,2,1) +  &
                 grid% xn (j1+nsp15,j2+nsp25,3,1)*grid% xn (j1+nsp16,j2+nsp26,3,1)
        zsp56  = ACOS (zsp56)
        zdxmin = MIN (zdxmin, zsp56)
        zdxmax = MAX (zdxmax, zsp56)
        zsp61  = grid% xn (j1+nsp16,j2+nsp26,1,1)*grid% xn (j1+nsp11,j2+nsp21,1,1) +  &
                 grid% xn (j1+nsp16,j2+nsp26,2,1)*grid% xn (j1+nsp11,j2+nsp21,2,1) +  &
                 grid% xn (j1+nsp16,j2+nsp26,3,1)*grid% xn (j1+nsp11,j2+nsp21,3,1)
        zsp61  = ACOS (zsp61)
        zdxmin = MIN (zdxmin, zsp61)
        zdxmax = MAX (zdxmax, zsp61)
      ENDDO
    ENDDO

    ! Communication for MPI

#ifndef NOMPI
    IF (lmpi) THEN
      zdxdhs(1) = zdxmax
      zdxdhs(2) = zdhmax
!      CALL MPI_Allreduce(zdxdhs,zdxdhr,2,MPI_REAL,MPI_MAX,  &
!           dc% comm ,mpi_err)
      zdxdhr = p_max (zdxdhs, dc% comm)
      zdxmax = zdxdhr(1)
      zdhmax = zdxdhr(2)
      zdxdhs(1) = zdxmin
      zdxdhs(2) = zdhmin
!      CALL MPI_Allreduce(zdxdhs,zdxdhr,2,MPI_REAL,MPI_MIN,  &
!           dc% comm ,mpi_err)
      zdxdhr = p_min (zdxdhs, dc% comm)
      zdxmin = zdxdhr(1)
      zdhmin = zdxdhr(2)
    ENDIF
#endif

    ! Calculation of the minimum height of the triangles; as a good
    ! estimate, we consider a spherical triangle with three sides of
    ! equal length "zdxmin"; the same is done for "zdhmax".

    zdhmin = ASIN (SIN (zdxmin)*SIN (2*Pid5))
    zdhmax = ASIN (SIN (zdxmax)*SIN (2*Pid5))

    grid% dhmin = zdhmin
    grid% dhmax = zdhmax
    grid% dxmin = zdxmin
    grid% dxmax = zdxmax

    WRITE (6,'(a,f12.3,a)') &
         'dxdhmin: minimum mesh width          = ', zdxmin*Re*1e-3, ' km'
    WRITE (6,'(a,f12.3,a)') &
         'dxdhmin: maximum mesh width          = ', zdxmax*Re*1e-3, ' km'

    WRITE (6,'(a,f12.3,a)') &
         'dxdhmin: minimum triangle height     = ', zdhmin*Re*1e-3, ' km'
    WRITE (6,'(a,f12.3,a)') &
         'dxdhmin: maximum triangle height     = ', zdhmax*Re*1e-3, ' km'

  END SUBROUTINE dxdhmin
  !---------------------------------------------------------------------------
  SUBROUTINE setmrp (grid)

    TYPE (triangular_grid), INTENT(inout) :: grid

    grid% ni1mrp(:) = HUGE(0)
    grid% ni2mrp(:) = HUGE(0)

    ! Define the four mirrored points

    IF (grid% ng1s == grid% ngg1s .AND. grid% ng2s == grid% ngg2s) THEN
      grid% ni1mrp (1) = grid% ng1s
      grid% ni2mrp (1) = grid% ng2s - 1
      grid% ni1mrp (5) = grid% ng1s - 1
      grid% ni2mrp (5) = grid% ng2s
    ENDIF

    IF (grid% ng1e == grid% ngg1e .AND. grid% ng2e == grid% ngg2e) THEN
      grid% ni1mrp (2) = grid% ng1e
      grid% ni2mrp (2) = grid% ng2e + 1
      grid% ni1mrp (6) = grid% ng1e + 1
      grid% ni2mrp (6) = grid% ng2e
    ENDIF

    IF (grid% ng1e == grid% ngg1e .AND. grid% ng2s == grid% ngg2s) THEN
      grid% ni1mrp (3) = grid% ng1e + 1
      grid% ni2mrp (3) = grid% ng2s - 1
      grid% ni1mrp (7) = grid% ng1e
      grid% ni2mrp (7) = grid% ng2s - 1
    ENDIF

    IF (grid% ng1s == grid% ngg1s .AND. grid% ng2e == grid% ngg2e) THEN
      grid% ni1mrp (4) = grid% ng1s - 1
      grid% ni2mrp (4) = grid% ng2e
      grid% ni1mrp (8) = grid% ng1s - 1
      grid% ni2mrp (8) = grid% ng2e + 1
    ENDIF

  END SUBROUTINE setmrp
  !---------------------------------------------------------------------------
  SUBROUTINE local_coordinates (grid, kerror)

    TYPE (triangular_grid), INTENT(inout) :: grid
    INTEGER,                INTENT(out)   :: kerror

    REAL(wp) :: zd12, & ! scalar product of two vectors
                zd13, & ! scalar product of two vectors
                zd22, & ! scalar product of two vectors
                zd23, & ! scalar product of two vectors
                z1  , & ! local temporary
                z2  , & ! local temporary
                zn  , & ! local norm
                zsn , & ! hemisphere discriminator
                zchi, & ! pi/5
                zeta, & ! meridian angle of diamond
                zd1 , & ! scalar product of two vectors
                zd2 , & ! scalar product of two vectors
                zd3     ! scalar product of two vectors

    INTEGER :: j1, j2, jd, jk ,jm, js  ! Loop indices
    INTEGER :: js1, js2                ! Indices of a gridpoint

    kerror = 0

    ! 1. Use later "nspoke" and "nspokes", the offsets of the 6 (5) neigh-
    !    bouring gridpoints relative to the central node for 2-d array
    !    addressing. "nspokes" will be used at the four mirrored points
    !    of arrays which have been extended by one row/column since
    !    there the normal spokes "nspoke" will take wrong values due to
    !    the pentagonal structure at the corners of the diamonds.
    !    Define the i1-, i2-indices of the four mirrored points of the
    !    extended array. - Definition now as parameter fields in this module.

    CALL setmrp (grid)

    ! 2. Compute the Cartesian coordinates of the local unit vectors
    !    "rlon" (aligned to the global east direction) and "rlat"
    !    (aligned to the global north direction)

    DO jd = 1, grid% nd                 ! Loop over all diamonds

      DO j2 = grid% ng2sm2, grid% ng2ep2     ! Loop over the extended array
        DO j1 = grid% ng1sm2, grid% ng1ep2   ! Loop over the extended array
          z1    = -grid% xn (j1,j2,2,jd)   ! y-coordinate of node vector
          z2    =  grid% xn (j1,j2,1,jd)   ! x-coordinate of node vector

!LK Comment 1.e-30_wp should be epsilon(wp) ?

          zn    = 1.0_wp/SQRT (z1**2 + z2**2 + 1.e-30_wp)

          grid% erlon(j1,j2,1,jd) = z1*zn ! Normalize to unit length
          grid% erlon(j1,j2,2,jd) = z2*zn ! Normalize to unit length

          ! z-coordinate of erlon is 0
          ! erlat is the cross product of grid% xn and erlon

          grid% erlat(j1,j2,1,jd) = -grid% xn (j1,j2,3,jd)*z2*zn
          grid% erlat(j1,j2,2,jd) =  grid% xn (j1,j2,3,jd)*z1*zn
          grid% erlat(j1,j2,3,jd) =  grid% xn (j1,j2,1,jd)*z2*zn &
                                    -grid% xn (j1,j2,2,jd)*z1*zn
        ENDDO
      ENDDO

      ! Define the local unit vectors for the north and south pole
      ! separately (zn is undefined there).
      ! Attention:  The local unit vectors erlon, erlat for the
      !             poles differ from one diamond to the next!

      zchi = 0.8_wp*ASIN(1.0_wp)       ! 2*pi/5 BETTER USE 2*Pid5
      zsn  = 1.0_wp                    ! hemisphere descriminator
      IF (jd >= 6) zsn = -1.0_wp
      zeta = (jd-1.5_wp)*zchi
      IF (jd >= 6) zeta = (jd-6)*zchi

      IF (grid% ng1s == 0 .AND. grid% ng2s == 1) THEN   ! Pole points
        grid% erlon(0,1,1,jd) = - SIN (zeta)
        grid% erlon(0,1,2,jd) =   COS (zeta)
        grid% erlat(0,1,1,jd) = - COS (zeta)*zsn
        grid% erlat(0,1,2,jd) = - SIN (zeta)*zsn
        grid% erlat(0,1,3,jd) = 0.0_wp
      ENDIF

#ifdef DEBUG
      field_name = 'erlon(*)'
      DO j1 = 1,2
        WRITE (field_name(7:7),'(i1)') j1
        CALL print_array ( &
             grid% erlon(grid% ng1sm2:grid% ng1ep2,grid% ng2sm2:grid% ng2ep2,j1,jd), &
             grid% ng1sm2, grid% ng1ep2, grid% ng2sm2, grid% ng2ep2, &
             grid% ng1sm2, grid% ng1ep2, grid% ng2sm2, grid% ng2ep2, jd, 1, 0, &
             field_name(1:9), '-', 'loc_coor', 'h', &
             1000._wp, 0._wp,  kerror)
      ENDDO
      field_name = 'erlat(*)'
      DO j1 = 1,3
        WRITE (field_name(7:7),'(i1)') j1
        CALL print_array ( &
             grid% erlat(grid% ng1sm2:grid% ng1ep2,grid% ng2sm2:grid% ng2ep2,j1,jd), &
             grid% ng1sm2, grid% ng1ep2, grid% ng2sm2, grid% ng2ep2, &
             grid% ng1sm2, grid% ng1ep2, grid% ng2sm2, grid% ng2ep2, jd, 1, 0, &
             field_name(1:9), '-', 'loc_coor', 'h', &
             1000._wp, 0._wp,  kerror)
      ENDDO
#endif

    ENDDO     ! End loop over the diamonds

    ! 3. Compute the local spherical coordinates "peta" (local longi-
    !    tude) and "pchi" (local latitude) of the 6 (5) surrounding
    !    gridpoints for each node in diamond 1 (plus one extra row and
    !    column around it); compute the sine and cosine of the rotation
    !    angle ("pspsi", "pcpsi") between the local coordinate systems
    !    of the 6 (5) neighbouring gridpoints and the central node.
    !
    !    All computations are performed for diamond 1 only since the
    !    relative distances are the same for the other 9 diamonds except
    !    at the pole points of the diamonds (see section 8.)

    DO jm = 1,6          ! Loop over the 6 (5) neighbours
      DO j2 = grid% ng2sm1, grid% ng2ep1
        js2 = j2 + nspoke(jm+6)
        DO j1 = grid% ng1sm1, grid% ng1ep1
          js1 = j1 + nspoke(jm)

          ! Scalar product of the location vectors

          zd1 = grid% xn   (j1,j2,1,1)*grid% xn (js1,js2,1,1) +    &
                grid% xn   (j1,j2,2,1)*grid% xn (js1,js2,2,1) +    &
                grid% xn   (j1,j2,3,1)*grid% xn (js1,js2,3,1)

          ! Scalar product of the local longitude vector and the location
          ! vector

          zd2 = grid% erlon(j1,j2,1,1)*grid% xn (js1,js2,1,1) +  &
                grid% erlon(j1,j2,2,1)*grid% xn (js1,js2,2,1)

          ! Scalar product of the local latitude vector and the location
          ! vector

          zd3 = grid% erlat(j1,j2,1,1)*grid% xn (js1,js2,1,1) +  &
                grid% erlat(j1,j2,2,1)*grid% xn (js1,js2,2,1) +  &
                grid% erlat(j1,j2,3,1)*grid% xn (js1,js2,3,1)

          ! Avoid values GT. 1 and LT. -1 which are due to round off errors

          zd3 = SIGN ( MIN (1.0_wp, ABS (zd3)), zd3)

          ! Local spherical coordinates "eta" and "chi" for the array
          ! core; take care of pole points!

          grid% eta (j1,j2,jm) = ATAN (zd2/(zd1 + 1.e-20_wp))
          IF (j1 == 0 .AND. j2 == 1) grid% eta (j1,j2,jm) = ASIN (zd2)
          grid% chi (j1,j2,jm) = ASIN (zd3)

          ! Rotation angles "pspsi" and "pcpsi"
          ! Scalar product of the location vector and the local longitude
          ! vector at the surrounding gridpoints

          zd12 = grid% xn    (j1,j2,1,1)*grid% erlon(js1,js2,1,1) +  &
                 grid% xn    (j1,j2,2,1)*grid% erlon(js1,js2,2,1)

          ! Scalar product of the location vector and the local latitude
          ! vector at the surrounding gridpoints

          zd13 = grid% xn    (j1,j2,1,1)*grid% erlat(js1,js2,1,1) +  &
                 grid% xn    (j1,j2,2,1)*grid% erlat(js1,js2,2,1) +  &
                 grid% xn    (j1,j2,3,1)*grid% erlat(js1,js2,3,1)

          ! Scalar product of the local longitude vectors at the home node
          ! and at the surrounding gridpoints

          zd22 = grid% erlon (j1,j2,1,1)*grid% erlon(js1,js2,1,1) +  &
                 grid% erlon (j1,j2,2,1)*grid% erlon(js1,js2,2,1)

          ! Scalar product of the local longitude vector at the home node
          ! and the latitude vector at the surrounding gridpoints

          zd23 = grid% erlon (j1,j2,1,1)*grid% erlat(js1,js2,1,1) +  &
                 grid% erlon (j1,j2,2,1)*grid% erlat(js1,js2,2,1)

          ! Rotation angles "spsi" and "cpsi"

          grid% cpsi (j1,j2,jm) = COS (grid% eta(j1,j2,jm))*zd22 -  &
                                  SIN (grid% eta(j1,j2,jm))*zd12
          grid% spsi (j1,j2,jm) = COS (grid% eta(j1,j2,jm))*zd23 -  &
                                  SIN (grid% eta(j1,j2,jm))*zd13

        ENDDO
      ENDDO

    ENDDO     ! End of loop over the 6 (5) neighbours

    ! 4. Correct "eta", "chi", "cpsi", "spsi" at the four
    !    mirrored points of the extended diamonds

    DO js = 1, 4    ! Loop over the four mirrored points
      j1 = grid% ni1mrp(js)
      j2 = grid% ni2mrp(js)
      IF (j1 < -10) CYCLE  ! If we do not own the mirrored point

      DO jm = 1, 6  ! Loop over the six neigbours
        js1 = j1 + nspokes(jm  ,js)
        js2 = j2 + nspokes(jm+6,js)

        ! Scalar product of the location vectors

        zd1 = grid% xn   (j1,j2,1,1)*grid% xn (js1,js2,1,1) +  &
              grid% xn   (j1,j2,2,1)*grid% xn (js1,js2,2,1) +  &
              grid% xn   (j1,j2,3,1)*grid% xn (js1,js2,3,1)

        ! Scalar product of the local longitude vector and the location
        ! vector

        zd2 = grid% erlon(j1,j2,1,1)*grid% xn (js1,js2,1,1) +  &
              grid% erlon(j1,j2,2,1)*grid% xn (js1,js2,2,1)

        ! Scalar product of the local latitude vector and the location
        ! vector

        zd3 = grid% erlat(j1,j2,1,1)*grid% xn (js1,js2,1,1) +  &
              grid% erlat(j1,j2,2,1)*grid% xn (js1,js2,2,1) +  &
              grid% erlat(j1,j2,3,1)*grid% xn (js1,js2,3,1)

        ! Avoid values > 1 and < -1 which are due to round off errors

        zd3 = SIGN ( MIN (1.0_wp, ABS (zd3)), zd3)

        ! Local spherical coordinates "eta" and "chi"

        grid% eta (j1,j2,jm) = ATAN (zd2/(zd1 + 1.e-20_wp))
        grid% chi (j1,j2,jm) = ASIN (zd3)

        ! Rotation angles "spsi" and "cpsi"
        ! Scalar product of the location vector and the local longitude
        ! vector at the surrounding gridpoints

        zd12 = grid% xn    (j1,j2,1,1)*grid% erlon(js1,js2,1,1) +  &
               grid% xn    (j1,j2,2,1)*grid% erlon(js1,js2,2,1)

        ! Scalar product of the location vector and the local latitude
        ! vector at the surrounding gridpoints

        zd13 = grid% xn    (j1,j2,1,1)*grid% erlat(js1,js2,1,1) +  &
               grid% xn    (j1,j2,2,1)*grid% erlat(js1,js2,2,1) +  &
               grid% xn    (j1,j2,3,1)*grid% erlat(js1,js2,3,1)

        ! Scalar product of the local longitude vectors at the home node
        ! and at the surrounding gridpoints

        zd22 = grid% erlon (j1,j2,1,1)*grid% erlon(js1,js2,1,1) +  &
               grid% erlon (j1,j2,2,1)*grid% erlon(js1,js2,2,1)

        ! Scalar product of the local longitude vector at the home node
        ! and the latitude vector at the surrounding gridpoints

        zd23 = grid% erlon (j1,j2,1,1)*grid% erlat(js1,js2,1,1) +  &
               grid% erlon (j1,j2,2,1)*grid% erlat(js1,js2,2,1)

        ! Rotation angles "spsi" and "cpsi"

        grid% cpsi (j1,j2,jm) = COS (grid% eta(j1,j2,jm))*zd22 -  &
                                SIN (grid% eta(j1,j2,jm))*zd12
        grid% spsi (j1,j2,jm) = COS (grid% eta(j1,j2,jm))*zd23 -  &
                                SIN (grid% eta(j1,j2,jm))*zd13

        ! Copy results to the mirror point

        grid% eta (grid% ni1mrp(js+4),grid% ni2mrp(js+4),jm) = grid% eta (j1,j2,jm)
        grid% chi (grid% ni1mrp(js+4),grid% ni2mrp(js+4),jm) = grid% chi (j1,j2,jm)
        grid% cpsi(grid% ni1mrp(js+4),grid% ni2mrp(js+4),jm) = grid% cpsi(j1,j2,jm)
        grid% spsi(grid% ni1mrp(js+4),grid% ni2mrp(js+4),jm) = grid% spsi(j1,j2,jm)

      ENDDO  ! End of loop over the six neighbours

    ENDDO    ! End of loop over the four mirrored points

    ! 5. Compute the factor "bary" which is needed for the calculation
    ! of the barycentric coordinates of the gridpoints which are
    ! used for the interpolation routines

    DO jm = 1,6
      jk = jm + 1
      IF (jk == 7) jk = 1
      DO j2 = grid% ng2s, grid% ng2e
        DO j1 = grid% ng1s, grid% ng1e
          z1 = grid% eta(j1,j2,jm)*grid% chi(j1,j2,jk) - grid% eta(j1,j2,jk)*grid% chi(j1,j2,jm)
          IF (z1 == 0.0_wp) THEN
            grid% bary (j1,j2,jm) = 0.0_wp
          ELSE
            grid% bary (j1,j2,jm) = 1.0_wp/z1
          ENDIF
        ENDDO
      ENDDO
    ENDDO    ! End of loop over the 6 (5) neighbours

    ! 6. Repeat the value for the neighbour "1" (lower left) at index
    !    "7" for the arrays "eta" and "chi"

    DO j2 = grid% ng2sm1, grid% ng2ep1
      DO j1 = grid% ng1sm1, grid% ng1ep1
        grid% eta  (j1,j2,7) = grid% eta  (j1,j2,1)
        grid% chi  (j1,j2,7) = grid% chi  (j1,j2,1)
      ENDDO
    ENDDO

#ifdef DEBUG
    field_name  = 'eta(*)'
    DO j1 = 1,6
      WRITE (field_name(5:5),'(i1)') j1
      CALL print_array ( &
           grid% eta  (grid% ng1sm1:grid% ng1ep1,grid% ng2sm1:grid% ng2ep1,j1), &
           grid% ng1sm1, grid% ng1ep1, grid% ng2sm1, grid% ng2ep1, &
           grid% ng1sm1, grid% ng1ep1, grid% ng2sm1, grid% ng2ep1, 1, 1, 0, &
           field_name(1:9), 'degr.', 'loc_coor', 'h', &
           573._wp, 0._wp,  kerror)
    ENDDO

    field_name  = 'chi(*)'
    DO j1 = 1,6
      WRITE (field_name(5:5),'(i1)') j1
      CALL print_array ( &
           grid% chi  (grid% ng1sm1:grid% ng1ep1,grid% ng2sm1:grid% ng2ep1,j1), &
           grid% ng1sm1, grid% ng1ep1, grid% ng2sm1, grid% ng2ep1, &
           grid% ng1sm1, grid% ng1ep1, grid% ng2sm1, grid% ng2ep1, 1, 1, 0, &
           field_name(1:9), 'degr.', 'loc_coor', 'h', &
           573._wp, 0._wp,  kerror)
    ENDDO

    field_name  = 'bary(*)'
    DO j1 = 1,6
      WRITE (field_name(6:6),'(i1)') j1
      CALL print_array ( &
           grid% bary(grid% ng1s:grid% ng1e,grid% ng2s:grid% ng2e,j1), &
           grid% ng1s, grid% ng1e, grid% ng2s, grid% ng2e, &
           grid% ng1s, grid% ng1e, grid% ng2s, grid% ng2e, 1, 1, 0, &
           field_name(1:9), ' - ', 'loc_coor', 'h', &
           100._wp, 0._wp,  kerror)
    ENDDO
#endif

  END SUBROUTINE local_coordinates
  !---------------------------------------------------------------------------
  SUBROUTINE triangle_area (grid, kerror)

    USE mo_constants, ONLY: Re

    TYPE (triangular_grid), INTENT(inout) :: grid
    INTEGER,                INTENT(out)   :: kerror

    REAL(wp) :: zxne(3,3), & ! Cartesian coordinates of a triangle
                zt1,       &
                zt2,       &
                zt3,       &
                zs,        &
                zw

    LOGICAL :: lmask(grid% ng1s:grid% ng1ep1,grid% ng2sm1:grid% ng2e,2)

    INTEGER :: j1, j2      ! DO loop indices

    kerror = 0

    DO j2 = grid% ng2sm1,grid% ng2e
      DO j1= grid% ng1s,grid% ng1ep1

        ! Upward pointing triangles (index "1")

        zxne(1,1) = grid% xn (j1  ,j2  ,1,1)
        zxne(1,2) = grid% xn (j1  ,j2  ,2,1)
        zxne(1,3) = grid% xn (j1  ,j2  ,3,1)
        zxne(2,1) = grid% xn (j1-1,j2+1,1,1)
        zxne(2,2) = grid% xn (j1-1,j2+1,2,1)
        zxne(2,3) = grid% xn (j1-1,j2+1,3,1)
        zxne(3,1) = grid% xn (j1-1,j2  ,1,1)
        zxne(3,2) = grid% xn (j1-1,j2  ,2,1)
        zxne(3,3) = grid% xn (j1-1,j2  ,3,1)

        zt1 = zxne(1,1)*zxne(2,1) + zxne(1,2)*zxne(2,2) + zxne(1,3)*zxne(2,3)
        zt2 = zxne(2,1)*zxne(3,1) + zxne(2,2)*zxne(3,2) + zxne(2,3)*zxne(3,3)
        zt3 = zxne(3,1)*zxne(1,1) + zxne(3,2)*zxne(1,2) + zxne(3,3)*zxne(1,3)

        ! Avoid values > 1 and < -1 which are due to round off errors

        zt1 = SIGN ( MIN (1.0_wp, ABS (zt1)), zt1)
        zt2 = SIGN ( MIN (1.0_wp, ABS (zt2)), zt2)
        zt3 = SIGN ( MIN (1.0_wp, ABS (zt3)), zt3)

        zt1 = 0.5_wp*ACOS(zt1)
        zt2 = 0.5_wp*ACOS(zt2)
        zt3 = 0.5_wp*ACOS(zt3)
        zs  = 0.5_wp*(zt1 + zt2 + zt3)
        zw  = TAN(zs)*TAN(zs-zt1)*TAN(zs-zt2)*TAN(zs-zt3)
        grid% area(j1,j2,1) = 4*ATAN(SQRT(zw))

        ! Downward pointing triangles (index "2")

        zxne(1,1) = grid% xn (j1  ,j2  ,1,1)
        zxne(1,2) = grid% xn (j1  ,j2  ,2,1)
        zxne(1,3) = grid% xn (j1  ,j2  ,3,1)
        zxne(2,1) = grid% xn (j1  ,j2+1,1,1)
        zxne(2,2) = grid% xn (j1  ,j2+1,2,1)
        zxne(2,3) = grid% xn (j1  ,j2+1,3,1)
        zxne(3,1) = grid% xn (j1-1,j2+1,1,1)
        zxne(3,2) = grid% xn (j1-1,j2+1,2,1)
        zxne(3,3) = grid% xn (j1-1,j2+1,3,1)

        zt1 = zxne(1,1)*zxne(2,1) + zxne(1,2)*zxne(2,2) + zxne(1,3)*zxne(2,3)
        zt2 = zxne(2,1)*zxne(3,1) + zxne(2,2)*zxne(3,2) + zxne(2,3)*zxne(3,3)
        zt3 = zxne(3,1)*zxne(1,1) + zxne(3,2)*zxne(1,2) + zxne(3,3)*zxne(1,3)

        ! Avoid values > 1 and < -1 which are due to round off errors

        zt1 = SIGN ( MIN (1.0_wp, ABS (zt1)), zt1)
        zt2 = SIGN ( MIN (1.0_wp, ABS (zt2)), zt2)
        zt3 = SIGN ( MIN (1.0_wp, ABS (zt3)), zt3)

        zt1 = 0.5_wp*ACOS(zt1)
        zt2 = 0.5_wp*ACOS(zt2)
        zt3 = 0.5_wp*ACOS(zt3)
        zs  = 0.5_wp*(zt1 + zt2 + zt3)
        zw  = TAN(zs)*TAN(zs-zt1)*TAN(zs-zt2)*TAN(zs-zt3)
        grid% area(j1,j2,2) = 4*ATAN(SQRT(zw))

      ENDDO
    ENDDO

    ! Due to the 4 special points at the vertices of diamond #1 with
    ! only 5 neighbours, the following 6 triangles are undefined and
    ! their areas set to 0.

    IF (grid% ng1s   == grid% ngg1s   .AND. grid% ng2sm1 == grid% ngg2s-1) THEN
      grid% area(grid% ng1s  ,grid% ng2sm1,1) = 0.0_wp
      grid% area(grid% ng1s  ,grid% ng2sm1,2) = 0.0_wp
    ENDIF
    IF (grid% ng1ep1 == grid% ngg1e+1 .AND. grid% ng2e   == grid% ngg2e  ) THEN
      grid% area(grid% ng1ep1,grid% ng2e  ,1) = 0.0_wp
      grid% area(grid% ng1ep1,grid% ng2e  ,2) = 0.0_wp
    ENDIF
    IF (grid% ng1ep1 == grid% ngg1e+1 .AND. grid% ng2sm1 == grid% ngg2s-1) THEN
      grid% area(grid% ng1ep1,grid% ng2sm1,1) = 0.0_wp
    ENDIF
    IF (grid% ng1s   == grid% ngg1s   .AND. grid% ng2e   == grid% ngg2e  ) THEN
      grid% area(grid% ng1s  ,grid% ng2e  ,1) = 0.0_wp
    ENDIF

    WHERE (grid% area == 0.0_wp)
      lmask = .FALSE.
    ELSEWHERE
      lmask = .TRUE.
    ENDWHERE

    WRITE (6,'(a,f12.3,a)') &
         'triangle_area: minimum triangle area = ', &
         MINVAL(grid% area, MASK = lmask)*Re*Re*1e-6, ' km**2'
    WRITE (6,'(a,f12.3,a)') &
         'triangle_area: maximum triangle area = ', &
         MAXVAL(grid% area)*Re*Re*1e-6, ' km**2'

#ifdef DEBUG
    CALL print_array ( &
         grid% area (grid% ng1s:grid% ng1ep1,grid% ng2sm1:grid% ng2e,1), &
         grid% ng1s  , grid% ng1ep1, grid% ng2sm1, grid% ng2e  , &
         grid% ng1s  , grid% ng1ep1, grid% ng2sm1, grid% ng2e  ,  1, 1, 0, &
         'area(1)', ' -    ','tri_area', 'h', &
         1000._wp  , 0._wp,  kerror)
    CALL print_array ( &
         grid% area (grid% ng1s:grid% ng1ep1,grid% ng2sm1:grid% ng2e,2), &
         grid% ng1s  , grid% ng1ep1, grid% ng2sm1, grid% ng2e  , &
         grid% ng1s  , grid% ng1ep1, grid% ng2sm1, grid% ng2e  ,  1, 1, 0, &
         'area(2)', ' -    ','tri_area', 'h', &
         1000._wp  , 0._wp,  kerror)
#endif

  END SUBROUTINE triangle_area
  !---------------------------------------------------------------------------
  SUBROUTINE hexagon_area (grid, dc, kerror)

    USE mo_constants, ONLY: Re, Pi

    TYPE (triangular_grid), INTENT(inout) :: grid
    TYPE (t_atm_dec)      , INTENT(in)    :: dc
    INTEGER,                INTENT(out)   :: kerror

    REAL(wp) :: zrarn(grid% ng1s-2:grid% ng1e+2, grid% ng2s-2:grid% ng2e+2, grid% nd)

    INTEGER :: j1, j2, jd                  ! DO loop indices

    kerror = 0

    ! Compute the area of the hexagons for the core of diamond #1

    DO j2 = grid% ng2s,grid% ng2e
      DO j1= grid% ng1s,grid% ng1e

        grid% rarn(j1,j2) = 3.0_wp/( grid% area(j1+1,j2  ,1) + grid% area(j1  ,j2-1,2) +  &
                                     grid% area(j1  ,j2  ,1) + grid% area(j1  ,j2  ,2) +  &
                                     grid% area(j1+1,j2-1,1) + grid% area(j1+1,j2-1,2) )
      ENDDO
    ENDDO

    ! The 4 special points at the vertices of diamond #1 have only 5
    ! neighbours, therefore only 5 triangles have to be considered

    IF (grid% ng1s == grid% ngg1s .AND. grid% ng2s == grid% ngg2s) THEN
      grid% rarn(grid% ng1s, grid% ng2s) = 3.0_wp/(5*grid% area(grid% ng1s+1, grid% ng2s  ,1))
    ENDIF
    IF (grid% ng1e == grid% ngg1e .AND. grid% ng2s == grid% ngg2s) THEN
      grid% rarn(grid% ng1e, grid% ng2s) = 3.0_wp/(5*grid% area(grid% ng1e  , grid% ng2s  ,1))
    ENDIF
    IF (grid% ng1e == grid% ngg1e .AND. grid% ng2e == grid% ngg2e) THEN
      grid% rarn(grid% ng1e, grid% ng2e) = 3.0_wp/(5*grid% area(grid% ng1e  , grid% ng2e-1,2))
    ENDIF
    IF (grid% ng1s == grid% ngg1s .AND. grid% ng2e == grid% ngg2e) THEN
      grid% rarn(grid% ng1s, grid% ng2e) = 3.0_wp/(5*grid% area(grid% ng1s+1, grid% ng2e-1,1))
    ENDIF

    ! Now extend the array rarn by one row/column around the original
    ! core of diamond #1; set_boundaries will do this, it will be
    ! provided with 10 diamonds which are all the same

    DO jd = 1,grid% nd
      DO j2 = grid% ng2s, grid% ng2e
        DO j1 = grid% ng1s, grid% ng1e
          zrarn(j1,j2,jd) = grid% rarn(j1,j2)
        ENDDO
      ENDDO
    ENDDO

    CALL set_boundaries( zrarn, dc, &
         grid% ng1s-2, grid% ng1e+2, grid% ng2s-2, grid% ng2e+2, &
         1, 1, grid% nd, 2 )

    grid% rarn(:,:) = zrarn(grid% ng1s-1:grid% ng1e+1,grid% ng2s-1:grid% ng2e+1,1)

    WRITE (6,'(a,f12.3,a)') &
         'hexagon_area: minimum hexagon area   = ', &
         MINVAL(1e-6/grid% rarn*Re*Re), ' km**2'
    WRITE (6,'(a,f12.3,a)') &
         'hexagon_area: maximum hexagon area   = ', &
         MAXVAL(1e-6/grid% rarn*Re*Re), ' km**2'

    ! Calculate the weight assigned to each gridpoint (which is the
    ! area of the hexagon related to this point. At the pole point,
    ! the area is divided by 5 since the pole is contained in 5
    ! diamonds. Only the inner core of the diamond # 1 is taken to
    ! avoid counting gridpoints twice, i.e. for j1=0,ni and j2=ni+1
    ! as well as for j1=0 and j2=1,ni+1, the weight is set to 0.

    DO j2 = grid% ng2s, grid% ng2e
      DO j1 = grid% ng1s, grid% ng1e
        grid% hexwgt(j1,j2) = 1.0_wp/grid% rarn(j1,j2)
        IF (j1 == grid% ngg1s .OR.  j2 == grid% ngg2e) THEN
          grid% hexwgt(j1,j2) = 0.0_wp
        ENDIF
        IF (j1 == grid% ngg1s .AND. j2 == grid% ngg2s) THEN
          grid% hexwgt(j1,j2) = 1.0_wp/(5*grid% rarn(j1,j2))
        ENDIF
      ENDDO
    ENDDO

    WRITE (6,'(a,f12.3,a)') &
         'hexagon_area: global area            = ', &
         100*10*SUM(grid% hexwgt)/(4*Pi), ' %'


#ifdef DEBUG
    CALL print_array ( grid% rarn, &
         grid% ng1sm1, grid% ng1ep1, grid% ng2sm1, grid% ng2ep1, &
         grid% ng1sm1, grid% ng1ep1, grid% ng2sm1, grid% ng2ep1,  1, 1, 0, &
         'rarn', ' -    ','hex_area', 'h', &
         10._wp  , 0._wp,  kerror)
    CALL print_array ( grid% hexwgt, &
         grid% ng1s  , grid% ng1e  , grid% ng2s  , grid% ng2e  , &
         grid% ng1s  , grid% ng1e  , grid% ng2s  , grid% ng2e  ,  1, 1, 0, &
         'rhexwgt', ' -    ','hex_area', 'h', &
         1.e5_wp , 0._wp,  kerror)
#endif

  END SUBROUTINE hexagon_area
  !---------------------------------------------------------------------------
  SUBROUTINE differentiation_operators (grid, kerror)

    USE mo_constants, ONLY: Re

    TYPE (triangular_grid), INTENT(inout) :: grid
    INTEGER,                INTENT(out)   :: kerror

    INTEGER, PARAMETER :: nm = 6, nn = 5

#ifdef SVD
    REAL(wp) :: zc(nm,nn), zu(nm,nn), zut(nn,nm), zv(nn,nn), zvt(nn,nn), zw(nn)
    REAL(wp) :: zci(nn,nm), zwork(21)
#else
    REAL(wp) :: za(nn,nn) ! matrix for the calculation of the gradient op.
    REAL(wp) :: zb(nn)    ! r.h.s. of equation
#endif

#ifdef DEBUG2
    REAL(wp) :: psi(nm), psi_0 = 5
#endif

#ifdef SVD
    INTEGER :: j1, j2, jn, jm  ! Loop indices
#else
    INTEGER :: j1, j2, jm  ! Loop indices
#endif

    ! Preset gradient and Laplacian operators with 0

    grid% grd (:,:,:,:) = 0.0_wp
    grid% rlap(:,:,:)   = 0.0_wp

    ! Compute the gradient and Laplacian operators for the core of
    ! the diamond #1 plus one row and column around

#ifndef SVD
    DO j2 = grid% ng2sm1, grid% ng2ep1

      next_triangle: DO j1 = grid% ng1sm1, grid% ng1ep1

        ! Preset za with 0

        za(:,:) = 0.0_wp

        ! Loop over the 6 (5) neighbours

         next_neighbour: DO jm = 1, 6

          ! Take care of the vertices of diamond #1 with only 5 neighbours

          IF (j1 == grid% ngg1s .AND. j2 == grid% ngg2s .AND. jm == 5) CYCLE next_neighbour
          IF (j1 == grid% ngg1e .AND. j2 == grid% ngg2s .AND. jm == 6) CYCLE next_neighbour
          IF (j1 == grid% ngg1s .AND. j2 == grid% ngg2e .AND. jm == 4) CYCLE next_neighbour
          IF (j1 == grid% ngg1e .AND. j2 == grid% ngg2e .AND. jm == 1) CYCLE next_neighbour

          ! Take care of the two undefinded points in the extended array

          IF (j1 == grid% ngg1s-1 .AND. j2 == grid% ngg2s-1) CYCLE next_triangle
          IF (j1 == grid% ngg1e+1 .AND. j2 == grid% ngg2e+1) CYCLE next_triangle

          ! Only the upper triangle plus the diagonal of the symmetric
          ! matrix "za" has to be calculated

          za(1,1) = za(1,1) + grid% eta(j1,j2,jm)**2
          za(1,2) = za(1,2) + grid% eta(j1,j2,jm)   *grid% chi(j1,j2,jm)
          za(1,3) = za(1,3) + grid% eta(j1,j2,jm)**3
          za(1,4) = za(1,4) + grid% eta(j1,j2,jm)**2*grid% chi(j1,j2,jm)
          za(1,5) = za(1,5) + grid% eta(j1,j2,jm)   *grid% chi(j1,j2,jm)**2
          za(2,2) = za(2,2) +                        grid% chi(j1,j2,jm)**2
          za(2,5) = za(2,5) +                        grid% chi(j1,j2,jm)**3
          za(3,3) = za(3,3) + grid% eta(j1,j2,jm)**4
          za(3,4) = za(3,4) + grid% eta(j1,j2,jm)**3*grid% chi(j1,j2,jm)
          za(3,5) = za(3,5) + grid% eta(j1,j2,jm)**2*grid% chi(j1,j2,jm)**2
          za(4,5) = za(4,5) + grid% eta(j1,j2,jm)   *grid% chi(j1,j2,jm)**3
          za(5,5) = za(5,5) +                        grid% chi(j1,j2,jm)**4

        ENDDO  next_neighbour

        za(2,3) = za(1,4)
        za(2,4) = za(1,5)
        za(4,4) = za(3,5)

        ! Invert the matrix

        kerror = 0

        CALL dpotrf('U',nn,za,nn,kerror)

        IF (kerror /= 0) THEN
            WRITE (0,*) &
                 ' Error in subroutine differentiation_operators calling ', &
                 ' LAPACK/dpotrf'
          RETURN
        ENDIF

        ! Compute the r.h.s. of the equations

        DO jm = 1, 6

          ! Take care of the vertices of diamond #1 with only 5 neighbours

          IF (j1 == grid% ngg1s .AND. j2 == grid% ngg2s .AND. jm == 5) CYCLE
          IF (j1 == grid% ngg1e .AND. j2 == grid% ngg2s .AND. jm == 6) CYCLE
          IF (j1 == grid% ngg1s .AND. j2 == grid% ngg2e .AND. jm == 4) CYCLE
          IF (j1 == grid% ngg1e .AND. j2 == grid% ngg2e .AND. jm == 1) CYCLE

          zb(1) = grid% eta(j1,j2,jm)
          zb(2) = grid% chi(j1,j2,jm)
          zb(3) = grid% eta(j1,j2,jm)**2
          zb(4) = grid% eta(j1,j2,jm)*grid% chi(j1,j2,jm)
          zb(5) = grid% chi(j1,j2,jm)**2

          ! Solve linear system

          CALL dpotrs ('U',nn,1,za,nn,zb,nn,kerror)

          IF (kerror /= 0) THEN
            WRITE (0,*) &
                 ' Error in subroutine differentiation_operators calling ', &
                 ' LAPACK/dpotrs'
            RETURN
          ENDIF

          ! The gradient operator at the 6 surrounding nodes

          ! better /(Re*COS(chi(j1,j2,jm+1))) ???

          grid% grd (j1,j2,jm+1,1) = zb(1)/(Re*COS(grid% chi(j1,j2,jm+1)))
          grid% grd (j1,j2,jm+1,1) = zb(1)/Re
          grid% grd (j1,j2,jm+1,2) = zb(2)/Re

          ! The Laplacian operator at the 6 surrounding nodes

          ! better /(Re**2*COS(chi(j1,j2,jm+1))) ???

          grid% rlap(j1,j2,jm+1) = 2*(zb(3)+zb(5))/(Re*COS(grid% chi(j1,j2,jm+1)))**2
          grid% rlap(j1,j2,jm+1) = 2*(zb(3)+zb(5))/Re**2

        ENDDO
      ENDDO next_triangle
    ENDDO

#else

    DO j2 = grid% ng2sm1, grid% ng2ep1
      DO j1 = grid% ng1sm1, grid% ng1ep1

        zc(:,:) = 0.0_wp

        ! Take care of the two undefinded points in the extended array

        IF (j1 == grid% ngg1s-1 .AND. j2 == grid% ngg2s-1) CYCLE
        IF (j1 == grid% ngg1e+1 .AND. j2 == grid% ngg2e+1) CYCLE

        ! Build the design matrix taking into account that the
        ! local coordinates of neighbour 6 for pentagons are 0.

        IF (j1 == grid% ngg1s .AND. j2 == grid% ngg2s ) THEN
          grid% eta(j1,j2,5) = 0.0_wp
          grid% chi(j1,j2,5) = 0.0_wp
        ELSE IF (j1 == grid% ngg1e .AND. j2 == grid% ngg2s ) THEN
          grid% eta(j1,j2,6) = 0.0_wp
          grid% chi(j1,j2,6) = 0.0_wp
        ELSE IF (j1 == grid% ngg1s .AND. j2 == grid% ngg2e ) THEN
          grid% eta(j1,j2,4) = 0.0_wp
          grid% chi(j1,j2,4) = 0.0_wp
        ELSE IF (j1 == grid% ngg1e .AND. j2 == grid% ngg2e ) THEN
          grid% eta(j1,j2,1) = 0.0_wp
          grid% chi(j1,j2,1) = 0.0_wp
        END IF

        zc(1:6,1) = grid% eta(j1,j2,1:6)
        zc(1:6,2) = grid% chi(j1,j2,1:6)
        zc(1:6,3) = grid% eta(j1,j2,1:6)**2
        zc(1:6,4) = grid% eta(j1,j2,1:6)*grid% chi(j1,j2,1:6)
        zc(1:6,5) = grid% chi(j1,j2,1:6)**2

        ! Solve SVD problem

        CALL dgesvd ('A','A',nm,nn,zc,nm,zw,zu,nm,zvt,nn,zwork,21,kerror)

        IF (kerror /= 0) THEN
          WRITE (0,*) &
               ' Error in subroutine differentiation_operators calling ', &
               ' LAPACK/dgesvd'
          RETURN
        ENDIF

        ! Transpose decomposed matrices

        zut = TRANSPOSE(zu)
        zv  = TRANSPOSE(zvt)

        ! Multiply (U-transpose) 1/w V for the inverse of the design matrix
        ! the columns of this matrix are the operator coefficients

        DO jn = 1, nn
          zv(:,jn) = zv(:,jn) / zw(jn)
        END DO
        zci = MATMUL ( zv, zut)

        ! Take care of the vertices of diamond #1 with only 5 neighbours

        IF (j1 == grid% ngg1s .AND. j2 == grid% ngg2s) zci(:,5) = 0.0_wp
        IF (j1 == grid% ngg1e .AND. j2 == grid% ngg2s) zci(:,6) = 0.0_wp
        IF (j1 == grid% ngg1s .AND. j2 == grid% ngg2e) zci(:,4) = 0.0_wp
        IF (j1 == grid% ngg1e .AND. j2 == grid% ngg2e) zci(:,1) = 0.0_wp

        grid% grd(j1,j2,2:7,1) = zci(1,1:6)/(Re*COS(grid% chi(j1,j2,2:7)))
!        grid% grd(j1,j2,2:7,1) = zci(1,1:6)/Re
        grid% grd(j1,j2,2:7,2) = zci(2,1:6)/Re
        grid% rlap(j1,j2,2:7) = 2*zci(3,1:6)/(Re*COS(grid% chi(j1,j2,2:7)))**2 &
                               +2*zci(5,1:6)/(Re*Re)                     &
                               -zci(2,1:6)*TAN(grid% chi(j1,j2,2:7))/(Re*Re)
!        grid% rlap(j1,j2,2:7) = 2*(zci(3,1:6)+zci(5,1:6))/(Re*Re)

        ! Recover local coordinates in original setup

        IF (j1 == grid% ngg1s .AND. j2 == grid% ngg2s ) THEN
          grid% eta(j1,j2,5) = grid% eta(j1,j2,4)
          grid% chi(j1,j2,5) = grid% chi(j1,j2,4)
        ELSE IF (j1 == grid% ngg1e .AND. j2 == grid% ngg2s ) THEN
          grid% eta(j1,j2,6) = grid% eta(j1,j2,5)
          grid% chi(j1,j2,6) = grid% chi(j1,j2,5)
        ELSE IF (j1 == grid% ngg1s .AND. j2 == grid% ngg2e ) THEN
          grid% eta(j1,j2,4) = grid% eta(j1,j2,3)
          grid% chi(j1,j2,4) = grid% chi(j1,j2,3)
        ELSE IF (j1 == grid% ngg1e .AND. j2 == grid% ngg2e ) THEN
          grid% eta(j1,j2,1) = grid% eta(j1,j2,7)
          grid% chi(j1,j2,1) = grid% chi(j1,j2,7)
        END IF

      ENDDO
    ENDDO

#endif

    ! Compute the gradient and Laplacian terms associated with the
    ! 'home' node.

    DO jm = 1, 6
      grid% grd (:,:,1,1) = grid% grd (:,:,1,1) - grid% grd (:,:,jm+1,1)
      grid% grd (:,:,1,2) = grid% grd (:,:,1,2) - grid% grd (:,:,jm+1,2)
      grid% rlap(:,:,1  ) = grid% rlap(:,:,1  ) - grid% rlap(:,:,jm+1  )
    ENDDO


#ifdef DEBUG2
    ! Check results, assume gradients are 1.

    WRITE (0,'(/,a)') 'Gradient weight check:'
    DO j2 = grid% ng2s, grid% ng2e
      DO j1 = grid% ng1s, grid% ng1e
        DO jn = 1, nm
          psi(jn) = psi_0 + grid% eta(j1,j2,jn) + grid% chi(j1,j2,jn) &
                          + grid% eta(j1,j2,jn)**2 + grid% eta(j1,j2,jn)*grid% chi(j1,j2,jn) &
                          + grid% chi(j1,j2,jn)**2
        ENDDO
        WRITE (0,'(2i4,3x,a,f10.7)') &
             j1, j2, 'dpsideta (= 1 ?) = ', &
             Re*SUM(grid% grd(j1,j2,2:7,1)*(psi(1:6)-psi_0))
        WRITE (0,'(2i4,3x,a,f10.7)') &
             j1, j2, 'dpsidchi (= 1 ?) = ', &
             Re*SUM(grid% grd(j1,j2,2:7,2)*(psi(1:6)-psi_0))
      ENDDO
    ENDDO
    WRITE (0,*)
#endif

#ifdef DEBUG
    field_name = 'pgrd(*,*)'
    DO j1 = 1,2
      WRITE (field_name(8:8), '(i1)') j1
      DO jm = 1,7
        WRITE (field_name(6:6), '(i1)') jm

        CALL print_array ( &
             grid% grd(grid% ng1sm1: grid% ng1ep1,grid% ng2sm1:grid% ng2ep1,jm,j1), &
             grid% ng1sm1, grid% ng1ep1, grid% ng2sm1, grid% ng2ep1, &
             grid% ng1sm1, grid% ng1ep1, grid% ng2sm1, grid% ng2ep1, 1, 1, 0, &
             field_name(1:9), '1/m ','loc_coor', 'h', &
             1.e10_wp, 0._wp,  kerror)
      ENDDO
    ENDDO

    field_name = 'prlap(*)'
    DO jm = 1,7
      WRITE (field_name(7:7), '(i1)') jm

      CALL print_array ( &
           grid% rlap(grid% ng1sm1:grid% ng1ep1,grid% ng2sm1:grid% ng2ep1,jm), &
           grid% ng1sm1, grid% ng1ep1, grid% ng2sm1, grid% ng2ep1, &
           grid% ng1sm1,grid%  ng1ep1, grid% ng2sm1, grid% ng2ep1, 1, 1, 0, &
           field_name(1:9), 'm**2/s ','loc_coor', 'h', &
           1.e15_wp, 0._wp,  kerror)
    ENDDO
#endif

    kerror = 0

  END SUBROUTINE differentiation_operators
  !---------------------------------------------------------------------------
  SUBROUTINE ll2xyz(long, lat, x, y, z)

    ! To convert longitude and latitude to cartesian coordinates
    ! on the unit sphere

    USE mo_kind, ONLY: wp

    REAL (wp), INTENT(in)  :: long, lat
    REAL (wp), INTENT(out) :: x, y, z

    REAL (wp) :: cln, sln, clt, slt

    sln = SIN(long)
    cln = COS(long)
    slt = SIN(lat)
    clt = COS(lat)

    x = cln*clt
    y = sln*clt
    z = slt

  END SUBROUTINE ll2xyz
  !---------------------------------------------------------------------------
  PURE SUBROUTINE xyz2ll(x,y,z,long,lat)

    ! To convert cartesian coordinates to longitude and latitude

    USE mo_kind,      ONLY: wp
    USE mo_constants, ONLY: pi

    REAL (wp), INTENT(in)  :: x, y, z
    REAL (wp), INTENT(out) :: long, lat

    REAL (wp) :: tlt, r !, tln

    IF (x == 0.0_wp) THEN
      IF (y >= 0.0_wp) THEN
        long = 0.5_wp*pi
      ELSE
        long = 1.5_wp*pi
      END IF
    ELSE
      ! Work around pgi runtime bug with atan(0):
      long = atan2 (y, x)
      ! Original code:
!     tln = y/x
!     long = ATAN(tln)
!     IF (x < 0.0_wp) THEN
!       long = long + pi
!     END IF
      IF (long < 0.0_wp) THEN
        long = long + 2.0_wp*pi
      END IF
    END IF

    r = SQRT(x*x+y*y)
    IF (r == 0.0_wp) THEN
      IF (z > 0.0_wp) THEN
        lat = 0.5_wp*pi
      ELSE
        lat = -0.5_wp*pi
      END IF
    ELSE
      tlt = z/r
      lat = ATAN(tlt)
    END IF

  END SUBROUTINE xyz2ll
  !---------------------------------------------------------------------------
  SUBROUTINE print_array (px, kil1s, kil1e, kil2s, kil2e, &
       ki1sc, ki1ec, ki2sc, ki2ec, kjd, kj3, &
       kntstep, yvarnam, yvardim, yscall , ypform, &
       pfactor, pbias, kerror)

    USE mo_kind, ONLY: wp

    INTEGER, INTENT(in) :: kil1s, kil1e, kil2s, kil2e, kjd, kj3
    INTEGER, INTENT(in) :: kntstep, ki1sc, ki1ec, ki2sc, ki2ec

    REAL(wp), INTENT(in) :: px(kil1s:kil1e, kil2s:kil2e)
    REAL(wp), INTENT(in) :: pfactor, pbias

    INTEGER, INTENT(out) :: kerror

    CHARACTER(len=*) :: yvarnam, yvardim, yscall, ypform

    INTEGER :: mpxscal (23)   ! auxiliary vector with scaled values

    INTEGER :: mi1s, mi1e, & ! local start and end addresses of printing
               mi2s, mi2e, & ! in the first and second array dimension
               j1, j2        ! DO loop variables

    REAL(wp) :: zpxmin, zpxmax      ! Minimum/maximum value of field "px"

    kerror = 0

    WRITE (0,*)  'SUBROUTINE print_array, variable: ', yvarnam, ' callee: ', yscall

    ! Check the print ranges in relation to the array dimensions

    IF (ki1sc < kil1s) THEN
      WRITE (0,*) 'print_array error, ki1sc= ', ki1sc, &
           '  is less than kil1s= ', kil1s
      kerror = -1
      RETURN
    ELSE IF (ki2sc < kil2s) THEN
      WRITE (0,*) 'print_array error, ki2sc= ', ki2sc, &
           '  is less than kil2s= ', kil2s
      kerror = -1
      RETURN
    ELSE IF (ki1ec > kil1e) THEN
      WRITE (0,*) 'print_array error, ki1ec= ', ki1ec, &
           '  is less than kil1e= ', kil1e
      kerror = -1
      RETURN
    ELSE IF (ki2ec > kil2e) THEN
      WRITE (0,*) 'print_array error, ki2ec= ', ki2ec, &
           '  is less than kil2e= ', kil2e
      kerror = -1
      RETURN
    ENDIF

    ! Calculate minimum and maximum value of field

    zpxmin = MINVAL (px(ki1sc:ki1ec,ki2sc:ki2ec))
    zpxmax = MAXVAL (px(ki1sc:ki1ec,ki2sc:ki2ec))

    ! Horizontal cross section, i.e. 2-dimensional field

    IF (ypform == 'h') THEN
      WRITE(0,'(a,i3)') ' Horizontal cross section for diamond #: ', kjd
      WRITE(0,'(a,i3,a,i8)') ' level number: ', kj3, '   time step: ', kntstep

      ! Print diamond 1 to 5 from north to south

      IF (kjd < 6) THEN
        mi2s = ki2sc

        DO
          WRITE (0,'(4a)') ' variable: ', yvarnam,'  unit: ', yvardim
          WRITE (0,'(2(a,e10.4),/)') ' factor: ', pfactor, '  bias: ', pbias
          WRITE (0,'(2(a,e10.4),/)') ' minimum: ', zpxmin,  '  maximum: ', zpxmax
          mi2e = MIN (mi2s+22, ki2ec)
          WRITE (0,'(a, 23i5,//)') '   j1, j2',(j2, j2=mi2s, mi2e)
          IF (MOD (mi2e, 23) /= 0) WRITE (0,'(a,/)') ' '
          DO j1 = ki1sc, ki1ec
            DO j2 = mi2s, mi2e
              mpxscal(j2-mi2s+1) = NINT ((px(j1,j2)-pbias)*pfactor)
            ENDDO
            WRITE (0,'(i5, 5x, 23i5)') j1, (mpxscal(j2), j2 = 1, mi2e-mi2s+1)
          ENDDO
          mi2s = mi2s + 23
          IF (mi2s > ki2ec) THEN
            WRITE (0,'(78("="))')
            RETURN
          ENDIF
        ENDDO

        ! Print diamond 6 to 10 from south to north

      ELSE IF (kjd > 5) THEN
        mi2s = ki2sc

        DO
          WRITE (0,'(4a)') '  variable: ', yvarnam,'  unit: ', yvardim
          WRITE (0,'(2(a,e10.4),/)') '  factor: ', pfactor, '   bias: ', pbias
          WRITE (0,'(2(a,e10.4),/)') ' minimum: ', zpxmin, '  maximum: ', zpxmax
          mi2e = MIN (mi2s+22, ki2ec)
          WRITE (0,'(a, 23i5,//)') '   j1, j2',(j2, j2=mi2s, mi2e)
          IF (MOD (mi2e, 23) /= 0) WRITE (0,'(a,/)') ' '
          DO j1 = ki1ec, ki1sc, -1
            DO j2 = mi2s, mi2e
              mpxscal(j2-mi2s+1) = NINT ((px(j1,j2)-pbias)*pfactor)
            ENDDO
            WRITE (0,'(i5, 5x, 23i5)') j1, (mpxscal(j2), j2 = 1, mi2e-mi2s+1)
          ENDDO
          mi2s = mi2s + 23
          IF (mi2s .GT. ki2ec) THEN
            WRITE (0,'(78("="))')
            RETURN
          ENDIF
        ENDDO
      ENDIF

      ! Vertical cross section

    ELSE IF (ypform == 'v') THEN

      WRITE (0,*) ' Vertical cross section for diamond #: ', kjd
      mi1s = ki1sc

      DO
        WRITE (0,'(4a)') '  variable: ', yvarnam,'  unit: ', yvardim
        WRITE (0,'(2(a,e10.4),/)') '  factor: ', pfactor, '   bias: ', pbias
        WRITE (0,'(2(a,e10.4),/)') ' minimum: ', zpxmin, '  maximum: ', zpxmax
        mi1e = MIN (mi1s+22, ki1ec)
        WRITE (0,'(a, 23i5,//)') '    k, j1',(j1, j1=mi1s, mi1e)
        IF (MOD (mi1e, 23) /= 0) WRITE (0,'(a,/)') ' '
        DO j2 = ki2sc, ki2ec
          DO j1 = mi1s, mi1e
            mpxscal(j1-mi1s+1) = NINT ((px(j1,j2)-pbias)*pfactor)
          ENDDO
          WRITE (0,'(i5, 5x, 23i5)') j2, (mpxscal(j1), j1 = 1, mi1e-mi1s+1)
        ENDDO
        mi1s = mi1s + 23
        IF (mi1s > ki1ec) THEN
          WRITE (0,'(78("="))')
          RETURN
        ENDIF
      ENDDO
    ENDIF

  END SUBROUTINE print_array

END MODULE mo_ico_grid
