!
!+ GNSS Radio occultation operator: Generation of icosahedral-hexagonal grid
!
MODULE ICO_grid
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Generation of icosahedral-hexagonal grid.
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
! V1_6         2009/06/10 Harald Anlauf
!  Fix precision of FP constants
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_22        2013-02-13 Harald Anlauf
!  add vertical coordinate type, geopotential height
! V1_27        2013-11-08 Harald Anlauf
!  for ICON: horizontal interpolation, remove ak,bk, waterload
! V1_42        2015-06-08 Harald Anlauf
!  horint_mode
! V1_43        2015-08-19 Harald Anlauf
!  preparations for alternative refractivity models
! V1_46        2016-02-05 Andreas Rhodin
!  base decisions on new flag 'vct', not 'ivctype'
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Reference:
!   Michael E. Gorbunov and Luis Kornblueh
!   Principles of variational assimilation of GNSS radio occultation data.
!   Max-Planck-Institut fuer Meteorologie, Hamburg, Report No. 350 (2003)
!
! Author:
! Michael E. Gorbunov  2004  original code
! Changes:
! Andreas Rhodin             adapted to DWD 3D-VAR
!==============================================================================
!
! Module ICO_grid
!
! Generation of icosahedral-hexagonal grid.
!----------------------------------------------------------
! (C) Copyright 2001-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 29 Nov 2001 | Adaptation for use
!         |             | with ECHAM library.
!   2.0   | 14 Apr 2002 | Multiple_Points.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, Pi, &
    WorkPr
!
Use mo_exception, only: &
    finish

!
Use ICO_boundary, only: &
! Imported Routines:
    Set_Boundaries
!
Use mo_atm_grid, only: &
    t_grid             ! Derived type for grid metadata
!----------------------------------------------------------
Implicit None
Private
!----------------------------------------------------------
! Public Parameters:

public :: gf       ! global fields
public :: t_global ! global fields data type
public :: t_single ! component of t_global (single level fields)
public :: nspoke   ! Offsets of 6(5) neighbours
!
Integer, Parameter :: &
   nd = 10            ! Number of diamonds
!
Integer, Parameter :: &
   nspoke(2,6) = &    ! Offsets of 6(5) neighbours
                      ! [1:6,j1/j2]
   Reshape(Source = (/ 1,  0, -1, -1,  0,  1,     &
                       0,  1,  1,  0, -1, -1 /),  &
           Shape  = (/ 2, 6 /), Order = (/ 2, 1 /) )
!
! Public Scalars:
!
Integer :: ni, ni2, nir       ! ni=2**ni2*nir triangles / diamond side
                              ! nir = 1 or 3
!
! --- First and second index intervals for diamond
!
 Integer :: ng1s,   ng1e       ! Core interval of 1st index 0..ni
!Integer :: ng1sm1, ng1ep1     ! 1-enlarged index interval -1..ni+1
 Integer :: ng1sm2, ng1ep2     ! 2-enlarged index interval -2..ni+2]
!Integer :: ngg1s,  ngg1e      ! = ng1s, ng1e (for compatibility)
 Integer :: ng2s,   ng2e       ! Core interval of 2nd index 1..ni+1
!Integer :: ng2sm1, ng2ep1     ! 1-enlarged index interval [ 0:ni+2]
 Integer :: ng2sm2, ng2ep2     ! 2-enlarged index interval [-1:ni+3]
!Integer :: ngg2s,  ngg2e      ! = ng2s, ng2e (for compatibility)
!
! Public Arrays:
!
!Integer :: &
!   ni1mrp(8), ni2mrp(8)  ! Indices of the mirrored points

Real(Double), Allocatable :: &
   xnglob(:,:,:,:),    & ! Cartesian  coordinates of gridpoints on unit-sphere.
                         !   rank 1: i index on ith diamond [0:ni]
                         !   rank 2: j index on ith diamond [1:ni+1]
                         !   rank 3: x, y, z coordinates    [1:3]
                         !   rank 4: number of diamond      [1:10]
   rlon(:,:,:),        & ! Gridpoint longitudes [rad]
   rlat(:,:,:),        & ! Gridpoint latitudes [rad]
   xn(:,:,:,:)           ! Local vertical unit vector = node vector
!----------------------------------------------------------
!
!==============================================================================
#if   defined (_CRAYFTN) \
  || (defined (__GFORTRAN__) && ((__GNUC__ * 100 + __GNUC_MINOR__) >= 406)  \
                             && ((__GNUC__ * 100 + __GNUC_MINOR__) <  800)) \
  || (defined (__INTEL_COMPILER) && (__INTEL_COMPILER >= 1500)) \
  || (defined (__NEC__) && (__NEC_VERSION__ >= 30004)) \
  || (defined (__PGI) && (__PGIC__ >= 15) && (__PGIC__ != 17))
#define _CONTIGUOUS         ,Contiguous
#define _POINTER    ,Pointer,Contiguous
#else
#define _CONTIGUOUS
#define _POINTER    ,Pointer
#endif
!==============================================================================

!--------------------
! single level fields
!--------------------
Type t_single
   integer      :: i1    ! first  (horizontal) index in global field
   integer      :: i2    ! second (horizontal) index in global field
   integer      :: id    ! fourth (diamond)    index in global field
   real(workpr) :: Hsur  ! Surface geopotential [gpm]
   real(workpr) :: Gundu ! geoid undulation       [m]
   real(workpr) :: Psur  ! Surface pressure      [Pa]
   real(workpr) :: XLon  ! Longitude grid       [deg]
   real(workpr) :: XLat  ! Latitude grid        [deg]
   real(Double) :: GCLat ! Geocentric latitudes [rad]
   real(Double) :: Zsur  ! Surface altitude      [km]
end Type t_single

!----------------------------------------------------
! global fields (single level and multi level fields)
!----------------------------------------------------
Type t_global
   !----------
   ! meta data
   !----------
   Integer                 :: Grid_Type =  0     ! Horizontal grid type
   Integer                 :: vctype    =  0     ! Vertical coordinate type
   Integer                 :: ngp       =  0     ! Subgrid dimension
   Integer                 :: ncol      =  0     ! used columns/profile
   Integer                 :: hint_mode =  1     ! Horizontal interpolation mode
   Integer                 :: ref_model =  1     ! Refractivity model
   Integer                 :: n         =  0     ! number of columns stored
   Integer                 :: nz        =  0     ! number of levels
   Integer                 :: nd        =  0     ! number of diamonds
   Integer                 :: ni        =  0     !
   Integer                 :: ni2       =  0     !
   Integer                 :: nir       =  1     !
   Integer                 :: lbg(3)    =  0     ! lower bound of fields
   Integer                 :: ubg(3)    = -1     ! upper ... (x1,x2,diam)
   Integer                 :: yyyymmdd  =  0     ! date
   Integer                 :: hhmmss    =  0     ! time
   Real(Double)            :: p0_msis   =  0     ! MSIS base pressure level [Pa]
   Real(Double)            :: ref_scale =  1     ! refractivity scaling factor
   integer        _POINTER :: i (:,:,:) =>NULL() ! index array (i1,i2,id)
   real(workpr)   _POINTER :: glon  (:) =>NULL() ! Longitude(only Gaussgr)[deg]
   real(workpr)   _POINTER :: glat  (:) =>NULL() ! Latitude (only Gaussgr)[deg]
!  Real(Double)   _POINTER :: A     (:) =>NULL() ! Vertical coordinate coeffs.
!  Real(Double)   _POINTER :: B     (:) =>NULL() ! Vertical coordinate coeffs.
   Real(Double)   _POINTER :: B_ad  (:) =>NULL() ! B used for adjoints
   Real(Double)   _POINTER :: Z   (:,:) =>NULL() ! Geopot.height(half lev.) [m]
   Real(Double)   _POINTER :: xnglob(:,:,:,:) =>NULL() ! cartesian coordinate
   integer        _POINTER :: marr  (:,:,:,:) =>NULL() ! (4,i,j,d) index field
   type(t_grid)   ,pointer :: grid      =>NULL() ! Grid metadata (e.g. for ICON)
   !-----
   ! data
   !-----
   Type(t_single) _POINTER :: s   (:)   =>NULL() ! single level fields
   real(workpr)   _POINTER :: T (:,:)   =>NULL() ! Temperature              [K]
   real(workpr)   _POINTER :: Q (:,:)   =>NULL() ! Specific humidity    [kg/kg]
   real(workpr)   _POINTER :: Qx(:,:)   =>NULL() ! Condensate(water/ice)[kg/kg]
   real(workpr)   _POINTER :: P (:,:)   =>NULL() ! Pressure(full levels)   [Pa]
   real(workpr)   _POINTER :: Ph(:,:)   =>NULL() ! Pressure(half levels)   [Pa]
end Type t_global

type (t_global) ,save :: gf



Contains


!==========================================================
Subroutine Factorize_ni &
  (ErrStat)       ! --> Error status
!
! Factorization of ni as 2**ni2*nir
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 27 Nov 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Output arguments:
!
Integer, Intent(Out) :: ErrStat  ! Error status
                                 !  0 - ni factorized
                                 ! -1 - ni cannot be factorized
!----------------------------------------------------------
! Local Scalars:
!
Integer  :: mx        ! auxiliary variable, set to kni initially
!----------------------------------------------------------
! Global variables used:
!
!  ni, ni2, nir       ! ni=2**ni2*nir triangles / diamond side
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------

ErrStat = 0

mx      = ni
ni2     = 0
nir     = 1


!----------------------------------------------------------
! 2. FACTORIZATION
!----------------------------------------------------------

Do While (mx > 1)
   If (MOD(mx,2) == 0) then
      ni2  = ni2+1
      mx    = mx/2
   Else
      nir   = mx
      mx    = 1
   End If
End Do


!----------------------------------------------------------
! 3. CHECK NIR
!----------------------------------------------------------

If (nir > 3) then
   ErrStat = -1
End If

End Subroutine Factorize_ni



!==========================================================
Subroutine Triangle_Center &
  (jd)     ! <-- Diamond number
!
! Computation of coordinates of icosahedral
! triangle centers.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 27 Nov 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(in)  :: &
   jd       ! Diamond number
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: zxnorm    ! Vector norm
Integer      :: j         ! Loop index
Integer      :: mi1       ! Index of center point
Integer      :: mi2       ! Index of top or bottom diamond corner
!----------------------------------------------------------
! Global variables used:
!
! xnglob, ni, nd
!----------------------------------------------------------


!----------------------------------------------------------
! 1. COMPUTATION OF TRIANGLE CENTER
!----------------------------------------------------------

Do j = 1, 2

   !--- 1.1. Setting indexes

   mi1  = j*ni/3
   mi2  = 1 + (j - 1)*ni

   !--- 1.2. Computation of vector sum

   xnglob(mi1,mi1+1,:,jd) = xnglob(mi2-1, mi2, :, jd) +   &
                            xnglob(ni,   1,    :, jd) +   &
                            xnglob(0,    ni+1, :, jd)

   !--- 1.3. Normalize to unit-sphere

   zxnorm = 1.0_Double/Sqrt(Sum(xnglob(mi1,mi1+1,:,jd)**2))
   xnglob(mi1,mi1+1,:,jd) = zxnorm*xnglob(mi1,mi1+1,:,jd)

End Do


End Subroutine Triangle_Center



!==========================================================
Subroutine Global_Coordinates  &
  (ErrStat)      ! --> Error status
!
! Computation of icosahedric triangular grid.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 29 Nov 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Output arguments:
!
Integer, Intent(Out) :: &
   ErrStat   ! Error status
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: zw      ! the spherical angle in an icosahedron
                        ! subtended by two vertices.
Real(Double) :: zcosw   ! cosine(zw)
Real(Double) :: zsinw   ! sine  (zw)
Real(Double) :: zsgn    ! zsgn is a hemisphere factor.
                        ! zsgn =  1.0 is north  (diamonds 1- 5)
                        ! zsgn = -1.0 is south  (diamonds 6-10)
Real(Double) :: zrlon   ! longitude of diamond vertices
Real(Double) :: zgamma  ! fraction of great circle angle
Real(Double) :: zchord  ! Cartesian distance between two points
Real(Double) :: ztheta  ! Great circle angle between two points
Real(Double) :: zalpha  ! Weighting factor
Real(Double) :: zbeta   ! Weighting factor
Integer :: mcosv(nd)    ! meridian angle locations of the 10
                        ! non-polar vertices in units of Pi/5
Integer :: ml           ! recursive index interval
Integer :: ml2          ! recursive bisected index interval
Integer :: ml3          ! trisected index interval
Integer :: mi1          ! recursive row index of new node
Integer :: mi2          ! recursive column index of new node
Integer :: mm           ! recursive number of subdivisions
Integer :: j1, j2, jd   ! Gridpoint indices
Integer :: jb           ! Bisecting interval index
!----------------------------------------------------------
! Global variables used:
!
!   xnglob(:,:,:,:)     ! Cartesian  coordinates of gridpoints
!   rlon(:,:,:)         ! Gridpoint longitudes [rad]
!   rlat(:,:,:)         ! Gridpoint latitudes [rad]
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------

ErrStat = 0


!----------------------------------------------------------
! 1. GLOBAL COORDINATE COMPUTATION
!----------------------------------------------------------


!--- 1.1. Compute angles associated with the icosahedron.

zw      = 2*ACOS(1.0_Double/(2*SIN(Pi/5.0_Double)))
zcosw   = COS(zw)
zsinw   = SIN(zw)


!--- 1.2. Compute meridian angle locations

Do jd = 1, nd
   If (MOD(jd,2) == 1) then
      mcosv((jd+1)/2) = -1 + (jd - 1) - nd*((jd - 1)/7)
   Else
      mcosv(jd/2+5)   = -1 + (jd - 1) - nd*((jd - 1)/7)
   End If
End Do

! Loop over the ten diamonds computing diamond vertices (x,y,z)
! coordinates and then iteratively bisecting them ni2 times.
! First a trisection is performed, if required (nir=3).


!--- 1.3. Computing global coordinates for diamonds


Diamonds: Do jd = 1, nd

   !--- 1.3.1. Toggle hemisphere

   If (jd >= 6) then
      zsgn = -1.0_Double    ! southern
   Else
      zsgn =  1.0_Double    ! northern
   End If

   !--- 1.3.2. Compute meridian angle for diamond home vertices

   zrlon = mcosv(jd)*Pi/5.0_Double

   !--- 1.3.3. Compute polar vertices (0,1,jd)

   xnglob(0,1,1:3, jd) =  &
     (/ 0.0_Double, 0.0_Double, zsgn /)

   !--- 1.3.4. Compute home vertices (ni,1,jd).

   xnglob(ni,1,1:3,jd) =  &
     (/ zsinw*COS(zrlon), zsinw*SIN(zrlon), zcosw*zsgn /)

   !--- 1.3.5. Compute vertex with same latitude as home

   xnglob(0,ni+1,1:3,jd) = &
     (/ zsinw*COS(zrlon + 2*Pi/5.0_Double), &
        zsinw*SIN(zrlon + 2*Pi/5.0_Double), &
        zcosw*zsgn /)

   !--- 1.3.6. Compute las vertex, which in opposite hemisphere (ni,ni+1,,)

   xnglob(ni,ni+1,1:3,jd) =  &
     (/ zsinw*COS(zrlon + Pi/5.0_Double), &
        zsinw*SIN(zrlon + Pi/5.0_Double), &
       -zcosw*zsgn /)

   !--- 1.3.7. Trisection if required (ni3=1).

   If (nir > 3) call finish('Global_Coordinates','nir > 3')

   Trisection: If (nir == 3) then

      ml3 = ni/3

      !--- 1.3.7.1 Trisect the rows of the diamond.

      Do j1 = 1,2
         Do j2 = 1,2

            mi1    = (j1-1)*ni
            mi2    = j2*ml3 + 1

            zgamma = Real(j2,Double)/3.0_Double
            zchord = Sqrt(Sum((xnglob(mi1,ni+1,:,jd) - xnglob(mi1,1,:,jd))**2))
            ztheta = 2*Asin(0.5_Double*zchord)

            zbeta  = Sin(zgamma*ztheta)/Sin(ztheta)
            zalpha = Sin((1.0_Double-zgamma)*ztheta)/Sin(ztheta)

            xnglob(mi1,mi2,:,jd) = zalpha*xnglob(mi1,1   ,:,jd) + &
                                   zbeta *xnglob(mi1,ni+1,:,jd)

         End Do
      End Do

      !--- 1.3.7.2. Trisect the columns of the diamond.

      Do j1 = 1,2
         Do j2 = 1,2

           mi1    = j2*ml3
           mi2    = (j1-1)*ni + 1

           zgamma = Real(j2,Double)/3.0_Double

           zchord =  Sqrt(Sum((xnglob(ni,mi2,:,jd) - xnglob(0,mi2,:,jd))**2))
           ztheta = 2*Asin(0.5_Double*zchord)

           zbeta  = Sin(zgamma*ztheta)/Sin(ztheta)
           zalpha = Sin((1.0_Double-zgamma)*ztheta)/Sin(ztheta)

           xnglob(mi1,mi2,:,jd) = zalpha*xnglob(0 ,mi2,:,jd) + &
                                  zbeta *xnglob(ni,mi2,:,jd)

         End Do
      End Do

      !--- 1.3.7.3. Trisect the diagonal of the diamond.

      Do j2 = 1,2

         mi1 = ni - j2*ml3
         mi2 =   1 + j2*ml3

         zgamma = Real(j2,Double)/3.0_Double

         zchord = Sqrt(Sum((xnglob(0,ni+1,:,jd) - xnglob(ni,1,:,jd))**2))
         ztheta = 2*Asin(0.5_Double*zchord)

         zbeta  = Sin(zgamma*ztheta)/Sin(ztheta)
         zalpha = Sin((1.0_Double-zgamma)*ztheta)/Sin(ztheta)

         xnglob(mi1,mi2,:,jd) = zalpha*xnglob(ni,1    ,:,jd) +  &
                                zbeta *xnglob(0  ,ni+1,:,jd)

      End Do

      !--- 1.3.7.4. Compute coordinates of icosahedral triangle centers.

      Call Triangle_Center &
         (jd)     ! <-- Diamond number

   End If Trisection

   !--- 1.3.8. Bisections

   Bisection: Do jb = 0, ni2-1

      mm  = (nir)*(2**jb)
      ml  = ni/mm
      ml2 = ml/2

      !--- 1.3.8.1. Compute the rows of the diamond.

      Do j1 = 1,mm+1
         Do j2 = 1,mm

            mi1    = (j1-1)*ml
            mi2    = (j2-1)*ml + ml2 + 1
            zgamma = 0.5_Double

            zchord = Sqrt(Sum((xnglob(mi1,mi2+ml2,:,jd) - &
                               xnglob(mi1,mi2-ml2,:,jd))**2))
            ztheta = 2*Asin(0.5_Double*zchord)

            zbeta  = Sin(zgamma*ztheta)/Sin(ztheta)
            zalpha = Sin((1.0_Double-zgamma)*ztheta)/Sin(ztheta)

            xnglob(mi1,mi2,:,jd) = zalpha*xnglob(mi1,mi2-ml2,:,jd) +  &
                                   zbeta *xnglob(mi1,mi2+ml2,:,jd)

         End Do
      End Do

      !--- 1.3.8.2. Compute the columns of diamond.

      Do j1 = 1,mm+1
         Do j2 = 1,mm

            mi1 = (j2-1)*ml + ml2
            mi2 = (j1-1)*ml + 1
            zgamma = 0.5_Double

            zchord = Sqrt(Sum((xnglob(mi1+ml2,mi2,:,jd) - &
                               xnglob(mi1-ml2,mi2,:,jd))**2))
            ztheta = 2*Asin(0.5_Double*zchord)

            zbeta  = Sin(zgamma*ztheta)/Sin(ztheta)
            zalpha = Sin((1.0_Double-zgamma)*ztheta)/Sin(ztheta)

            xnglob(mi1,mi2,:,jd) = zalpha*xnglob(mi1-ml2,mi2,:,jd) +  &
                                   zbeta *xnglob(mi1+ml2,mi2,:,jd)

         End Do
      End Do

      !--- 1.3.8.3. Compute the diagonals of the diamond.

      Do j1 = 1,mm
         Do j2 = 1,mm

            mi1    = (j1-1)*ml + ml2
            mi2    = (j2-1)*ml + ml2 + 1
            zgamma = 0.5_Double

            zchord = Sqrt(Sum((xnglob(mi1+ml2,mi2-ml2,:,jd) -       &
                               xnglob(mi1-ml2,mi2+ml2,:,jd))**2))
            ztheta = 2*Asin(0.5_Double*zchord)

            zbeta  = Sin(zgamma*ztheta)/Sin(ztheta)
            zalpha = Sin((1.0_Double-zgamma)*ztheta)/Sin(ztheta)

            xnglob(mi1,mi2,:,jd) =  &
                zalpha*xnglob(mi1-ml2,mi2+ml2,:,jd) +  &
                zbeta *xnglob(mi1+ml2,mi2-ml2,:,jd)

         End Do
      End Do

   End Do Bisection


   !--- 1.3.9. Correction for round-off errors

   Where (ABS(xnglob(:,:,:,jd)) < 2.5e-14_Double)
      xnglob(:,:,:,jd) = 0.0_Double
   End Where


End Do Diamonds


!--- 1.4. Extend nodal array by two rows/colums around

xn(ng1s:ng1e,ng2s:ng2e,:,:) = xnglob(ng1s:ng1e,ng2s:ng2e,:,:)

Call Set_Boundaries &
  (xn,              & ! <-> Array to be enlarged on boundaries
   ng1sm2, ng1ep2,  & ! <-- First grid index limits
   ng2sm2, ng2ep2,  & ! <-- Second grid index limits
   1, 3,            & ! <-- Component inex limits
   2)                 ! <-- Number of additional rows


!--- 1.5. Compute longitude and latitude of core gridpoints.

Do jd = 1, nd

   Do j2 = ng2s, ng2e
      Do j1 = ng1s, ng1e
         rlon (j1,j2,jd) = Atan2(xn(j1,j2,2,jd), &
                                 xn(j1,j2,1,jd) + 1.e-20_Double)
         rlat (j1,j2,jd) = Asin(xn(j1,j2,3,jd))
      End Do
   End Do

End Do


End Subroutine Global_Coordinates



!==========================================================
Subroutine Icosahedral_Grid  &
  (kni)       ! <-- Number of triangles per diamond side
!
! Generating icosahedral triangular grid.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 29 Nov 2001 | Adaptation for use
!         |             | with ECHAM library.
!   2.0   | 08 Dec 2001 | Joint with Setup_Grid.
!   2.1   | 13 Apr 2002 | Deallocation of old data.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(in) :: kni ! Number of triangles per diamond side
!----------------------------------------------------------
! Local Scalars:
!
Integer :: ErrStat = 0    ! Error flag
!----------------------------------------------------------
! Global variables used:
!
!   ng1s,   ng1e       ! First grid index in diamond [0:ni]
!   ng1sm1, ng1ep1     ! 1-enlarged index interval [-1:ni+1]
!   ng1sm2, ng1ep2     ! 2-enlarged index interval [-2:ni+2]
!   ng2s,   ng2e       ! Second grid index in diamond [1:ni+1]
!   ng2sm1, ng2ep1     ! 1-enlarged index interval [ 0:ni+2]
!   ng2sm2, ng2ep2     ! 2-enlarged index interval [-1:ni+3]
!   ngg1s,  ngg1e      !
!   ngg2s,  ngg2e      !
!   xnglob(:,:,:,:),    & ! Cartesian  coordinates of gridpoints on unit-sphere.
!   rlon(:,:,:),        & ! Gridpoint longitudes [rad]
!   rlat(:,:,:),        & ! Gridpoint latitudes [rad]
!   xn(:,:,:,:),        & ! Local vertical unit vector = node vector
!----------------------------------------------------------


!----------------------------------------------------------
! 1. DIMENSION FACTORIZATION
!----------------------------------------------------------

ni = kni

Call Factorize_ni (ErrStat)

If (ErrStat /= 0) then
   Write (0,*) 'setup_grid: ni = ', ni, ' cannot be factorized correctly!'
   Stop
End If


!----------------------------------------------------------
! 2. SET BOUNDARY INDICES
!----------------------------------------------------------

ng1s = 0
ng1e = kni

ng2s = 1
ng2e = kni+1

!ng1sm1 = ng1s - 1
ng1sm2 = ng1s - 2
!ng1ep1 = ng1e + 1
ng1ep2 = ng1e + 2

!ng2sm1 = ng2s - 1
 ng2sm2 = ng2s - 2
!ng2ep1 = ng2e + 1
 ng2ep2 = ng2e + 2

!ngg1s = 0
!ngg1e = kni

!ngg2s = 1
!ngg2e = kni+1


!----------------------------------------------------------
! 3. ALLOCATE GRID INFORMATION FIELDS
!----------------------------------------------------------

Deallocate(xnglob, Stat = ErrStat)
Deallocate(xn,     Stat = ErrStat)
Deallocate(rlon,   Stat = ErrStat)
Deallocate(rlat,   Stat = ErrStat)

Allocate(xnglob (0     :ni    , 1     :ni+1  , 3, nd))
Allocate(xn     (ng1sm2:ng1ep2, ng2sm2:ng2ep2, 3, nd))
Allocate(rlon   (ng1s  :ng1e  , ng2s  :ng2e  ,    nd))
Allocate(rlat   (ng1s  :ng1e  , ng2s  :ng2e  ,    nd))


!----------------------------------------------------------
! 4. CALCULATE GLOBAL GRID
!----------------------------------------------------------

Call Global_Coordinates (ErrStat)

If (ErrStat /= 0) then
   Write (0,*) '  Error in Subroutine global_coordinates!'
   Stop
End If


End Subroutine Icosahedral_Grid



!==========================================================
Subroutine Multiple_Points &
  (IDX,       & ! <-- Grid point indices
   NMP,       & ! --> Number of multiple points
   MIDX)        ! --> Multiple point indices
!
! Finiding multiple grid points
!----------------------------------------------------------
! (C) Copyright 2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 14 Apr 2002 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   IDX(3)      ! Grid point indices
!
! Output arguments:
!
Integer, Intent(Out) :: &
   NMP        ! Number of multiple points
!
Integer, Intent(Out) :: &
   MIDX(5,3)  ! Multiple point grid indices
!----------------------------------------------------------
! Local Parameters:
!
Integer, Parameter :: &
   NDI(0:6)  = (/  5, 1, 2, 3, 4,  5, 1 /), &
   SDI(5:11) = (/ 10, 6, 7, 8, 9, 10, 6 /)
!
! Local Scalars:
!
Integer :: IMP  ! Multiple point index
!
! Local Arrays:
!
Integer :: W(3) ! Work array for exchange
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------

NMP       = 1
MIDX(1,:) = IDX(:)


!----------------------------------------------------------
! 1. MULTIPLE POINT PROCESSING
!----------------------------------------------------------


!--- 1.1. North hemisphere

If (IDX(3) <= 5) then

   !--- 1.1.1. North pole

   If (IDX(1) == 0 .and. IDX(2) == 1) then

      NMP       = 5
      MIDX(:,1) = 0
      MIDX(:,2) = 1
      MIDX(:,3) = (/ 1, 2, 3, 4, 5 /)

   !--- 1.1.2. Counter-polar point

   Else If (IDX(1) == NI .and. IDX(2) == NI + 1) then

      NMP       = 3
      MIDX(2,:) = (/  0, NI + 1, SDI(IDX(3) + 4) /)
      MIDX(3,:) = (/ NI,      1, SDI(IDX(3) + 5) /)

   !--- 1.1.3. Left angle

   Else If (IDX(1) == 0 .and. IDX(2) == NI + 1) then

      NMP       = 3
      MIDX(2,:) = (/ NI, NI + 1, SDI(IDX(3) + 5) /)
      MIDX(3,:) = (/ NI,      1, NDI(IDX(3) + 1) /)

   !--- 1.1.4. Right angle

   Else If (IDX(1) == NI .and. IDX(2) == 1) then

      NMP       = 3
      MIDX(2,:) = (/ NI, NI + 1, SDI(IDX(3) + 4) /)
      MIDX(3,:) = (/  0, NI + 1, NDI(IDX(3) - 1) /)

   !--- 1.1.5. Polar right side

   Else If (IDX(2) == 1) then
      NMP       = 2
      MIDX(2,:) = (/ 0, IDX(1) + 1, NDI(IDX(3) - 1) /)

   !--- 1.1.6. Polar left side

   Else If (IDX(1) == 0) then
      NMP       = 2
      MIDX(2,:) = (/ IDX(2) - 1, 1, NDI(IDX(3) + 1) /)

   !--- 1.1.7. Counter-polar right side

   Else If (IDX(1) == NI) then
      NMP       = 2
      MIDX(2,:) = (/ NI + 1 - IDX(2), NI + 1, SDI(IDX(3) + 4) /)

   !--- 1.1.8. Counter-polar left side

   Else If (IDX(2) == NI + 1) then
      NMP       = 2
      MIDX(2,:) = (/ NI, NI + 1 - IDX(1), SDI(IDX(3) + 5) /)

   End If


!--- 1.2. South hemisphere

Else

   !--- 1.2.1. South pole

   If (IDX(1) == 0 .and. IDX(2) == 1) then

      NMP       = 5
      MIDX(:,1) = 0
      MIDX(:,2) = 1
      MIDX(:,3) = (/ 6, 7, 8, 9, 10 /)

   !--- 1.2.2. Counter-polar point

   Else If (IDX(1) == NI .and. IDX(2) == NI + 1) then

      NMP       = 3
      MIDX(2,:) = (/  0, NI + 1, NDI(IDX(3) - 5) /)
      MIDX(3,:) = (/ NI,      1, NDI(IDX(3) - 4) /)

   !--- 1.2.3. Left angle

   Else If (IDX(1) == NI .and. IDX(2) == 1) then

      NMP       = 3
      MIDX(2,:) = (/ NI, NI + 1, NDI(IDX(3) - 5) /)
      MIDX(3,:) = (/  0, NI + 1, SDI(IDX(3) - 1) /)

   !--- 1.2.4. Right angle

   Else If (IDX(1) == 0 .and. IDX(2) == NI + 1) then

      NMP       = 3
      MIDX(2,:) = (/ NI, NI + 1, NDI(IDX(3) - 4) /)
      MIDX(3,:) = (/ NI,      1, SDI(IDX(3) + 1) /)

   !--- 1.2.5. Polar left side

   Else If (IDX(2) == 1) then
      NMP       = 2
      MIDX(2,:) = (/ 0, IDX(1) + 1, SDI(IDX(3) - 1) /)

   !--- 1.2.6. Polar right side

   Else If (IDX(1) == 0) then
      NMP       = 2
      MIDX(2,:) = (/ IDX(2) - 1, 1, SDI(IDX(3) + 1) /)

   !--- 1.2.7. Counter-polar left side

   Else If (IDX(1) == NI) then
      NMP       = 2
      MIDX(2,:) = (/ NI + 1 - IDX(2), NI + 1, NDI(IDX(3) - 5) /)

   !--- 1.2.8. Counter-polar right side

   Else If (IDX(2) == NI + 1) then
      NMP       = 2
      MIDX(2,:) = (/ NI, NI + 1 - IDX(1), NDI(IDX(3) - 4) /)

   End If


End If


!----------------------------------------------------------
! 3. REORDERING
!----------------------------------------------------------

IMP         = Sum(MinLoc(MIDX(1:NMP,3)))
W(:)        = MIDX(IMP,:)
MIDX(IMP,:) = MIDX(1,:)
MIDX(1,:)   = W(:)


End Subroutine Multiple_Points



End Module ICO_grid
