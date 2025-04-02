!
!+ Interpolation Routines gathered from (Michael Gorbunovs) raytracing-program
!
MODULE mo_grid_intpol
!
! Description:
!   Interpolation Routines gathered from different modules of the
!   GPS Radio-occultation raytracer-program
!   (Author: Michael Gorbunov) and adapted to the 3DVAR/PSAS grid
!   data type.
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
! V1_4         2009/03/26 Hendrik Reich
!  changes for rotated grids
! V1_5         2009/05/25 Harald Anlauf
!  Disable ftrace-regions
! V1_8         2009/12/09 Harald Anlauf
!  Symplectic_Weights_Small: improve vectorization for SX-9 by expanding SUM()
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  adjust warning messages; increase tolerance for triangle search
! V1_22        2013-02-13 Andreas Rhodin
!  changed interface (pass ngp ! number of neighbour grid-points)
! V1_23        2013-03-26 Harald Anlauf
!  Add horizontal interpolation for ICON grid
! V1_26        2013/06/27 Harald Anlauf
!  Add cases for ICON hor.grid
! V1_27        2013-11-08 Harald Anlauf
!  GPSRO: horizontal interpolation for ICON
! V1_28        2014/02/26 Andreas Rhodin
!  move subroutine fix_grid_indices from mo_letkf to mo_grid_intpol
! V1_42        2015-06-08 Andreas Rhodin
!  ICON LETKF coarse grid; MEC temporal interpolation; GPSRO n.n. interpolation
! V1_43        2015-08-19 Andreas Rhodin
!  refine nearest neighbour interpolation for flake variables
! V1_45        2015-12-15 Andreas Rhodin
!  parallelise soil/surface interpolation
! V1_47        2016-06-06 Andreas Rhodin
!  subroutine hor_intp: allow for MPI partitioned output filed
! V1_50        2017-01-09 Andreas Rhodin
!  adapt COSMO-MEC to ICON-LAM: return the 6 surrounding grid-points
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Michael Gorbunov                      original code
! Andreas Rhodin    DWD/MPI  2002-2007  adapted to 3DVAR/PSAS, extended
! Harald Anlauf     DWD      2007-2008  optimizations for SX8
!============================================================================

! Compiler options: (Hopefully, this will work in the future)
!!NEC$ options "-finline-file=../../../../occ_lib/interpolation.f90"

!
! Diagnostics for vectorization
!
#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION) && 0
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!============================================================================
  !--------------
  ! Modules used:
  !--------------
  Use mo_kind,          only: wp, i8
  Use mo_exception,     only: finish, message
  Use mo_wmo_tables,    only: WMO6_GAUSSIAN,   &!
                              WMO6_LATLON,     &!
                              WMO6_ROTLL,      &!
                              DWD6_ICOSAHEDRON,&!
                              DWD6_ICON         !
  use mo_mpi_dace,      only: dace,            &!
                              p_sum,           &!
                              p_max,           &!
                              p_allgather       !
  use mo_run_params,    only: nproc1, nproc2
  !----------------------------------------------------------------------------
  Use mo_physics,       only: dtr => d2r
  Use mo_physics,       only: r2d, d2r, rearth, gacc
! Use ECHAM_fields,     only: NG1
  Use mo_ico_grid,      only: nspoke2
  Use mo_t_obs,         only: t_obs, t_mcols, t_mcol, t_hic, t_coord, t_icol, &
                              t_imcol, set_xuv
  Use mo_atm_grid,      only: t_grid, phi2phirot, rla2rlarot, &
                              same_horizontal_grid ! check for same grid
  Use mo_icon_grid,     only: search_icon_global
  use mo_dace_string,   only: tolower,         &! string to lower case
                              split             ! string (words) to char-array
  Implicit None
  !----------------
  ! Public entities
  !----------------
  private
  public :: Grid_Indices        ! get indices of neighbor gridpoints
  public :: idx_init            ! determine model indices for given location
  public :: ECHAM_rays_idx_init ! determine model indices for a ray
  public :: ECHAM_1d_idx_init   ! determine model indices for 1d operator
  public :: add_index           ! utility to set  indices
  public :: alloc_imcol         ! (re)allocate pointer argument imcol
  public :: hor_intp            ! interpolate field horizontally
  public :: hor_intp_coeff      ! calculate coefficients for 'hor_intp'
  public :: t_h_idx             ! derived type for indices,weights
  public :: mp                  ! max.no. of points for hor. interpolation
! public :: NG1                 ! No. grid points for lat/lon interpolation
  public :: fix_grid_indices    ! fix interpolation coefficients to GME grid
  public :: distribute_icon_gme ! distribute ICON grid for GME partitioning
  public :: HORINT_DEFAULT      ! interpolate between multiple gp
  public :: HORINT_SURROUND     ! interpolate between surrounding gp
  public :: HORINT_NEAREST      ! select horizontally nearest gp
  public :: HORINT_ZDIST        ! select smallest vertical distance
  public :: HORINT_FRLAND       ! select largest land fraction
  public :: HORINT_FRSEA        ! select largest sea  fraction
  public :: HORINT_LAND         ! prefer landpoints if available
  public :: HORINT_SEA          ! prefer seapoints  if available
  public :: HORINT_LANDONLY     ! allow  only landpoints
  public :: HORINT_SEAONLY      ! allow  only seapoints
  public :: text2mode           ! decode   hor.interpolation mode specification
  public :: mode2text           ! describe hor.interpolation mode specification
  public :: check_horint        ! check    hor.interpolation mode specification
  public :: filter_idx          ! select gridpoints for hor.interpolation mode
  !---------
  ! datatype
  !---------
  integer ,parameter :: mp = 4  ! max.no.points (linear interp. currently)
  type t_h_idx                  ! derived type for indices,weights
    real(wp) :: w    (mp)       ! weights
    integer  :: ijdp (mp,4)     ! indices
    integer  :: n               ! number of points used
  end type t_h_idx

  !-------------------------------------------------------------
  ! Steering of gridpoint selection for horizontal interpolation
  !-------------------------------------------------------------
  integer, parameter :: HORINT_DEFAULT  = 2**0  ! interpolate between multiple gp
  integer, parameter :: HORINT_SURROUND = 2**1  ! interpolate between surround.gp
  integer, parameter :: HORINT_NEAREST  = 2**2  ! select horizontally nearest gp
  integer, parameter :: HORINT_ZDIST    = 2**3  !   ~ smallest vertical distance
  integer, parameter :: HORINT_FRLAND   = 2**4  !   ~ largest land fraction
  integer, parameter :: HORINT_FRSEA    = 2**5  !   ~ largest sea  fraction
  integer, parameter :: HORINT_LAND     = 2**6  ! prefer landpoints if available
  integer, parameter :: HORINT_SEA      = 2**7  ! prefer seapoints  if available
  integer, parameter :: HORINT_LANDONLY = 2**8  ! allow  only landpoints
  integer, parameter :: HORINT_SEAONLY  = 2**9  ! allow  only seapoints
  integer, parameter :: HORINT_CUSTOM   = 2**10 ! custom method (reserved)
  integer, parameter :: HORINT_OTHER    = 2**11 ! other / not yet implemented

  type t_horint
     integer       :: mode      ! bitmask
     character(12) :: desc      ! description of gridpoint selection
  end type t_horint

#ifndef __PGI
  type(t_horint), parameter :: horint_modes(*) = &
#else   /* Work around implied-shape bug in nvidia 23.11 */
  type(t_horint), parameter :: horint_modes(13) = &
#endif
     [ t_horint (HORINT_DEFAULT,  "default")    ,&
       t_horint (HORINT_SURROUND, "surrounding"),&
       t_horint (HORINT_NEAREST,  "nearest")    ,&
       t_horint (HORINT_ZDIST,    "zdist")      ,&
       t_horint (HORINT_FRLAND,   "frland")     ,&
       t_horint (HORINT_FRSEA,    "frsea")      ,&
       t_horint (HORINT_LAND,     "preferland") ,&
       t_horint (HORINT_SEA,      "prefersea")  ,&
       t_horint (HORINT_LAND,     "land")       ,&! same as preferland
       t_horint (HORINT_SEA,      "sea")        ,&! same as prefersea
       t_horint (HORINT_LANDONLY, "landonly")   ,&
       t_horint (HORINT_SEAONLY,  "seaonly")    ,&
       t_horint (HORINT_CUSTOM,   "RESERVED")    ]
  !-----------
  ! Interfaces
  !-----------
  interface Grid_Indices
    module procedure Grid_Indices_geo  ! grid indices from geodetic lon,lat
    module procedure Grid_Indices_cart ! grid indices from cartesian coord.
  end interface Grid_Indices

  interface hor_intp
    module procedure hor_intp_field
  end interface hor_intp
  !----------
  ! constants
  !----------
  integer, parameter :: mimcol = 200   ! estimated max.number of columns
  !
  ! Number of grid points for lat/lon interpolation (formerly ECHAM_fields):
  Integer, Parameter :: NG1    = 2     ! 2 - linear interpolation
Contains
!==============================================================================

Subroutine idx_init (cor, hic, mc, iatm, itrac, grid, ts, tw, order)
!---------------------------------------------------------------------
! Determine model grid indices for a given location given by 'cor'.
! The horizontal interpolation coefficients 'hic' and the list
! of required model columns 'mc' (actual arguments usually passed from
! 'spot% col% h' and obs% o% mc) and are updated accordingly.
!---------------------------------------------------------------------
type(t_coord) ,intent(in)    :: cor     ! coordinates
type(t_hic)   ,intent(inout) :: hic     ! horizontal interpolation coefficients
type(t_mcols) ,intent(inout) :: mc      ! model column descriptor
integer(i8)   ,intent(in)    :: iatm    ! ids of atmospheric parameters required
integer       ,intent(in)    :: itrac   ! tracers required
type(t_grid)  ,intent(in)    :: grid    ! model grid description
integer       ,intent(in)    :: ts      ! time slot
real(wp)      ,intent(in)    :: tw      ! time interpolation weight
integer       ,intent(in)    :: order   ! 1:nn 2:linear 4:higher order
               optional      :: order

  integer               :: idx (16,4)   !(point, ijdp)
  integer               :: i, it, ix, inn
  integer               :: tsi
  integer               :: np
  integer               :: iw12
  type(t_mcol) ,pointer :: c, tmp(:)
  type(t_hic)  ,save    :: hic0
  logical               :: nn

  iw12 = hic% iw12
  hic  = hic0
  nn   = .false.; if (present (order)) nn = (order==1)
  if (nn .and. iw12 == 1) call finish ('idx_init','nn .and. iw12 == 1')

  !----------------------------
  ! determine neighbour columns
  !----------------------------
  select case (iw12)
  case default
    call Grid_Indices_geo &
      (        cor% dlon,& ! <-- geodetic longitude
               cor% dlat,& ! <-- geodetic latitude
               grid,     & ! <-- grid data type
               IDX,      & ! --> Grid point indices [Point, index]
               hic% w,   & ! --> Weight
               np,       & ! --> number of points returned
        order= order     ) ! <~~ order of interpolation
  case (1)
    call Grid_Indices_geo &
      (        cor% dlon,& ! <-- geodetic longitude
               cor% dlat,& ! <-- geodetic latitude
               grid,     & ! <-- grid data type
               IDX,      & ! --> Grid point indices [Point, index]
               hic% w,   & ! --> Weight
               np,       & ! --> number of points returned
               hic% w1,  & ! ~~> gradient: d weight / d x
               hic% w2,  & ! ~~> gradient: d weight / d y
        order= order     ) ! <~~ order of interpolation
  end select

  !----------------------------
  ! just keep nearest neighbour
  !----------------------------
  inn = 1
  if (np > 0) then
    inn = sum(maxloc(hic% w(1:np)))
    if (nn) then
      idx    (1,:) = idx    (inn,:)
      hic% w (1)   = 1._wp
      inn          = 1
      np           = 1
    endif
  end if

  !-------------------------------------------------
  ! copy indices and weights to subroutine arguments
  !-------------------------------------------------
  hic% imc (:,:)   = 0
  hic% w   (np+1:) = 0._wp

  tsi = ts
  do it = 1, 2
    if (it == 2) then
      if (tw == 0._wp) exit
      tsi = ts + 1
    endif
    do i=1,np
      ix = mc% idx (idx(i,1),idx(i,2),idx(i,3),tsi)
      if (ix == 0) then
        mc% n = mc% n + 1
        ix = mc% n
        mc% idx (idx(i,1),idx(i,2),idx(i,3),tsi) = ix
        if (size (mc% c) < ix) then
          tmp => mc% c
          allocate (mc% c (int(1.5 * ix)))
          mc% c (1:size(tmp)) = tmp
          deallocate (tmp)
        endif
        c => mc% c(ix)
        c% icol  = ix
        c% iatm  = iatm
        c% itrac = itrac
        c% ijdtp = (/idx(i,:3),tsi,idx(i,4)/)
      else
        c => mc% c(ix)
        c% iatm  = ior (iatm,  c% iatm)
        c% itrac = ior (itrac, c% itrac)
      endif

      hic% imc (i,it) = ix                                  ! <<<<<<<

    end do
  end do

  if (np > 0) then
    hic% ijdp(1:3) = mc% c (hic% imc (inn,1))% ijdtp(1:3)
    hic% ijdp(4)   = mc% c (hic% imc (inn,1))% ijdtp(  5)
  end if

end Subroutine idx_init
!==============================================================================
Subroutine ECHAM_rays_idx_init &
  (XT,        & ! <-- Transmitter position
   XR,        & ! <-- Receiver position
   XLC,       & ! <-- Curvature center
   grid,      & ! <-- grid description variable
   iatm,      & ! <-- parameters required (T+Q)
   natm,      & ! <-- number of parameters required (2)
   ts,        & ! <-- time slot
   obs,       & ! <-> Observation data type
   imcol)       ! --> List of model columns required

!
! Initialization of dynamical data indices.
! Determine model indices required along a ray of a
! radio occultation.
!----------------------------------------------------------
! Method:
!   Tracing the surface projection of
!   XT-XR ray and marking necessary meshes.
!----------------------------------------------------------
! (C) Copyright 2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Apr 2002 | Original code.
!----------------------------------------------------------
  ! Modules used:
  !--------------
  Use Coordinates, only: Cartesian,     &! Imported Type Definitions:
                         Vector_Normed, &! Imported Routines:
                         Vector_Norm,   &
                         Operator(*),   &! Imported Operators:
                         Operator(+),   &
                         Operator(-),   &
                         Operator(.x.), &
                         Operator(/)
  Use Earth,       only: R_Earth,       &! Imported Parameters:
                         H_atm
  Implicit None
  !-----------
  ! arguments:
  !-----------
  Type(Cartesian), Intent(In)    :: XT       ! Transmitter position
  Type(Cartesian), Intent(In)    :: XR       ! Receiver position
  Type(Cartesian), Intent(In)    :: XLC      ! Curvature center
  type(t_grid),    Intent(In)    :: grid     ! grid description variable
  integer(i8),     intent(in)    :: iatm     ! parameters required
  integer,         intent(in)    :: natm     ! number of atm. parameters req.
  integer,         intent(in)    :: ts       ! time slot
  type(t_obs),     intent(inout) :: obs      ! observations data type
  type(t_imcol),   pointer       :: imcol(:) ! list of model columns required

  !----------------------------------------------------------
  ! Local Parameters:
  !
  Real(wp), Parameter :: &
     R_atm = R_Earth + H_atm + 60.0_wp  ! Effective radius of atmosphere
  !
  Real(wp), Parameter :: &
     DX = 15.0_wp                       ! Discretization step
  !
  ! Local Scalars:
  !
  Type(Cartesian) :: XTC      ! Transmitter coordinates relative to XLC
  Type(Cartesian) :: XRC      ! Receiver coordinates relative to XLC
  Type(Cartesian) :: XN       ! Normal to occultation plane
  Type(Cartesian) :: X        ! Current point of main ray
  Type(Cartesian) :: XS       ! Side point
  Real(wp)        :: D        ! Discriminant of square equation
  Real(wp)        :: S1       ! Beginning of ray part inside atmosphere
  Real(wp)        :: S2       ! End of ray part inside atmosphere
  Type(Cartesian) :: X1       ! Beginning of ray part inside atmosphere
  Type(Cartesian) :: X2       ! End of ray part inside atmosphere
  Integer         :: NS       ! Number of discrete points
  Integer         :: IS       ! Index of discrete points
  Integer         :: JS       ! Index of side points
  Integer         :: NGP      ! Subgrid dimension
  Integer         :: np       ! Subgrid dimension (returned from Grid_Indices)
  Integer         :: nimcol   ! actual number of columns required
  !
  ! New (A.Rhodin)
  !
  Integer  :: Stat            ! Error status
  real(wp) :: W(9)            ! weights (horizontal interpolation)
  integer  :: ixmcol(9)       ! indices to imcol
  !
  ! Local Arrays:
  !
  Integer, Allocatable ::  &
     IDX(:,:)              ! Subgrid indices [point, index]
  !----------------------------------------------------------


  !----------------------------------------------------------
  ! 1. INITIALIZATION
  !----------------------------------------------------------


  !--- 1.1. Relative coordinates

  XTC = XT - XLC
  XRC = XR - XLC


  !--- 1.2. Normal to occultation plane

  XN  = Vector_Normed(XTC .x. XRC)


  !--- 1.3. Interpolation subgrid size

  Select Case (grid% gridtype)
     Case (WMO6_GAUSSIAN, WMO6_LATLON, WMO6_ROTLL)
        NGP = NG1**2
     Case (DWD6_ICOSAHEDRON, DWD6_ICON)
        NGP = 3
  End Select

  Allocate(IDX(NGP,4))


  !--- 1.4. Index field initialization

  Stat          = 0

  !----------------------------------------------------------
  ! 2. DETERMINATION OF RAY PART INSIDE ATMOSPHERE
  !----------------------------------------------------------

  !--- 2.1. Discriminant

  D = (XRC*XTC)**2 - (XRC*XRC)*(XTC*XTC) + &
       R_atm**2*((XTC-XRC)*(XTC-XRC))

  If (D <= 0) then
     Stat = 1
     call finish ('ECHAM_rays_idx_init','Discriminant <= 0')
     Return
  End If


  !--- 2.2. Roots

  S1 = (- XRC*(XTC-XRC) + Sqrt(D))/((XTC-XRC)*(XTC-XRC))
  S2 = (- XRC*(XTC-XRC) - Sqrt(D))/((XTC-XRC)*(XTC-XRC))

  X1 = S1*XTC + (1-S1)*XRC
  X2 = S2*XTC + (1-S2)*XRC


  !--- 2.3. Determination of discrete step

  NS = Ceiling(Vector_Norm(X2-X1)/DX)


  !------------------------------------------------------------------------
  ! determine number of gridpoints already associated with this occultation
  !------------------------------------------------------------------------
  if (associated(imcol)) then
    nimcol = size(imcol)
  else
    nimcol = 0
    call alloc_imcol (imcol, mimcol)
  endif

  !----------------------------------------------------------
  ! 3. SCANNING RAY
  !----------------------------------------------------------

  Do Is = -20, NS+20
     X = ((NS-IS)*X1 + IS*X2)/real(NS,wp)
     Do JS = -2, 2
        XS = X + JS*DX*XN
        Call Grid_Indices &
          (XS,       & ! <-- Cartesian coordinates of point
           grid,     & ! <-- grid data type
           IDX,      & ! --> Grid point indices [Point, index]
           W,        & ! --> Interpolation weights (not used here)
           np)         ! --> number of points returned
        if(np==0) call finish ('ECHAM_rays_idx_init',&
                               'cannot handle out of domain condition')
        !------------------------------------------------
        ! loop over grid indices associated with this ray
        !------------------------------------------------
        call add_index (idx, np, obs, iatm, natm, imcol, nimcol, ixmcol, grid,&
                         ts, 0._wp) ! +++ currently no multiple time slices
     end do
  end do

  !----------------------------------------------------------
  ! 4. MEMORY DEALLOCATION
  !----------------------------------------------------------

  call alloc_imcol (imcol, nimcol)
  Deallocate(IDX)

end subroutine ECHAM_rays_idx_init
!------------------------------------------------------------------------------

subroutine add_index (idx, ngp, obs, iatm, natm, imcol, nimcol, ixmcol, grid, &
                      ts, ws                                                  )
!--------------------------------------------------------------------------
! Insert given model column indices (specified by 'idx', 'ngp', 'ts', 'ws')
! in the data structures describing the model columns required by the
! observation operators (for operators which work directly on the model
! columns, i.e. GPSRO, STD/ZTD, and COSMO observation operstors.
! The fields 'obs% mc% ..' and 'imcol' are updated accordingly.
!--------------------------------------------------------------------------
integer      ,intent(in)    :: IDX(:,:) ! Subgrid indices [point, index]
integer      ,intent(in)    :: ngp      ! number of neighbour points
type(t_obs)  ,intent(inout) :: obs      ! observations data type
integer(i8)  ,intent(in)    :: iatm     ! parameters required
integer      ,intent(in)    :: natm     ! number of atm. parameters req.
type(t_imcol),pointer       :: imcol(:) ! list of model columns required
Integer      ,intent(inout) :: nimcol   ! actual number of columns required
Integer      ,intent(out)   :: ixmcol(:)! indices to imcol
type(t_grid) ,Intent(in)    :: grid     ! grid description variable
integer      ,intent(in)    :: ts       ! time slots
real(wp)     ,intent(in)    :: ws       ! time interpolation weight

  Integer               :: IGP      ! Subgrid index
  Integer               :: I1, I2   ! Grid indices
  Integer               :: ID       ! Diamond index
  Integer               :: ix       ! index variable
  type(t_mcol) ,pointer :: c, tmp(:)
  logical               :: new      ! flag for new entry in index array imcol
  integer               :: i        ! index variable
  integer               :: it       ! time slot index
  integer               :: tsi

  ixmcol = 0

  !----------------------
  ! loop over grid points
  !----------------------
  Do IGP = 1, NGP
    I1 = IDX(IGP,1)
    I2 = IDX(IGP,2)
    ID = IDX(IGP,3)

    !---------------------
    ! loop over time slots
    !---------------------
    tsi = ts
    do it = 1, 2
      if (it == 2) then
        if (ws == 0._wp) exit
        tsi = ts + 1
      endif

      !--------------------------------------------------------------
      ! grid point/time slot not yet associated with any occultation:
      !--------------------------------------------------------------
      If (obs% mc% idx(I1,I2,ID,TSI) == 0) then
        obs% mc% n            = obs% mc% n + 1
        ix                    = obs% mc% n
        obs% mc% idx(I1,I2,ID,TSI) = ix

        !------------------------------------------
        ! extend index array (for all occultations)
        !------------------------------------------
        if (size (obs% mc% c) < ix) then
          tmp => obs% mc% c
          allocate (obs% mc% c (int(1.5 * ix)))
          obs% mc% c (1:size(tmp)) = tmp
          deallocate (tmp)
        endif
        c => obs% mc% c(ix)
        c% icol  = ix
        c% iatm  = iatm
        c% itrac = 0
        c% ijdtp(:3) = idx(igp,:3)
        c% ijdtp( 4) = tsi
        c% ijdtp( 5) = idx(igp, 4)
        !------------------------------------------
        ! extend index array (for this occultation)
        !------------------------------------------
        call add_idx
        ixmcol(igp) = nimcol
      !---------------------------------------------------
      ! grid point already associated with an occultation:
      !---------------------------------------------------
      else
        ix = obs% mc% idx(I1,I2,ID,TSI)
        c => obs% mc% c(ix)
        c% iatm  = ior (iatm,  c% iatm)
        new = .true.
        do i=1,nimcol
          if (imcol(i)% imc(it) == ix) then
            new = .false.
            ixmcol(igp) = i
            exit
          endif
        end do
        if (new) then
          call add_idx
          ixmcol(igp) = nimcol
        endif
      End If
    end do
  End Do

contains

  !--------------------------------------
  ! add an index to the index array imcol
  !--------------------------------------
  subroutine add_idx
    if (it == 1) then
      nimcol = nimcol + 1
      if (nimcol > size(imcol)) call alloc_imcol (imcol, 2*nimcol)
      imcol (nimcol)%   iatm    = iatm
      imcol (nimcol)%   natm    = natm
      imcol (nimcol)%   itrac   = 0
      imcol (nimcol)%   nsl     = 1
      imcol (nimcol)%c% dlat    = grid% rlat (i1,i2,1,id) * r2d
      imcol (nimcol)%c% dlon    = grid% rlon (i1,i2,1,id) * r2d
      call set_xuv (imcol (nimcol))
    endif
    imcol   (nimcol)%   imc(it) = ix
  end subroutine add_idx

End Subroutine add_index
!------------------------------------------------------------------------------
!------------------------------------------------
! (re)allocate pointer argument imcol with size n
! copy content of old pointer
!------------------------------------------------
subroutine alloc_imcol (imcol, n)
type (t_imcol) ,pointer    :: imcol(:)
integer        ,intent(in) :: n
  integer                 :: m
  type (t_imcol) ,pointer :: tmp(:)
  allocate (tmp(n))
  if (associated (imcol)) then
    m = min (n,size(imcol))
    tmp (:m) = imcol(:m)
    deallocate (imcol)
  endif
  imcol => tmp
end subroutine alloc_imcol

!==============================================================================
Subroutine ECHAM_1d_idx_init (col, obs, iatm, natm, grid, ts, nneighb, imcol)
type(t_icol) ,intent(inout) :: col      ! interpolated column description
type(t_obs)  ,intent(inout) :: obs      ! observations data type
integer(i8)  ,intent(in)    :: iatm     ! atmospheric parameters required (T+Q)
integer      ,intent(in)    :: natm     ! number of atmospheric parameters
type(t_grid) ,intent(in)    :: grid     ! model grid description
integer      ,intent(in)    :: ts       ! time slot
logical      ,intent(in)    :: nneighb  ! nearest neighbour mode
type(t_imcol),pointer       :: imcol(:) ! list of model columns required

!-------------------------------------------------------
! determine data indices for 1d operator (tangent point)
!-------------------------------------------------------
  integer              :: ngp      ! subgrid size
  integer              :: igp      ! subgrid index
  integer              :: i1,i2,id,ix
  integer              :: order

  order = 2; if (nneighb) order = 1

  call idx_init (col% c, col% h, obs% mc, iatm, 0, grid, ts, 0._wp, order)

  select case (grid% gridtype)
     case (WMO6_GAUSSIAN, WMO6_LATLON, WMO6_ROTLL)
       ngp = ng1**2
     case (DWD6_ICOSAHEDRON, DWD6_ICON)
       ngp = 3
  end select

  if (nneighb) ngp = 1

  !------------------------------------------
  ! extend index array (for this occultation)
  !------------------------------------------
  allocate (imcol(ngp))
  do igp = 1, ngp
    ix = col% h%         imc (igp,1)
    i1 = obs% mc% c(ix)% ijdtp(1)
    i2 = obs% mc% c(ix)% ijdtp(2)
    id = obs% mc% c(ix)% ijdtp(3)
    imcol (igp)%   imc(1)= col% h% imc(igp,1) ! index to model column
    imcol (igp)%   iatm  = iatm               ! atmospheric parameter ids
    imcol (igp)%   natm  = natm               ! number of atmospheric parameters
    imcol (igp)%   itrac = 0                  ! bit field for tracers
    imcol (igp)%   nsl   = 1                  ! number of single level parameters
    imcol (igp)%c% dlat  = grid% rlat (i1,i2,1,id) * r2d
    imcol (igp)%c% dlon  = grid% rlon (i1,i2,1,id) * r2d
    call set_xuv (imcol (igp))
  end do

end Subroutine ECHAM_1d_idx_init
!==============================================================================
Subroutine Grid_Indices_cart &
  (X,        & ! <-- Cartesian coordinates of point
   grid,     & ! <-- grid data type
   IDX,      & ! --> Grid point indices [Point, index]
   W,        & ! --> Weights            [Point]
   np)         ! --> number of points returned
!----------------------------------------------------------
! Calculation of interpolated (1 + N)*Grad(N) for
! ECHAM global fields.
!----------------------------------------------------------
! Method:
!   Calculation of gradient in geodetic coordinates and
!   transform to Cartesian coordinates.
!----------------------------------------------------------
! (C) Copyright 2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Apr 2002 | Original code.
!         !  2 May 2002 | special routine for cartesian coord. (A.Rhodin)
!----------------------------------------------------------
  !
  ! Modules used:
  !
  Use Coordinates, only: Cartesian
  Use Earth,       only: Geodetic,                 &! Type Definition
                         Geod_from_Cart             ! Imported Routine
  Implicit None
  !-----------
  ! arguments:
  !-----------
  Type(Cartesian), Intent(In)  :: X        ! Cartesian coordinates of point
  Type(t_grid),    Intent(In)  :: grid     ! grid description variable
  Integer,         Intent(Out) :: IDX(:,:) ! Grid point indices [Point, index]
  Real(wp),        Intent(Out) :: w (:)    ! weights
  Integer,         Intent(Out) :: np       ! number of points returned
!----------------------------------------------------------
! Local Scalars:
!
Type(Geodetic) :: G       ! Geodetic coordinates of point
Real(wp)       :: PLon    ! Longitude of point
Real(wp)       :: PLat    ! Latitude of point
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------


!--- 1.1. Transform to geodetic coordinates

G    = Geod_from_Cart(X)
PLon = G%Lambda
PLat = G%Phi

call Grid_Indices_geo &
  (plon,     & ! <-- geodetic longitude
   plat,     & ! <-- geodetic latitude
   grid,     & ! <-- grid data type
   IDX,      & ! --> Grid point indices [Point, index]
   w,        & ! --> Weights
   np)         ! --> number of points returned

end Subroutine Grid_Indices_cart
!------------------------------------------------------------------------------
Subroutine Grid_Indices_geo &
  (plon,     & ! <-- geodetic longitude
   plat,     & ! <-- geodetic latitude
   grid,     & ! <-- grid data type
   IDX,      & ! --> Grid point indices [Point, index]
   w,        & ! --> Weights
   np,       & ! --> number of points returned
   w1,       & ! ~~> gradient: d weight / d x
   w2,       & ! ~~> gradient: d weight / d y
   lunique,  & ! <~~ set unique index
   order,    & ! <~~ 1:nn, 2:linear, 4:higher order
   mode      ) ! <~~ mode for ICON grid
!----------------------------------------------------------
! Calculation of interpolated (1 + N)*Grad(N) for
! ECHAM global fields.
!----------------------------------------------------------
! Method:
!   Calculation of gradient in geodetic coordinates and
!   transform to Cartesian coordinates.
!----------------------------------------------------------
! (C) Copyright 2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Apr 2002 | Original code.
!         |  2 May 2002 | special routine for geodetic coord. (A.Rhodin)
!         | 31 Jan 2007 | optional arguments w1, w2           (A.Rhodin)
!----------------------------------------------------------
  !
  ! Modules used:
  !
  Implicit None
  !-----------
  ! arguments:
  !-----------
  Real(wp),         Intent(In)  :: PLon     ! Geodetic Longitude of point
  Real(wp),         Intent(In)  :: PLat     ! Geodetic Latitude of point
  Type(t_grid),     Intent(In)  :: grid     ! grid description variable
  Integer,          Intent(Out) :: IDX(:,:) ! Grid point indices [Point, index]
  Real(wp),         Intent(Out) :: w  (:)   ! weights
  Real(wp),Optional,Intent(Out) :: w1 (:)   ! gradient: d weight / d x
  Real(wp),Optional,Intent(Out) :: w2 (:)   ! gradient: d weight / d y
  Integer,          Intent(Out) :: np       ! number of points returned
  Logical, Optional,Intent(in)  :: lunique  ! set unique index
  Integer, Optional,Intent(in)  :: order    ! interpolation order
  Integer, Optional,Intent(in)  :: mode     ! mode for ICON grid
!----------------------------------------------------------
! Local Scalars:
!
Integer        :: IGP     ! Subgrid point index
Integer        :: KLon    ! Longitude subgrid index
Integer        :: KLat    ! Latitude subgrid index
Integer        :: m1, m2  ! Interpolation triangle indices
Logical        :: lu      ! set unique index
logical        :: nn      ! Copy of nneighb
integer        :: i
integer        :: lmode   ! Copy of mode (ICON grid only)
logical        :: l_triangle ! Use fixed triangle interpolation
!
! Local Arrays:
!
! --- Interpolation on rectangular grid
!

Integer, Parameter :: MG = 4 ! max.# of interpol. points (in 1d, for GME )
Integer, Parameter :: MM = 6 ! max.# of interpol. points (in 2d, for ICON)

real(wp),parameter :: TOL = 1.e-12_wp ! Rounding tolerance for ROT_LL

Integer :: ng           ! number of grid points (order)
Integer :: ILon(MG)     ! Longitude interpolation subgrid
Integer :: ILat(MG)     ! Latitude cell indices
Real(wp) :: &
   WLon (MG),         & ! Weight for longitudinal interpolation
   WLat (MG),         & ! Weight for latitudinal interpolation
   WLon1(MG),         & ! Gradient
   WLat1(MG),         & ! Gradient
   PLon_loc,          & ! local PLon
   PLat_loc             ! local PLat
!
! --- Interpolation on icosahedral grid
!
Real(wp) :: &
   PRLon (1),          & ! Longitude of interpolation point
   PRLat (1),          & ! Latitude  of interpolation point
   SWS (3,1)             ! Symplectic weights for triangle vertices
Integer :: &
   KIDX(4,1),          & ! Interpolation triangle index
   ngp   (1)             ! number of neighbour grid-points returned
integer :: pe_i1_i2_id(4)
!
! --- Interpolation on ICON grid
!
Integer  :: jl (MM,1), & ! Line  indices
            jb (MM,1)    ! Block indices
Real(wp) :: WS (MM,1)    ! Weights for interpolation
!----------------------------------------------------------

!----------------------------
! process optional parameters
!----------------------------
  lu = .true. ; if (present (lunique)) lu = lunique
  nn = .false.
  ng = 2
  l_triangle=.false.
  if (present (order)) then
    select case (order)
    case (1)
      nn = .true.
    case (3)
      l_triangle=.true.
    case (4)
      ng = 4
    case (5)
      ng = 4
      l_triangle=.true.
    end select
  endif

!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATION WEIGHTS
!----------------------------------------------------------


Select Case (grid% gridtype)

   !--- 2.1. Interpolation on rectangular grid

   Case (WMO6_GAUSSIAN, WMO6_LATLON, WMO6_ROTLL)

      !--- special treatment for rotated grid (WMO6_ROTLL)
      if (grid% rot) then
         PLon_loc = rla2rlarot(PLat,PLon,grid% dlatr,grid% dlonr,0._wp)
         PLat_loc = phi2phirot(PLat,PLon,grid% dlatr,grid% dlonr)
         !------------------------------------------
         ! Use rounding tolerance near grid boundary
         ! (limited accuracy of transcendentals)
         !------------------------------------------
         if (abs (PLon_loc - grid% dlon(grid% nx)) < TOL) then
            PLon_loc = grid% dlon(grid% nx)
         end if
         if (abs (PLat_loc - grid% dlat(grid% ny)) < TOL) then
            PLat_loc = grid% dlat(grid% ny)
         end if
      else
         PLon_loc = PLon
         PLat_loc = PLat
      end if

      !--- 2.1.1. Longitudinal interpolation

      Call Longitude_Interpolation &
        (grid% dlon,    & ! <-- Longitude grid [deg]
         PLon_loc,      & ! <-- Longitude of point [deg]
         grid% cyc_x,   & ! <-- cyclic boundary conditions
         ILon(1:ng),    & ! --> Interpolation subgrid
         WLon,          & ! --> Weights of subgrid points
         W1)              ! ~~> Gradient

      !--- 2.1.2. Latitudinal interpolation

      Call Latitude_Interpolation &
        (grid% dlat,    & ! <-- Latitude grid [deg]
         PLat_loc,      & ! <-- Latitude of point [deg]
         grid% poly,    & ! <-- poles at bounds
         ILat(1:ng),    & ! --> Interpolation subgrid
         WLat,          & ! --> Weights of subgrid points
         W2)              ! ~~> Gradient

      !--- Check for out of domain condition

      If (Ilat(1)==0 .or. Ilon(1)==0) then
        idx = 0
        w   = 0._wp
        np  = 0
        return
      Endif


      !--- 2.1.3. Interpolation indices and weights

      if (PLat_loc < 89.9_wp .and. PLat_loc > -89.9_wp) then
        if (present(w1)) Wlon1(1:ng) = w1 (1:ng) / cos (d2r*PLat_loc)
        if (present(w2)) Wlat1(1:ng) = w2 (1:ng)
      else
        Wlon1 = 0._wp    !+++ currently no gradients at the poles
        Wlat1 = 0._wp
      endif

      IGP = 0

      Do KLat = 1,ng
         Do KLon = 1,ng
            IGP         = IGP + 1
            IDX(IGP,1)  = ILon(KLon)
            IDX(IGP,2)  = ILat(KLat)
            IDX(IGP,3)  = 1
            w  (igp)    = wlon(KLon)*wlat(KLat)
            if(present(w1)) w1 (igp) = wlon1(KLon) * wlat (KLat)
            if(present(w2)) w2 (igp) = wlon (KLon) * wlat1(KLat)
         End Do
      End Do

      np = igp

   !--- 2.2. Interpolation on icosahedral grid

   Case (DWD6_ICOSAHEDRON)

      PRLon(1) = PLon
      PRLat(1) = PLat

FTRACE_BEGIN("Grid_Indices_geo:ico")
      Call Setup_Interpolation &
        (PRLon,      & ! <-- Point longitudes [deg]
         PRLat,      & ! <-- Point latitudes [deg]
         grid,       & ! <-- grid data type
         KIDX,       & ! --> Index array for interpolation
         SWS)          ! --> Symplectic weights for triangle vertices
FTRACE_END  ("Grid_Indices_geo:ico")

      m1         = KIDX(4,1)
      m2         = Mod(m1,6) + 1
      IDX(1,1:3) = KIDX(1:3,1)
      IDX(2,1)   = KIDX(1,1) + nspoke2(1,m1)
      IDX(2,2)   = KIDX(2,1) + nspoke2(2,m1)
      IDX(2,3)   = KIDX(3,1)
      IDX(3,1)   = KIDX(1,1) + nspoke2(1,m2)
      IDX(3,2)   = KIDX(2,1) + nspoke2(2,m2)
      IDX(3,3)   = KIDX(3,1)
      w  (1:3)   = SWS (1:3,1)

      np = 3

   Case (DWD6_ICON)
      PRLon(1) = PLon
      PRLat(1) = PLat
      lmode    = 3
      if (present (order)) then
         if (order >= 4)  lmode = 6
      end if
      if (present (mode)) lmode = mode
      call search_icon_global (grid% icongrid, &! <-- Interpolation metadata
                               PRLon,          &! <-- Point longitude [deg]
                               PRLat,          &! <-- Point latitude  [deg]
                               jl, jb,         &! --> Line and block indices
                               WS,             &! --> Interpolation weights
                               ngp,            &! --> Number of points returned
                               lmode,          &! <~~ Return 3 or 6 neighbours
                               l_triangle      )! <-- fixed triangle interpolation
      np          = ngp(1)

      idx (1   ,1) = 0
      idx (1:np,1) = jl (1:np,1)
      idx (1:np,2) = jb (1:np,1)
      idx (1:np,3) = 1
      w   (1:np)   = WS (1:np,1)

      if (idx (1,1) <= 0) then
         idx = 0
         w   = 0._wp
         np  = 0
         return
      end if
   Case default
      call finish ('Grid_Indices','invalid gridtype')
End Select

   !---------------------------
   ! Nearest model column only?
   !---------------------------
   if (nn) then
      i = maxloc (w(1:np), dim=1)
      idx(1,1:3) = idx(i,1:3)
      w  (:)     = 0.0_wp       ! Trivial weights
      w  (1)     = 1.0_wp
      np         = 1
   end if

   !--- determine unique index and processor element
   if (lu) then
     do igp = 1, np
       pe_i1_i2_id  = grid% marr (:,IDX(igp,1),IDX(igp,2),IDX(igp,3))
       IDX(igp,1:3) = pe_i1_i2_id(2:4)
       IDX(igp,4)   = pe_i1_i2_id(1)
     end do
   endif

End Subroutine Grid_Indices_geo
!==============================================================================
Subroutine Setup_Interpolation &
  (slon,          & ! <-- Point longitudes [deg]
   slat,          & ! <-- Point latitudes [deg]
   g,             & ! <-- grid data type
   jg,            & ! --> Triangle mesh indices
   SWS)             ! --> Symplectic weights for triangle vertices
!
! Find triangle grid meshes containing points
! and corresponding symplectic weights.
!----------------------------------------------------------
! Method:
!   Binary search.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 06 Dec 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(wp), Intent(In) :: &
   slon(:)     ! Point longitude [deg]
!
Real(wp), Intent(In) :: &
   slat(:)     ! Point latitude [deg]
!
Type(t_grid), Intent(In) :: &
   g           ! grid description
!
! Output arguments:
!
Integer, Intent(Out) :: &
   jg(4,size(slon)) ! Triangle mesh indices
                    ! jg(1:3) - (j1,j2,jd) - top vertex indices
                    ! jg(4)   - mt - triangle number
!
Real(wp), Intent(Out) :: &
   SWS(3,size(slon))    ! Symplectic weights for triangle vertices
!----------------------------------------------------------
! Local Parameters:
!
!
! Local Scalars:
!
Integer :: ji         ! Binary search index
Integer :: IP         ! Interpolation points index
Integer :: ipa        ! Pole-/antipoleward index
Integer :: iv         ! Vertex index
!
! Derived from input parameters
!
Integer :: np                        ! Number of interpolation points
Integer :: ni                        ! Number of segments per diamond
Integer :: ni2                       ! Factorization of ni (factor 2)
Integer :: nir                       ! Factorization of ni (remaining factor)
Integer :: nd                        ! Number of diamonds
real(wp) ,pointer :: xnglob(:,:,:,:) ! locations of gridpoints on unitsphere
!
! Local Arrays:
!
Integer :: &
   JICO(2,3,2)                ! Icosahedral triangle vertices
                              ! [pole-/antipoleward,vertex,index]
Integer :: jd(size(slon))     ! Diamond index
Integer :: kpa(size(slon))    ! Pole-/antipoleward triangle discriminator
!
Real(wp) :: &
   zr(3,size(slon))           ! Cartesian coordinates of points
Integer      :: &
   jb(3,2,size(slon)),       &! Vertex indices of big triangle
   js(3,2,size(slon)),       &! Vertex indices of sub-triangle
                              !    (Top/Left/Right,j1/j2)
!  jtv(size(slon))            ! Top vertex index
   jtv                        ! Top vertex index
Real(wp) :: &
!  rtv(g%nd,2,3,3),          &! Icosahedral triangle vertices
!                             ! [diamond, pole-/antipoleward,vertex,component]
   rtv_(g%nd,2,3,3),         &! Icosahedral triangle vertices
                              ! [diamond, pole-/antipoleward,component,vertex]
   SWB(g%nd,2,3,size(slon)), &! Symplectic weights of vertices
                              ! [diamond, pole-/antipoleward,vertex,point]
   SWMB(g%nd,2,size(slon))    ! Minimum symplectic weight of triangle
                              ! [diamond,pole-/antipoleward,point]
Real(wp) :: &
   SWM(size(slon))            ! Minimum symplectic weights
                              ! of triangles containing points
Integer  :: &
!  ib(2,size(slon))           ! Triangle indices (jd, kpa)
   ib_s(2)                    ! Triangle indices (jd, kpa)
Real(wp) :: &
!  spv(3,size(slon))          ! Scalar product for vertices
   spv_s(3)                   ! Scalar product for vertices
Integer      :: &
   Spk(2:3,2,6)               ! Spokes for 2nd and 3rd vertices
!----------------------------------------------------------

!------------------------
! Derive input parameters
!------------------------
NP     =  size (slon)
ni     =  g% ni
ni2    =  g% ni2
nir    =  g% nir
nd     =  g% nd
xnglob => g% xnglob

!----------------------------------------------------------
! 1. ICOSAHEDRAL TRIANGLE SEARCH
!----------------------------------------------------------


!--- 1.1. Interpolation points in Cartesian coordinates

Do IP=1,NP
   zr(:,IP) = (/ COS(slon(IP)*dtr)*COS(slat(IP)*dtr), &
                 SIN(slon(IP)*dtr)*COS(slat(IP)*dtr), &
                 SIN(slat(IP)*dtr) /)
End Do


!--- 1.2. Icosahedral triangle vertices

JICO(:,:,:) = Reshape &
      (Source =(/ (/  0,    1 /), (/ ni,    1 /), (/  0, ni+1 /),     &
                  (/ ni, ni+1 /), (/  0, ni+1 /), (/ ni,    1 /) /),  &
       Shape = (/ 2, 3, 2 /), Order = (/ 3, 2, 1 /) )

!FTRACE_BEGIN("Setup_Interpolation:1")
Do ipa=1,2
   Do iv=1,3
      rtv_(:,ipa,:,iv) = Transpose(xnglob(JICO(ipa,iv,1),JICO(ipa,iv,2),:,:))
   End Do
End Do
!FTRACE_END  ("Setup_Interpolation:1")


!--- 1.3. Icosahedral triangle search

!FTRACE_BEGIN("Setup_Interpolation:Symplectic_Weights_ICO")
!CDIR IEXPAND(Symplectic_Weights_ICO)
Call Symplectic_Weights_ICO &
  (zr(:,1:NP),              & ! <-- Point in 3D space
   rtv_(1:nd,1:2,:,1),      & ! <-- 1st vertices
   rtv_(1:nd,1:2,:,2),      & ! <-- 2nd vertices
   rtv_(1:nd,1:2,:,3),      & ! <-- 3rd vertices
   SWB (1:nd,1:2,1:3,1:NP))   ! --> Array of symplectic weights
!FTRACE_END  ("Setup_Interpolation:Symplectic_Weights_ICO")

SWMB(1:nd,1:2,1:NP) = MinVal(SWB(1:nd,1:2,:,1:NP),Dim=3)

Do IP=1,NP
   ib_s(:)  = MaxLoc(SWMB(1:nd,1:2,IP))
   jd(IP)   = ib_s(1)
   kpa(IP)  = ib_s(2)
   jg(3,IP) = jd(IP)
End Do


!----------------------------------------------------------
! 2. TRIANGLE SEARCH
!----------------------------------------------------------

!--- 2.1. Initialization of search

Do IP=1,NP
   jb(:,:,IP) = JICO(kpa(IP),:,:)
End Do


!--- 2.2. Binary search

select case (nir)
case default
  call finish ('Setup_Interpolation','grid root is not 1 or 3')
case (1)
case (3)

FTRACE_BEGIN("Setup_Interpolation:Search_Triangle1")
   Call Search_Triangle &
     (NP,               & ! <-- Number of interpolation points
      zr,               & ! <-- Cartesian coordinates of point
      9,                & ! <-- Number of subtriangles
      jd,               & ! <-- Diamond index
      jb,               & ! <-- Vertex indices of big triangle
      xnglob,           & ! <-- gridpoint locations on unit sphere
      js,               & ! --> Vertex indices of sub-triangle
      SWS)                ! --> Symplectic weights for sub-triangle

   jb(:,:,:) = js(:,:,:)
FTRACE_END  ("Setup_Interpolation:Search_Triangle1")

end select

FTRACE_BEGIN("Setup_Interpolation:Search_Triangle2")
Do ji=1,ni2

   Call Search_Triangle &
     (NP,               & ! <-- Number of interpolation points
      zr,               & ! <-- Cartesian coordinates of point
      4,                & ! <-- Number of subtriangles
      jd,               & ! <-- Diamond index
      jb,               & ! <-- Vertex indices of big triangle
      xnglob,           & ! <-- gridpoint locations on unit sphere
      js,               & ! --> Vertex indices of sub-triangle
      SWS)                ! --> Symplectic weights for sub-triangle

   jb(:,:,:) = js(:,:,:)

End Do
FTRACE_END  ("Setup_Interpolation:Search_Triangle2")

SWM(:) = MinVal(SWS(:,:),Dim=1)

!If (MinVal(SWM(:)) < -16000*Epsilon(SWM)) then
 If (MinVal(SWM(:)) < -1.e-11_wp) then
   Write(*,'(A)') 'Warning: inaccurate triangle search!'
   Do IP=1,NP
!     If (SWM(IP) < -16000*Epsilon(SWM)) then
      If (SWM(IP) < -1.e-11_wp) then
         Write(*,'(I3,2(3X,A,F8.3))') &
            IP, 'Slon = ', Slon(IP), 'Slat = ', Slat(IP)
         Write(*,'(A,3(ES11.4,1X))') 'SWS = ', SWS(:,IP)
      End If
   End Do
   call message ('Setup_Interpolation (mo_grid_intpol)',&
                 'Warning: inaccurate triangle search')
End If


!--- 2.3. Setting top to nearest grid point.

Do IP=1,NP

   Do iv=1,3
      spv_s(iv) = Sum(zr(:,IP)*xnglob(jb(iv,1,IP),jb(iv,2,IP),:,jd(IP)))
   End Do

   jtv          = Sum(MaxLoc(spv_s(:)))
   if (jtv /= 1) then ! workaround for IBM xlf 8.1.1.1 bug (rus4)
     js(:,:,IP) = CShift(js(:,:,IP), Shift=jtv-1, Dim=1)
     jb(:,:,IP) = CShift(jb(:,:,IP), Shift=jtv-1, Dim=1)
     SWS(:,IP)  = CShift(SWS(:,IP),  Shift=jtv-1)
   endif
   jg(1:2,IP) = jb(1,1:2,IP)

End Do


!--- 2.4. Determination of triangle number

Spk(2,:,:) = nspoke2(:,:)
Spk(3,:,:) = CShift(nspoke2(:,:), Shift=1, Dim=2)

Do IP=1,NP
   jg(4,IP) = Sum(MinLoc(nspoke2(1,:), &
                Mask = ((Spk(2,1,:) == jb(2,1,IP) - jb(1,1,IP)) .and. &
                        (Spk(2,2,:) == jb(2,2,IP) - jb(1,2,IP)) .and. &
                        (Spk(3,1,:) == jb(3,1,IP) - jb(1,1,IP)) .and. &
                        (Spk(3,2,:) == jb(3,2,IP) - jb(1,2,IP)))))

End Do

End Subroutine Setup_Interpolation



!==========================================================
Subroutine Search_Triangle &
  (NP,               & ! <-- Number of interpolation points
   zr,               & ! <-- Cartesian coordinates of point
   NS,               & ! <-- Number of subtriangles
   jd,               & ! <-- Diamond index
   jb,               & ! <-- Vertex indices of big triangle
   xnglob,           & ! <-- gridpoint locations on unit sphere
   js,               & ! --> Vertex indices of sub-triangle
   SWS)                ! --> Symplectic weights for sub-triangle
!
! Determine sub-triangle containing point. The sub-triangle
! is one of four or nine which partition big triangle.
!----------------------------------------------------------
! Method:
!   Search for big triangle with nearest center;
!   search for triangle with maximum symplectic weight.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 05 Dec 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   NP               ! Number of interpolation points
!
Real(wp), Intent(In) :: &
   zr(3,NP)         ! Cartesian coordinates of point
!
Integer, Intent(In) :: &
   NS               ! Number of subtriangles
!
Integer, Intent(In) :: &
   jd(NP)           ! Diamond index
!
Integer, Intent(In) :: &
   jb(3,2,NP)       ! Vertex indices of big triangle
                    ! (Top/Left/Right,j1/j2,Point)
!
Real(wp), Intent(In) :: &
   xnglob(0:,:,:,:) ! locations of gridpoints on unit sphere
!
! Output arguments:
!
Integer, Intent(Out) :: &
   js(3,2,NP)       ! Vertex indices of sub-triangle
                    ! (Top/Left/Right,j1/j2,Point)
!
Real(wp), Intent(Out) :: &
   SWS(3,NP)        ! Symplectic weights for sub-triangle
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: IP    ! Point number
Integer      :: iv    ! Vertex index
Integer      :: id    ! Diamond index
!
! Local Arrays:
!
Integer      ::  &
   ide(NP)          ! Subtriangle index
Real(wp) :: &
   SW(NS,3,NP),   & ! Symplectic weights
!  SWM(NS,NP)       ! Minimum symplectic weight
   SWM(NS)          ! Minimum symplectic weight
Integer :: &
   j(10,2,NP)       ! Subtriangle vertex indices
Real(wp) ::  &
!  rtv(NS,3,3,NP)   ! Subtriangle vertices (original dimension order)
                    ! (Sub-triangle,vertex,component,point)
   rtv_x(NS,3,NP,3) ! Subtriangle vertices (modified dimension order [ha])
                    ! (Sub-triangle,component,point,vertex)
Integer, Pointer :: &
   SI(:,:)          ! Pointer for sub-triangle index array
!----------------------------------------------------------
Integer, Target :: &
   SI9(9,3) = &     ! Sub-triangle indices
      Reshape( (/ 1,  2,  3,       &  ! D1
                  2,  4,  5,       &  ! D2
                  4,  7,  8,       &  ! D3
                  5,  8,  9,       &  ! D4
                  6,  9, 10,       &  ! D5
                  3,  5,  6,       &  ! D6
                  5,  3,  2,       &  ! D7
                  8,  5,  4,       &  ! D8
                  9,  6,  5  /),   &  ! D9
               Shape = (/ 9, 3 /), &
               Order = (/ 2, 1 /) )
!
! Poleward triangle                 1(1)
!                                   / \
!                                  / D1\
!                                 /     \
!                                2-------3
!                               / \     / \
!                              / D2\ D7/ D6\
!                             /     \ /     \
!                            4-------5-------6
!                           / \     / \     / \
!                          / D3\ D8/ D4\ D9/ D5\
!                         /     \ /     \ /     \
!                       7(2)-----8-------9-----10(3)
!
!-----------------------------------------------------------------------
Integer, Target :: &
   SI4(4,3) = &     ! Sub-triangle indices
      Reshape( (/ 1,  2,  3,       &  ! D1
                  2,  4,  5,       &  ! D2
                  3,  5,  6,       &  ! D3
                  5,  3,  2 /),    &  ! D4
               Shape = (/ 4, 3 /), &
               Order = (/ 2, 1 /) )
!
! Poleward triangle                 1(1)
!                                   / \
!                                  / D1\
!                                 /     \
!                                2------ 3
!                               / \     / \
!                              / D2\ D4/ D3\
!                             /     \ /     \
!                           4(2)-----5------6(3)
!
!-----------------------------------------------------------------------
! Antipoleward triangle: top is oriented away from pole,
!                        lower left - eastern, lower right - western.
!=======================================================================

!CDIRR ON_ADB(j)
!CDIRR ON_ADB(rtv_x)

!----------------------------------------------------------
! 1. SUB-TRIANGLE DEFINITION
!----------------------------------------------------------

!--- 1.1. Vertices

Select Case (NS)

   Case (4)

      j( 1,:,:)  = jb(1,:,:)
      j( 4,:,:)  = jb(2,:,:)
      j( 6,:,:)  = jb(3,:,:)
      j( 2,:,:)  = (j(1,:,:) + j(4,:,:))/2
      j( 3,:,:)  = (j(1,:,:) + j(6,:,:))/2
      j( 5,:,:)  = (j(4,:,:) + j(6,:,:))/2

      SI => SI4

   Case (9)

      j( 1,:,:)  = jb(1,:,:)
      j( 7,:,:)  = jb(2,:,:)
      j(10,:,:)  = jb(3,:,:)
      j( 2,:,:)  =  j(1,:,:) +   (j( 7,:,:) - j(1,:,:))/3
      j( 4,:,:)  =  j(1,:,:) + 2*(j( 7,:,:) - j(1,:,:))/3
      j( 3,:,:)  =  j(1,:,:) +   (j(10,:,:) - j(1,:,:))/3
      j( 6,:,:)  =  j(1,:,:) + 2*(j(10,:,:) - j(1,:,:))/3
      j( 8,:,:)  =  j(7,:,:) +   (j(10,:,:) - j(7,:,:))/3
      j( 9,:,:)  =  j(7,:,:) + 2*(j(10,:,:) - j(7,:,:))/3
      j( 5,:,:)  = (j(8,:,:) + j(3,:,:))/2

      SI => SI9

   Case Default

      Write(*,'(A,I0/A)') 'Search_Triangle: ns = ', ns, &
                          '   must be 4 or 9'
      Stop

End Select

!CDIRR ON_ADB(SI)

!--- 1.2. Centres

Do IP=1,NP
   Do iv=1,3
      Do id=1,NS
         rtv_x(id,:,IP,iv) = &
            xnglob(j(SI(id,iv),1,IP),j(SI(id,iv),2,IP),:,jd(IP))
      End Do
   End Do
End Do


!----------------------------------------------------------
! 2. SUB-TRIANGLE SEARCH
!----------------------------------------------------------

!--- 2.1. Search for triangle with maximum symplectic weight

!CDIR IEXPAND(Symplectic_Weights_Small)
Call Symplectic_Weights_Small &
  (zr(:,:),         & ! <-- Point in 3D space
   rtv_x(:,:,:,1),  & ! <-- 1st vertices
   rtv_x(:,:,:,2),  & ! <-- 2nd vertices
   rtv_x(:,:,:,3),  & ! <-- 3rd vertices
   SW(:,:,:))         ! --> Array of symplectic weights

Do IP=1,NP
   SWM(:)    = MinVal(SW(:,:,IP),Dim=2)
   ide(IP)   = Sum(MaxLoc(SWM(:)))
   SWS(:,IP) = SW(ide(IP),:,IP)
End Do


!--- 2.2. Setting vertices

Do IP=1,NP
   js(:,:,IP) = j(SI(ide(IP),:),:,IP)
End Do

End Subroutine Search_Triangle



!==========================================================
Subroutine Symplectic_Weights_Small &
  (z,       & ! <-- Point in 3D space
   x1,      & ! <-- 1st vertices
   x2,      & ! <-- 2nd vertices
   x3,      & ! <-- 3rd vertices
   SW)        ! --> Array of symplectic weights
!
! Calculation of symplectic weights for
! multiple (point/small triangle set) pairs.
!----------------------------------------------------------
! Method:
!            (z,x2,x3)     (x1,z,x3)     (x1,x2,z)
! SW(1:3) = ------------, ------------, ------------
!                N             N             N
! N is norming constant: Sum(SW(1:3)) = 1
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 05 Dec 2001 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(wp), Intent(In) :: &
   z(:,:)     ! Points
              ! [Component,point number]
!
Real(wp), Intent(In) :: &
   x1(:,:,:)  ! 1st vertex direction
              ! [Triangle number,component,point number]
!
Real(wp), Intent(In) :: &
   x2(:,:,:)  ! 2nd vertex direction
!
Real(wp), Intent(In) :: &
   x3(:,:,:)  ! 3rd vertex direction
!
! Output arguments:
!
Real(wp) :: &
   SW(:,:,:)  ! Array of symplectic weights
              ! [Triangle number,vertex,point number]
!----------------------------------------------------------
! Local Scalars:
!
Integer :: IP  ! Point index
integer :: i
!
! Local Arrays:
!
Real(wp) :: &
   D(Size(SW,1))                ! Symplectic denominator
real(wp) :: y1(3), y2(3), y3(3) ! Temporaries
!----------------------------------------------------------

Do IP=1,Size(SW,3)

#if 0  /* original version */

   SW(:,1,IP) = (  z(1,IP)*x2(:,2,IP)*x3(:,3,IP)    &
                 - z(1,IP)*x2(:,3,IP)*x3(:,2,IP) +  &
                   z(2,IP)*x2(:,3,IP)*x3(:,1,IP)    &
                 - z(2,IP)*x2(:,1,IP)*x3(:,3,IP) +  &
                   z(3,IP)*x2(:,1,IP)*x3(:,2,IP)    &
                 - z(3,IP)*x2(:,2,IP)*x3(:,1,IP))

   SW(:,2,IP) = (  x1(:,1,IP)*z(2,IP)*x3(:,3,IP)    &
                 - x1(:,1,IP)*z(3,IP)*x3(:,2,IP) +  &
                   x1(:,2,IP)*z(3,IP)*x3(:,1,IP)    &
                 - x1(:,2,IP)*z(1,IP)*x3(:,3,IP) +  &
                   x1(:,3,IP)*z(1,IP)*x3(:,2,IP)    &
                 - x1(:,3,IP)*z(2,IP)*x3(:,1,IP))

   SW(:,3,IP) = (  x1(:,1,IP)*x2(:,2,IP)*z(3,IP)    &
                 - x1(:,1,IP)*x2(:,3,IP)*z(2,IP) +  &
                   x1(:,2,IP)*x2(:,3,IP)*z(1,IP)    &
                 - x1(:,2,IP)*x2(:,1,IP)*z(3,IP) +  &
                   x1(:,3,IP)*x2(:,1,IP)*z(2,IP)    &
                 - x1(:,3,IP)*x2(:,2,IP)*z(1,IP))

#else  /* alternative version for testing of sensitivity to rounding */

   !--------------------------------------------------------
   ! Rewrite determinants into numerically more stable form,
   ! using (z,x2,x3) == (z,x2-z,x3-z) etc.
   !--------------------------------------------------------
   do i = 1, size(sw,1)
     y1(:) = x1(i,:,IP) - z(:,IP)
     y2(:) = x2(i,:,IP) - z(:,IP)
     y3(:) = x3(i,:,IP) - z(:,IP)
     SW(i,1,IP) = (  y2(2) * y3(3)              &
                   - y2(3) * y3(2) ) * z (1,IP) &
                + (  y2(3) * y3(1)              &
                   - y2(1) * y3(3) ) * z (2,IP) &
                + (  y2(1) * y3(2)              &
                   - y2(2) * y3(1) ) * z (3,IP)

     SW(i,2,IP) = (  y1(1) * y3(3)              &
                   - y1(3) * y3(1) ) * z (2,IP) &
                + (  y1(2) * y3(1)              &
                   - y1(1) * y3(2))  * z (3,IP) &
                + (  y1(3) * y3(2)              &
                   - y1(2) * y3(3) ) * z (1,IP)

     SW(i,3,IP) = (  y1(1) * y2(2)              &
                   - y1(2) * y2(1) ) * z (3,IP) &
                + (  y1(2) * y2(3)              &
                   - y1(3) * y2(2) ) * z (1,IP) &
                + (  y1(3) * y2(1)              &
                   - y1(1) * y2(3) ) * z (2,IP)
   end do

#endif


#if defined (__SX__)
   ! Explicit expansion of sum() for better vectorization
   D(:)       = SW(:,1,IP) + SW(:,2,IP) + SW(:,3,IP)
#else
   D(:)       = Sum (SW(:,1:3,IP), Dim=2)
#endif

   SW(:,1,IP) = SW(:,1,IP)/D(:)
   SW(:,2,IP) = SW(:,2,IP)/D(:)
   SW(:,3,IP) = SW(:,3,IP)/D(:)

End Do


End Subroutine Symplectic_Weights_Small



!==========================================================
Subroutine Symplectic_Weights_ICO &
  (z,       & ! <-- Point in 3D space
   x1,      & ! <-- 1st vertices
   x2,      & ! <-- 2nd vertices
   x3,      & ! <-- 3rd vertices
   SW)        ! --> Array of symplectic weights
!
! Calculation of symplectic weights for icosahedral
! triangles and multiple points
!----------------------------------------------------------
! Method:
!            (z,x2,x3)     (x1,z,x3)     (x1,x2,z)
! SW(1:3) = ------------, ------------, ------------
!            (x1,x2,x3)    (x1,x2,x3)    (x1,x2,x3)
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 06 Dec 2001 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(wp), Intent(In) :: &
   z(:,:)       ! Points
                ! [Component,point number]
!
Real(wp), Intent(In) :: &
   x1(:,:,:)    ! 1st vertex direction
                ! [Diamond,pole-/antipoleward,component]
!
Real(wp), Intent(In) :: &
   x2(:,:,:)    ! 2nd vertex direction
!
Real(wp), Intent(In) :: &
   x3(:,:,:)    ! 3rd vertex direction
!
! Output arguments:
!
Real(wp) :: &
   SW(:,:,:,:)  ! Array of symplectic weights
                ! [Diamond,pole-/antipoleward,vertex,point number]
!----------------------------------------------------------
! Local Scalars:
!
Integer :: IP  ! Point index
!
! Local Arrays:
!
Real(wp) :: &
   D(Size(SW,1),Size(SW,2))              ! Symplectic denominator
!----------------------------------------------------------

D(:,:)    =  ( x1(:,:,1)*x2(:,:,2)*x3(:,:,3)    &
             - x1(:,:,1)*x2(:,:,3)*x3(:,:,2) +  &
               x1(:,:,2)*x2(:,:,3)*x3(:,:,1)    &
             - x1(:,:,2)*x2(:,:,1)*x3(:,:,3) +  &
               x1(:,:,3)*x2(:,:,1)*x3(:,:,2)    &
             - x1(:,:,3)*x2(:,:,2)*x3(:,:,1))

Do IP=1,Size(SW,4)

   SW(:,:,1,IP) = (  z(1,IP)*x2(:,:,2)*x3(:,:,3)    &
                   - z(1,IP)*x2(:,:,3)*x3(:,:,2) +  &
                     z(2,IP)*x2(:,:,3)*x3(:,:,1)    &
                   - z(2,IP)*x2(:,:,1)*x3(:,:,3) +  &
                     z(3,IP)*x2(:,:,1)*x3(:,:,2)    &
                   - z(3,IP)*x2(:,:,2)*x3(:,:,1)) / D(:,:)

   SW(:,:,2,IP) = (  x1(:,:,1)*z(2,IP)*x3(:,:,3)    &
                   - x1(:,:,1)*z(3,IP)*x3(:,:,2) +  &
                     x1(:,:,2)*z(3,IP)*x3(:,:,1)    &
                   - x1(:,:,2)*z(1,IP)*x3(:,:,3) +  &
                     x1(:,:,3)*z(1,IP)*x3(:,:,2)    &
                   - x1(:,:,3)*z(2,IP)*x3(:,:,1)) / D(:,:)

   SW(:,:,3,IP) = (  x1(:,:,1)*x2(:,:,2)*z(3,IP)    &
                   - x1(:,:,1)*x2(:,:,3)*z(2,IP) +  &
                     x1(:,:,2)*x2(:,:,3)*z(1,IP)    &
                   - x1(:,:,2)*x2(:,:,1)*z(3,IP) +  &
                     x1(:,:,3)*x2(:,:,1)*z(2,IP)    &
                   - x1(:,:,3)*x2(:,:,2)*z(1,IP)) / D(:,:)

End Do


End Subroutine Symplectic_Weights_ICO



!==========================================================
Subroutine Interpolation &
  (NP,        & ! <-- Number of interpolation points
   lb,        & ! <-- Lower bounds of gridded field
   px,        & ! <-- Gridded field
   SWS,       & ! <-- Symplectic weights for triangle vertices of points
   jg,        & ! <-- Index of the triangle containing point
   pxi)         ! --> Interpolated field
!
! Interpolation on icosahedral-hexagonal grid.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 05 Dec 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   NP                    ! Number of interpolation points
!
Integer, Intent(In) :: &
   lb(4)                 ! Lower bounds of gridded field
!
Real(wp), Intent(In) :: &
   px(lb(1):,         &  ! Gridded field
      lb(2):,         &
      lb(4):)
!
Real(wp), Intent(In) :: &
   SWS(3,NP)             ! Symplectic weights for triangle vertices
!
Integer, Intent(In) :: &
   jg(4,NP)              ! Index of the triangle containing point:
                         !  j1,j2,jd, mt
!
! Output arguments:
!
Real(wp), Intent(Out) :: &
   pxi(NP)               ! Interpolated field
!----------------------------------------------------------
! Local Scalars:
!
Real(wp) :: zx1, zx2, zx3   ! values at three gridpoints
Integer      :: IP              ! Interpolation point index
Integer      :: j1, j2, jd      ! Gridpoint indices
Integer      :: m1, m2          ! Triangle indices
!----------------------------------------------------------


Do IP = 1, NP

   j1    = jg(1,IP)
   j2    = jg(2,IP)
   jd    = jg(3,IP)

   m1    = jg(4,IP)
   m2    = MOD (m1,6) + 1

   zx1 = px(j1,j2,jd)
   zx2 = px(j1+nspoke2(1,m1), j2+nspoke2(2,m1), jd)
   zx3 = px(j1+nspoke2(1,m2), j2+nspoke2(2,m2), jd)

   pxi(IP) = SWS(1,IP)*zx1 + SWS(2,IP)*zx2 + SWS(3,IP)*zx3

End Do

End Subroutine Interpolation
!==============================================================================

Subroutine Longitude_Interpolation &
  (XLon,      & ! <-- Longitude grid [deg]
   PLon,      & ! <-- Longitude of point [deg]
   Cyc_x,     & ! <-- Cyclic boundary conditions
   ILon,      & ! --> Interpolation subgrid
   WLon,      & ! --> Weights of subgrid points
   WLon1,     & ! ~~> Weight derivatives
   WLon2)       ! ~~> Weight 2nd derivatives
!
! Longitude interpolation: calculation of interpolation
! weights of grid points and, optionally, the 1st and 2nd
! derivatives of weighting function.
!----------------------------------------------------------
! Method:
!   For N-point scheme, the central cell containing
!   interpolation point is found. Coefficients of
!   polynomial of N-1 power are found so as for
!   derivatives up to (N-2)th to be continuous.
!----------------------------------------------------------
! (C) Copyright 1998-99, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 13 Nov 1998 | Original code.
!   2.0   ! 20 Feb 1999 ! Faster index search.
!----------------------------------------------------------
! Modules used:
!
Use Interpolation, only: &
! Imported Routines:
    Seek_Index
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(wp), Intent(In) :: &
   XLon(1:)   ! Longitude grid [deg]
!
Real(wp), Intent(In) :: &
   PLon       ! Longitude of point [deg]
!
Logical, Intent(In) :: &
   Cyc_x      ! Cyclic boundary conditions
!
! Output arguments:
!
Integer, Intent(Out) :: &
   ILon(1:)   ! Indices of grid points used
              ! Size(ILon) must be even
!
Real(wp), Intent(Out) :: &
   WLon(1:)   ! Weights of grid points used
              ! Size(WLon) must be = Size(ILon)
!
Real(wp), Intent(Out), Optional :: &
   WLon1(1:)  ! Weight derivatives
              ! Size(WLon1) must be = Size(ILon)
!
Real(wp), Intent(Out), Optional :: &
   WLon2(1:)  ! Weight 2nd derivatives
              ! Size(WLon2) must be = Size(ILon)
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: NLon  ! Longitude grid dimension
Integer      :: DLon  ! Longitude grid direction
Integer      :: NI    ! Interpolation scheme dimension
Integer      :: i     ! Grid index
Integer      :: IC    ! Lower index of central cell of subgrid
Integer      :: B     ! Border indicator
Real(wp)     :: TLon  ! Longitude in 0-360 deg interval
!
! epsilon for out of bound check
!
Real(wp) ,parameter :: epsilon = 1.e-13_wp
!
! Local Arrays:
!
Real(wp) :: &
   SLon(Size(ILon))    ! Interpolation subgrid
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF GRID SIZE AND DIRECTION
!----------------------------------------------------------


!--- 1.1. Main grid

NLon = Size(XLon)
DLon = Nint(Sign(1.0_wp, XLon(NLon)-XLon(1)))


!--- 1.2. Interpolation subgrid

NI   = Size(ILon)


!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATION CELL
!----------------------------------------------------------


!--- 2.1. Transform of longitude to 0-360 interval

!--- Out Of Domain Check

If(Cyc_x) then
  TLon = XLon(1) + Modulo(PLon-XLon(1), 360.0_wp)
Else
   If (DLon*PLon > DLon*XLon(NLon) + epsilon) then
      ILon = 0
      WLon = 0._wp
      Return
   Else If (DLon*PLon < DLon*XLon(1) - epsilon) then
      ILon = 0
      WLon = 0._wp
      Return
   Endif
   TLon = PLon
Endif



!--- 2.2. Calculation of central cell of interpolation subgrid

IC       = NI/2

!ILon(IC) = Sum(MaxLoc(DLon*Xlon(:), &
!                      DLon*XLon(:) <= DLon*TLon))
If(Cyc_x) then
  If (DLon*TLon >= DLon*XLon(NLon)) then
     ILon(IC) = NLon
  Else
     ILon(IC) = Seek_Index (XLon, Real(TLon, wp))
  End If
  ILon(IC) = Max(1, ILon(IC))
Else
  If (DLon*TLon >= DLon*XLon(NLon)) then
     ILon(IC) = NLon - 1
  Else If (DLon*TLon <= DLon*XLon(1)) then
     ILon(IC) = 1
  Else
     ILon(IC) = Seek_Index (XLon, Real(TLon, wp))
  End If
  ILon(IC) = Min(Max(1, ILon(IC)), NLon-1)
Endif
SLon(IC) = XLon(ILon(IC))


!--- 2.3. Calculation of interpolation subgrid

If(Cyc_x) then
  Do i=1,NI
     ILon(i) = 1  + Modulo(ILon(IC)+i-IC-1, NLon)
     SLon(i) = SLon(IC)-180.0_wp + &
        Modulo(XLon(ILon(i))-SLon(IC)+180.0_wp, 360.0_wp)
  End Do
Else
  Do i=1,NI
     ILon(i) = Min(Max(1, ILon(IC)+i-IC), NLon)
     SLon(i) = XLon(ILon(i))
  End Do
Endif

!--- 2.4 Determination of border position

if (Cyc_x) then
   B = 0
else
   If (ILon(IC) == 1) then
      B = -1
   Else If (ILon(IC) == NLon-1) then
      B = 1
   Else
      B = 0
   END If
end if

!----------------------------------------------------------
! 3. CALCULATION OF INTERPOLATION WEIGHTS
!----------------------------------------------------------


Select Case (NI)


!--- 3.1. 2-point linear interpolation

Case(2)

   WLon(1) = (SLon(2) - TLon)/(SLon(2)-SLon(1))
   WLon(2) = (TLon - SLon(1))/(SLon(2)-SLon(1))

   If (Present(WLon1)) then
      WLon1(1) = -1._wp / (SLon(2)-SLon(1))
      WLon1(2) =  1._wp / (SLon(2)-SLon(1))
   End If

   If (Present(WLon2)) then
      WLon2(:) = 0
   End If


!--- 3.2. 4-point linear interpolation

Case(4)

   Call Weight4 &
     (TLon,      & ! <-- Interpolation point
      SLon,      & ! <-- Interpolation subgrid
      B,         & ! <-- Border indicator
      WLon,      & ! --> Interpolation weights
      WLon1)       ! ~~> Interpolation weight derivatives

   If (Present(WLon2)) then
      WLon2(:) = 0
   End If


!--- 3.3. Default

Case Default

   WLon(:) = 0

   If (Present(WLon1)) then
      WLon1(:) = 0
   End If

   If (Present(WLon2)) then
      WLon2(:) = 0
   End If


End Select


End Subroutine Longitude_Interpolation



!==========================================================
Subroutine Latitude_Interpolation &
  (XLat,      & ! <-- Latitude grid [deg]
   PLat,      & ! <-- Latitude of point [deg]
   Poly,      & ! <-- Poles at bounds
   ILat,      & ! --> Interpolation subgrid
   WLat,      & ! --> Weights of subgrid points
   WLat1,     & ! ~~> Weight derivatives
   WLat2)       ! ~~> Weight 2nd derivatives
!
! Latitude interpolation: calculation of interpolation
! weights of grid points and, optionally, the 1st and 2nd
! derivatives of weighting function.
!----------------------------------------------------------
! Method:
!   For N-point scheme, the central cell containing
!   interpolation point is found. Coefficients of
!   polynomial of N-1 power are found so as for
!   derivatives up to (N-2)th to be continuous.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 13 Nov 1998 | Original code.
!   2.0   ! 20 Feb 1999 ! Faster index search.
!----------------------------------------------------------
! Modules used:
!
Use Interpolation, only: &
! Imported Routines:
    Seek_Index
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(wp), Intent(In) :: &
   XLat(1:)   ! Latitude grid [deg]
!
Real(wp), Intent(In) :: &
   PLat       ! Latitude of point [deg]
!
Logical, Intent(In) :: &
   Poly       ! poles at bounds
!
! Output arguments:
!
Integer, Intent(Out) :: &
   ILat(1:)   ! Indices of grid points used
              ! Size(ILat) must be even
!
Real(wp), Intent(Out) :: &
   WLat(1:)   ! Weights of grid points used
              ! Size(WLat) must be = Size(ILat)
!
Real(wp), Intent(Out), Optional :: &
   WLat1(1:)  ! Weight derivatives
              ! Size(WLat1) must be = Size(ILat)
!
Real(wp), Intent(Out), Optional :: &
   WLat2(1:)  ! Weight 2nd derivatives
              ! Size(WLat2) must be = Size(ILat)
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: NLat  ! Latitude grid dimension
Integer      :: DLat  ! Latitude grid direction
Integer      :: NI    ! Interpolation scheme dimension
Integer      :: i     ! Grid index
Integer      :: IC    ! Lower index of central cell of subgrid
Integer      :: B     ! Border indicator
!
! epsilon for out of bound check
!
Real(wp) ,parameter :: epsilon = 1.e-13_wp
!
! Local Arrays:
!
Real(wp) :: &
   SLat(Size(ILat))   ! Interpolation subgrid
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF GRID SIZE AND DIRECTION
!----------------------------------------------------------


!--- 1.1. Main grid

NLat = Size(XLat)
DLat = Nint(Sign(1.0_wp, XLat(NLat)-XLat(1)))


!--- 1.2. Interpolation subgrid

NI   = Size(ILat)


!--- 1.3. Out Of Domain Check

If(.not.Poly) then
   If (DLat*PLat > DLat*XLat(NLat) + epsilon) then
      ILat = 0
      WLat = 0._wp
      Return
   Else If (DLat*PLat < DLat*XLat(1) - epsilon) then
      ILat = 0
      WLat = 0._wp
      Return
   Endif
Endif

!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATION CELL
!----------------------------------------------------------


!--- 2.1. Calculation of central cell of interpolation subgrid

IC       = NI/2

!ILat(IC) = Sum(MaxLoc(DLat*XLat(:), &
!                      DLat*XLat(:) <= DLat*PLat))
If (DLat*PLat >= DLat*XLat(NLat)) then
   ILat(IC) = NLat - 1
Else If (DLat*PLat <= DLat*XLat(1)) then
   ILat(IC) = 1
Else
   ILat(IC) = Seek_Index (XLat, Real(PLat, wp))
End If
ILat(IC) = Min(Max(1, ILat(IC)), NLat-1)
SLat(IC) = XLat(ILat(IC))


!--- 2.2. Calculation of interpolation subgrid

Do i=1,NI
   ILat(i) = Min(Max(1, ILat(IC)+i-IC), NLat)
   SLat(i) = XLat(ILat(i))
End Do


!--- 2.3. Determination of border position

If (ILat(IC) == 1) then
   B = -1
Else If (ILat(IC) == NLat-1) then
   B = 1
Else
   B = 0
End If


!----------------------------------------------------------
! 3. CALCULATION OF INTERPOLATION WEIGHTS
!----------------------------------------------------------


Select Case (NI)


!--- 3.1. 2-point linear interpolation

Case(2)

   WLat(1) = (SLat(2) - PLat)/(SLat(2)-SLat(1))
   WLat(2) = (PLat - SLat(1))/(SLat(2)-SLat(1))

   If (Present(WLat1)) then
      WLat1(1) = -1._wp /(SLat(2)-SLat(1))
      WLat1(2) =  1._wp /(SLat(2)-SLat(1))
   End If

   If (Present(WLat2)) then
      WLat2(:) = 0
   End If


Case(4)

   Call Weight4 &
     (PLat,    & ! <-- Interpolation point
      SLat,    & ! <-- Interpolation subgrid
      B,       & ! <-- Border indicator
      WLat,    & ! --> Interpolation weights
      WLat1)     ! ~~> Interpolation weight derivatives

   If (Present(WLat2)) then
      WLat2(:) = 0
   End If


!--- 3.x. Default

Case Default

   WLat(:) = 0

   If (Present(WLat1)) then
      WLat1(:) = 0
   End If

   If (Present(WLat2)) then
      WLat2(:) = 0
   End If


End Select


End Subroutine Latitude_Interpolation



!==========================================================
Subroutine Weight4 &
  (XP,        & ! <-- Interpolation point
   X,         & ! <-- Interpolation subgrid
   B,         & ! <-- Border indicator
   W,         & ! --> Interpolation weights
   DW)          ! ~~> Interpolation weight derivatives
!
! Calculation of weights for 4-point interpolation.
!----------------------------------------------------------
! Method:
!   Cubic interpolation between 2nd and 3rd points using
!   2 derivatives from finite differences and 2 values
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov, A. K. Steiner.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 09 Mar 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(wp), Intent(In) :: &
   XP       ! Interpolation point
!
Real(wp), Intent(In) :: &
   X(:)     ! Interpolation subgrid
!
Integer, Intent(In) :: &
   B        ! Border indicator
            !   -1 - left border
            !    0 - middle
            !    1 - right border
!
! Output arguments:
!
Real(wp), Intent(Out) :: &
   W(:)     ! Interpolation weights
!
Real(wp), Optional, Intent(Out) :: &
   DW(:)    ! Interpolation weights
!----------------------------------------------------------



Select Case(B)

   Case (0)
      W(1) = 0.5*(XP-X(2))*(XP-X(3))**2/((X(1)-X(2))*(X(2)-X(3))**2)

      W(2) = ((XP-X(3))/(X(2)-X(3))**3)*     &
             ((XP-X(3))*(3*X(2)-X(3)-2*XP) + &
                0.5*(XP-X(2))*((XP-X(3))*(X(3)-2*X(2)+X(1))/(X(1)-X(2)) + &
                XP-X(2)))

      W(3) = ((XP-X(2))/(X(3)-X(2))**3)*     &
             ((XP-X(2))*(3*X(3)-X(2)-2*XP) + &
                0.5*(XP-X(3))*((XP-X(2))*(X(4)-2*X(3)+X(2))/(X(4)-X(3)) + &
                XP-X(3)))

      W(4) = 0.5*(XP-X(3))*(XP-X(2))**2/((X(4)-X(3))*(X(3)-X(2))**2)

      If (Present(DW)) then

         DW(1) = ((2*X(2) + X(3) - 3*XP)*(X(3) - XP))/   &
                 (2*(X(1) - X(2))*(X(2) - X(3))**2)

         DW(2) = (-X(2)**3 + X(2)**2*(6*X(3) - 4*XP) +   &
                 X(2)*XP*(-4*X(3) + 3*XP) +              &
                 X(1)*(X(2)**2 - 8*X(2)*X(3) + X(3)**2 + &
                 6*X(2)*XP + 6*X(3)*XP - 6*XP**2) +      &
                 X(3)*(X(3)**2 - 4*X(3)*XP + 3*XP**2))/  &
                 (2*(X(1) - X(2))*(X(2) - X(3))**3)

         DW(3) = (X(2)**3 - X(3)**3 + X(2)**2*(X(4) - 4*XP) + &
                 X(3)**2*(X(4) - 4*XP) - 6*X(4)*XP**2 +       &
                 3*X(3)*XP*(2*X(4) + XP) + X(2)*(6*X(3)**2 -  &
                 4*X(3)*(2*X(4) + XP) + &
                 3*XP*(2*X(4) + XP)))/(2*(X(2) - X(3))**3*(X(3) - X(4)))

         DW(4) = -(((X(2) + 2*X(3) - 3*XP)*(X(2) - XP))/ &
                  (2*(X(2) - X(3))**2*(X(3) - X(4))))

      End If

   Case(-1)

      W(1) = 0

      W(2) = 0.5*(XP-X(3))*(X(2)-2*X(3)+XP)/(X(3)-X(2))**2

      W(3) = ((XP-X(2))/(X(3)-X(2))**2)*                &
             (2*X(3)-X(2)-XP + 0.5*(2*X(3)-X(4)-X(2))*  &
             (XP-X(3))/(X(3)-X(4)))

      W(4) = 0.5*(XP-X(3))*(XP-X(2))/((X(4)-X(3))*(X(3)-X(2)))

      If (Present(DW)) then

         DW(1) = 0

         DW(2) = (X(2) - 3*X(3) + 2*XP)/(2*(X(2) - X(3))**2)

         DW(3) = (X(2)**2 + 2*X(3)**2 - 3*X(3)*X(4) + X(2)* &
           (-X(3) + X(4) - 2*XP) + 2*X(4)*XP)/              &
           (2*(X(2) - X(3))**2*(X(3) - X(4)))

         DW(4) = -((X(2) + X(3) - 2*XP)/ &
                  (2*(X(2) - X(3))*(X(3) - X(4))))

      End If

   Case(1)

      W(1) = 0.5*(XP-X(2))*(XP-X(3))/((X(1)-X(2))*(X(2)-X(3)))

      W(2) = ((XP-X(3))/(X(2)-X(3))**2)*                &
             (2*X(2)-X(3)-XP + 0.5*(2*X(2)-X(3)-X(1))*  &
             (XP-X(2))/(X(2)-X(1)))

      W(3) = 0.5*(XP-X(2))*(X(3)-2*X(2)+XP)/(X(3)-X(2))**2

      W(4) = 0

      If (Present(DW)) then

         DW(1) = -((X(2) + X(3) - 2*XP)/(2*(X(1) - X(2))*(X(2) - X(3))))

         DW(2) = (-2*X(2)**2 + X(2)*X(3) + X(1)*                &
                 (3*X(2) - X(3) - 2*XP) - X(3)*(X(3) - 2*XP))/  &
                 (2*(X(1) - X(2))*(X(2) - X(3))**2)

         DW(3) = (-3*X(2) + X(3) + 2*XP)/(2*(X(2) - X(3))**2)

         DW(4) = 0

      End If

End Select

!DEBUG - NO HORIZONTAL GRADIENTS
!DW(:) = 0


End Subroutine Weight4
!==============================================================================
  subroutine hor_intp_coeff (gout, gin, idx, idxin,                 &
                             nn, sti, sto, di, do, stc1, stc2, dist2)
  type(t_grid) ,intent(in) :: gout         ! destination grid
  type(t_grid) ,intent(in) :: gin          ! source grid
  type(t_h_idx),pointer    :: idx   (:,:,:)! interpolation coefficients
  type(t_h_idx),intent(in) :: idxin (:,:,:)! coefficients to re-use
  logical      ,intent(in) :: nn           ! nearest neighbour
  integer      ,intent(in) :: sti (gin %lbg(1):,gin %lbg(2):,:)! generalised
  integer      ,intent(in) :: sto (gout%lbg(1):,gout%lbg(2):,:)! .. soil types
  real(wp)     ,intent(in) :: di  (gin %lbg(1):,gin %lbg(2):,:)! for 'distance'
  real(wp)     ,intent(in) :: do  (gout%lbg(1):,gout%lbg(2):,:)! .. calculation
  logical      ,intent(in) :: stc1  (:,:)  ! soil type classes allowed (1st try)
  logical      ,intent(in) :: stc2  (:,:)  ! soil type classes allowed (2nd try)
  real(wp)     ,intent(in) :: dist2        ! max.distance allowed (km) (3nd try)
  optional                 :: idxin, nn, sti, sto, di, do, stc1, stc2, dist2
  !------------------------------------------------------------
  ! calculate coefficients for grid-to grid interpolations with
  ! options for nearest neighbour, same soiltyp, etc.
  !------------------------------------------------------------
    integer  :: m1, m2, n1, n2, nd ! dimensions of destination grid
    integer  :: i1, i2, id         ! indices    of destination grid
    integer  :: lnn                ! nearest neighbour interpolation flag
    integer  :: j1, j2, jd         ! indices    of source grid
    integer  :: i, j               ! indices
    integer  :: pp,ii,jj,ll        ! indices
    real(wp) :: w, wi              ! weight
    integer  :: counts             ! total number of source gridpoints
    integer  :: misses             ! number of gridpoints not interpolated
    integer  :: found(4)           ! number of gridpoints     interpolated
    integer  :: sts, std           ! source/destination soil type
    real(wp) :: d                  ! distance
    real(wp) :: dmax               ! max. distance for nearest neighbour
    real(wp) :: dcut               ! max. distance allowed
    real(wp) :: dcut2              ! max. distance allowed squared
    real(wp) :: dist2_             ! temporary
    integer  :: idxi(4)            ! temporary for grid indices
    logical  :: reordered          ! same grid but re-ordered points

    lnn = 2
    if (present(nn)) then
      if (nn) lnn = 1
    endif

    !------------------------------------
    ! allocate index array if not present
    !------------------------------------
    m1 = gout% lb (1)
    m2 = gout% lb (2)
    n1 = gout% ub (1)
    n2 = gout% ub (2)
    nd = gout% ub (4)
    if (.not.associated(idx)) then
      allocate (idx(m1:n1,m2:n2,nd))
    endif

    !---------------
    ! derive indices (linear interpolation or simple nearest neighbour)
    !---------------
    if (present (idxin)) then
      idx = idxin
    else
      if (same_horizontal_grid (gin, gout, reordered=reordered)) then
        !------------------------------------
        ! for same grid the indices are known
        !------------------------------------
        if (.not. reordered) then
          do id = 1,nd
            do i2 = m2,n2
              do i1 = m1,n1
                idx(i1,i2,id)% n           = 1
                idx(i1,i2,id)% w    (1)    = 1._wp
                idx(i1,i2,id)% w    (2:)   = 0._wp
                idx(i1,i2,id)% ijdp (1,:)  = [i1,i2,id,dace% pe]
                idx(i1,i2,id)% ijdp (2:,:) = -1
              end do
            end do
          end do
        else  ! re-ordered
          !-------------------------------------------------------------
          ! same grid, but grid-points of the source grid are re-ordered
          !-------------------------------------------------------------
          if (gout% d_gme(1) < 0 .or. gin%  d_gme(1) >= 0)          &
            call finish ('hor_intp_coeff','inconsistent re-ordering')
          !------------------------------------------------------------------
          ! same grid, but grid-points of the destination grid are re-ordered
          !------------------------------------------------------------------
          do id = 1,nd
            do i2 = gout% lbg(2), gout% ubg (2)
              do i1 = gout% lbg(1), gout% ubg (1)
                if ( gout% marr (1,i1,i2,id) /= dace% pe) cycle
                pp = gin % marr (1,i1,i2,id)
                ii = gout% marr (2,i1,i2,id)
                jj = gout% marr (3,i1,i2,id)
                ll = gout% marr (4,i1,i2,id)
                idx(ii,jj,ll)% n           = 1
                idx(ii,jj,ll)% w    (1)    = 1._wp
                idx(ii,jj,ll)% w    (2:)   = 0._wp
                idx(ii,jj,ll)% ijdp (1 ,:) = [i1,i2,id,pp]
                idx(ii,jj,ll)% ijdp (2:,:) = -1
              end do
            end do
          end do
        endif
      else
FTRACE_BEGIN("hor_intp_coeff:grid_indices")
        !------------------------------------------------
        ! different grids: search gridpoints individually
        !------------------------------------------------
        do id = 1,nd
          do i2 = m2,n2
            do i1 = m1,n1
               call grid_indices (gout% rlon(i1,i2,1,id) * r2d, &! <-- geodetic longitude
                                  gout% rlat(i1,i2,1,id) * r2d, &! <-- geodetic latitude
                                  gin,                          &! <-- input grid meta data
                                  idx(i1,i2,id)% ijdp,          &! --> Grid point indices
                                  idx(i1,i2,id)% w,             &! --> Weights
                                  idx(i1,i2,id)% n,             &! number of points returned
                            order=lnn                           )! nearest neighbour interpolation
            end do
          end do
        end do
FTRACE_END  ("hor_intp_coeff:grid_indices")
      endif
    endif

    if (present (stc1)) then
      !--------------------------------------------------------
      ! stc1 present:
      ! search same soil type class within neighbour gridpoints
      !--------------------------------------------------------
      counts = 0
      misses = 0
      found  = 0
      dmax   = 0._wp
      ! Default distance 3 times resolution of courser grid
      dist2_ = 3._wp * max(gin% d_km, gout% d_km)
      if (present(dist2)) dist2_ = dist2
      dcut = dist2_ * 1000._wp / rearth
      dcut2  = dcut * dcut
      !---------------------------------
      ! loop over destination gridpoints
      !---------------------------------
      do id = 1,nd
      do i2 = m2,n2
      do i1 = m1,n1
        counts = counts + 1
        std = sto (i1, i2, id)
        j = 0
        w = huge(0._wp)
        do i = 1, idx(i1,i2,id)%n
          if (.not.present(di)) then
            wi = - idx(i1,i2,id)% w(i)
            if (wi >= w) cycle
          endif
          j1 = idx(i1,i2,id)% ijdp(i,1)
          j2 = idx(i1,i2,id)% ijdp(i,2)
          jd = idx(i1,i2,id)% ijdp(i,3)
          sts = sti (j1, j2, jd)
          if (sts == 9999._wp)       cycle ! missing value in ICON GRIB
          if (std == 9999._wp)       cycle ! missing value in ICON GRIB
          if (.not. stc1 (std, sts)) cycle
          if (present(di)) then
            wi = abs (do(i1,i2,id) - di(j1,j2,jd))
            if (wi == 0._wp) wi = - idx(i1,i2,id)% w(i)
            if (wi >= w) cycle
          endif
          j = i
          w = wi
        end do
        if (j/=0) found(1) = found(1) + 1
        if (j==0 .and. present (stc2)) then
          !------------------------------------------------------
          ! stc2 present:
          ! try other soil type class within neighbour gridpoints
          !------------------------------------------------------
          do i = 1, idx(i1,i2,id)%n
            if (.not.present(di)) then
              wi = - idx(i1,i2,id)% w(i)
              if (wi >= w) cycle
            endif
            j1 = idx(i1,i2,id)% ijdp(i,1)
            j2 = idx(i1,i2,id)% ijdp(i,2)
            jd = idx(i1,i2,id)% ijdp(i,3)
            sts = sti (j1, j2, jd)
            if (sts == 9999._wp)       cycle ! missing value in ICON GRIB
            if (std == 9999._wp)       cycle ! missing value in ICON GRIB
            if (.not. stc2 (std, sts)) cycle
            if (present(di)) then
              wi = abs (do(i1,i2,id) - di(j1,j2,jd))
              if (wi == 0._wp) wi = - idx(i1,i2,id)% w(i)
              if (wi >= w) cycle
            endif
            j = i
            w = wi
          end do
          if (j/=0) found(2) = found(2) + 1
        endif
        if (j/=0) then
          idx(i1,i2,id)% n         = 1
          idx(i1,i2,id)% w         = 0._wp
          idx(i1,i2,id)% w   (1)   = 1._wp
          idx(i1,i2,id)% ijdp(1,:) = idx(i1,i2,id)% ijdp(j,:)
        else if (present (stc2)) then
          !---------------------------------------------------
          ! stc2 present:
          ! search nearest neighbour for other soil type class
          !---------------------------------------------------
          w      = 9._wp
          idxi   = -1
          do jd = gin% lbg(4), gin% ubg(4)
          do j2 = gin% lbg(2), gin% ubg(2)
          do j1 = gin% lbg(1), gin% ubg(1)
            sts = sti (j1, j2, jd)
            if (sts == 9999._wp)       cycle ! missing value in ICON GRIB
            if (std == 9999._wp)       cycle ! missing value in ICON GRIB
            if (.not. stc2 (std, sts)) cycle

            if ((gout% xnglob(i1,i2,1,id) - gin% xnglob(j1,j2,1,jd))**2 > dcut2) cycle
            if ((gout% xnglob(i1,i2,2,id) - gin% xnglob(j1,j2,2,jd))**2 > dcut2) cycle
            if ((gout% xnglob(i1,i2,3,id) - gin% xnglob(j1,j2,3,jd))**2 > dcut2) cycle
            d = sqrt (sum ((gout% xnglob(i1,i2,:,id) - gin% xnglob(j1,j2,:,jd)) **2))
            if (d > dcut)                cycle
            if (d < w) then
              j         = 1
              w         = d
              idxi(1:3) = gin% marr (2:4,j1,j2,jd)
              idxi(  4) = gin% marr (  1,j1,j2,jd)
            endif
          end do
          end do
          end do
          if (j > 0) then
            found(3)                 = found(3) + 1
            idx(i1,i2,id)% n         = 1
            idx(i1,i2,id)% w         = 0._wp
            idx(i1,i2,id)% w   (1)   = 1._wp
            idx(i1,i2,id)% ijdp(1,:) = idxi
            dmax = max (dmax, w)
          endif
        endif
        if (j==0 .and. present (stc2) .and. .not.present(dist2) ) then
          !------------------------------------
          ! stc2 present and dist2 not present:
          ! take nearest neighbour
          !------------------------------------
          do i = 1, idx(i1,i2,id)%n
            j1 = idx(i1,i2,id)% ijdp(i,1)
            j2 = idx(i1,i2,id)% ijdp(i,2)
            jd = idx(i1,i2,id)% ijdp(i,3)
            sts = sti (j1, j2, jd)
            if (sts == 9999._wp)       cycle ! missing value in ICON GRIB
            if (std == 9999._wp)       cycle ! missing value in ICON GRIB
            wi = - idx(i1,i2,id)% w(i)
            if (wi >= w) cycle
            j = i
            w = wi
          end do
          if (j/=0) then
            idx(i1,i2,id)% n         = 1
            idx(i1,i2,id)% w         = 0._wp
            idx(i1,i2,id)% w   (1)   = 1._wp
            idx(i1,i2,id)% ijdp(1,:) = idx(i1,i2,id)% ijdp(j,:)
            found(4)                 = found(4) + 1
          endif
        endif
        if (j==0) then
          misses                = misses + 1
          idx(i1,i2,id)% n      = 0
          idx(i1,i2,id)% w      = 0._wp
        endif
      end do
      end do
      end do
      dmax = dmax * rearth / 1000._wp
      dmax   = p_max (dmax)
      counts = p_sum (counts)
      found  = p_sum (found)
      misses = p_sum (misses)
      if (dace% lpio) write(6,'(a,6i10,f10.3)') &
        'hor_intp_coeff: counts found misses dmax=',counts,found,misses,dmax
    endif

  end subroutine hor_intp_coeff
!------------------------------------------------------------------------------
  subroutine hor_intp_field (xout, xin, gout, gin, idx)
  type (t_grid)  ,intent(in)        :: gout
  type (t_grid)  ,intent(in)        :: gin
  real(wp)       ,intent(out)       :: xout (gout%lb (1):,gout%lb (2):,:,:)
  real(wp)       ,intent(in)        :: xin  (gin %lbg(1):,gin %lbg(2):,:,:)
  type (t_h_idx) ,pointer ,optional :: idx (:,:,:)
  !------------------------------------------
  !  horizontally interpolate a field
  !
  !! currently not implemented in parallel !!
  !------------------------------------------
    type (t_h_idx) ,pointer :: ix (:,:,:)
    integer :: m1, m2, n1, n2, n3, nd
    integer :: i1, i2, i3, id, i
    logical :: first
#if defined (__SX__)
    integer :: n3_, l1, u1, l2, u2, u3, u4
    real(wp), allocatable :: xout_(:,:,:,:), xin_(:,:,:,:)
#endif
    !------------------------------------
    ! allocate index array if not present
    !------------------------------------
    m1 = lbound(xout,1)
    m2 = lbound(xout,2)
    n1 = ubound(xout,1)
    n2 = ubound(xout,2)
    n3 = ubound(xout,3)
    nd = ubound(xout,4)
    if (present (idx)) then
      if (.not.associated(idx)) then
        first = .true.
      else
        first = .false.
      endif
      ix => idx
    else
      first = .true.
    endif
    !---------------
    ! derive indices
    !---------------
FTRACE_BEGIN("hor_intp_field:grid_indices")
    if (first) call hor_intp_coeff (gout, gin, ix)
FTRACE_END  ("hor_intp_field:grid_indices")
    !------------
    ! interpolate
    !------------
FTRACE_BEGIN("hor_intp_field:interpolate")
#if defined (__SX__)
    !-----------------------------------------------------------------------
    ! Attempt vectorization over vertical levels.
    ! Use auxiliary fields with odd first dimension to reduce bank conflicts
    !-----------------------------------------------------------------------
    n3_ = n3
    if (n3 > 1 .and. mod (n3,2) == 0) n3_ = n3+1
!print *, "hor_intp_field: n3_=", n3_

    l1 = lbound(xin,1)
    l2 = lbound(xin,2)
    u1 = ubound(xin,1)
    u2 = ubound(xin,2)
    u3 = ubound(xin,3)
    u4 = ubound(xin,4)
    allocate (xin_(n3_,l1:u1,l2:u2,u4))
    do i3 = 1, n3
       xin_(i3,:,:,:) = xin (:,:,i3,1:u4)
    end do

    allocate (xout_(n3_,m1:n1,m2:n2,nd))
    xout_(:,:,:,:) = 0._wp
    do id = 1,nd
      do i2 = m2,n2
        do i1 = m1,n1
          do i3 = 1,n3
!!!CDIR EXPAND=8
            do i=1,ix(i1,i2,id)%n
              xout_(i3,i1,i2,id) = xout_(i3,i1,i2,id)  &
                                 + ix(i1,i2,id)%w(i) * &
               xin_(i3,                                &
                    ix(i1,i2,id)% ijdp(i,1),           &
                    ix(i1,i2,id)% ijdp(i,2),           &
                    ix(i1,i2,id)% ijdp(i,3))
            end do
          end do
        end do
      end do
    end do
    do i3 = 1, n3
       xout (:,:,i3,:) = xout_(i3,:,:,:)
    end do
    deallocate (xin_, xout_)
#else
    do id = 1,nd
      do i2 = m2,n2
        do i1 = m1,n1
          do i3 = 1,n3
            xout (i1,i2,i3,id) = 0._wp
            do i=1,ix(i1,i2,id)%n
              xout (i1,i2,i3,id) = xout (i1,i2,i3,id)  &
                                 + ix(i1,i2,id)%w(i) * &
               xin (ix(i1,i2,id)% ijdp(i,1),           &
                    ix(i1,i2,id)% ijdp(i,2), i3,       &
                    ix(i1,i2,id)% ijdp(i,3))
            end do
          end do
        end do
      end do
    end do
#endif
FTRACE_END  ("hor_intp_field:interpolate")
    !--------------------------------------
    ! deallocate index array if not present
    !--------------------------------------
    if (present (idx)) then
      if (.not.associated(idx)) deallocate (ix)
    else
      deallocate (ix)
    endif
  end subroutine hor_intp_field
!==============================================================================
  subroutine fix_grid_indices (ix, grid, d)
  type (t_h_idx) ,intent(inout)     :: ix    ! interpolation coefficients
  type (t_grid)  ,intent(in)        :: grid  ! grid meta data
  integer        ,intent(in)        :: d     ! required diamond index
  !-----------------------------------------------
  ! fix interpolation coefficients to GME grid:
  ! (used for LETKF on coarse grid)
  ! remove interpolation indices with zero weights
  ! (nearly zero weights due to rounding errors)
  ! fix interpolation indices at diamond bounds
  ! (so that all references are to diamond 'd'
  !-----------------------------------------------

    integer             :: i, j
    real(wp) ,parameter :: eps = 1.e-10 ! 'small' weight
    !-----------------------------------------------
    ! remove interpolation indices with zero weights
    ! (nearly zero weights due to rounding errors)
    !-----------------------------------------------
    j = 0
    do i = 1, ix% n
      if (abs (ix% w(i)) > eps) then
        j = j + 1
        ix% w   (j)   = ix% w   (i)
        ix% ijdp(j,:) = ix% ijdp(i,:)
      endif
    end do
    ix% n             = j
    ix% w    (j+1:)   = 0._wp
    ix% ijdp (j+1:,:) = 0
    !--------------------------------------------
    ! fix interpolation indices at diamond bounds
    ! (so that all references are to diamond 1
    !--------------------------------------------
    do i = 1, ix% n
      !----------------------------------
      ! index already refers to diamond 1
      !----------------------------------
      if      (ix% ijdp (i,3) == d) cycle
      !----------------------------------
      ! the remainder only refers to d==1
      !----------------------------------
!     if (d == 1) then
!       !------------------------------------------------------
!       ! diamond 2: i2 = lbg(1) = 1: exchange i1, i2, offset 1
!       !------------------------------------------------------
!       if (ix% i (i,3) == 2 .and. ix% i (i,2) == grid% lbg(2)) then
!         ix% i (i,2) = ix% i (i,1) + 1
!         ix% i (i,1) = grid% lbg(1)
!         ix% i (i,3) = 1
!         cycle
!       endif
!       !--------------------------------------------------------------
!       ! diamond 6: i1 = ubg(1) = 1: exchange i1, i2, reverse counting
!       !--------------------------------------------------------------
!       if (ix% i (i,3) == 6 .and. ix% i (i,1) == grid% ubg(1)) then
!         ix% i (i,1) = grid% ubg(2) - ix% i (i,2)
!         ix% i (i,2) = grid% ubg(2)
!         ix% i (i,3) = 1
!         cycle
!       endif
!     endif
      !--------------------------------------
      ! other diamond, this should not happen
      !--------------------------------------
      write (0,*) 'fix_coarse_grid_indices: invalid diamond index'
      do j = 1, ix% n
        write (0,*) dace% pe,'d,i,j,d,w,bounds(i),bounds(j):',        &
                    d, ix% ijdp(j,1:3), ix% w(j),                     &
                    grid% lbg(1),grid% ubg(1),grid% lbg(2),grid% ubg(2)
      end do
      call finish('fix_coarse_grid_indices','invalid diamond index')
    end do

  end subroutine fix_grid_indices
!==============================================================================
  subroutine distribute_icon_gme (g_icon, g_gme)
  !--------------------------------------------------------------------------
  ! Index the ICON gridpoints so that they are distributed according to the
  ! partitioning of a given GME (or regular/rotated lat-lon) reference grid.
  ! ICON gridpoints 'x' are assigned to the PE holding the reference
  ! gridpoints 'Oi,j' with min(j) and min(i) taken over all 3 neighbours 'o'.
  ! The diamond is not changed, even if the gridpoint i,j resides on
  ! another diamond (which may be the case for the GME grid at the poles).
  !
  ! The permutation is returned by setting the array g_icon% marr.
  !
  ! Example:
  !            Oi,j-----oi+1,j
  !             |      / |
  !             |     /  |
  !             |    /   |
  !             |   / x  |
  !             |  /     |
  !             | /      |
  !            oi,j+1---oi+1,j+1
  !
  !--------------------------------------------------------------------------
  type (t_grid) ,intent(inout) :: g_icon ! ICON grid meta data
  type (t_grid) ,intent(in)    :: g_gme  ! GME  grid meta data

    type (t_h_idx) ,allocatable :: ix (:,:,:)     ! indices for interpolation
    integer                     :: lbi(4), ubi(4) ! local ICON grid bounds
!   integer                     :: lbg(4), ubg(4) ! local GME  grid bounds
    integer                     :: lb1, ub1       ! old bounds
    integer                     :: i,j,d,k        ! grid indices
    integer                     :: i1, i2
    integer                     :: n              ! total # of ICON gridpoints
    integer                     :: pe             ! processor index
    integer        ,allocatable :: ir (:,:,:,:)   ! reference indices
    integer                     :: ii    (4)
    integer        ,allocatable :: ix_pe (:)
    integer        ,allocatable :: ix_d  (:)
    integer                     :: ip (  -1   :dace% npe-1) ! # gridpoints per pe
    integer                     :: id (0:10, 0:dace% npe-1) ! # " per pe/diamond
!   integer        ,allocatable :: ic (:,:,:)     ! counter for testing

    select case (g_gme% gridtype)
    case (WMO6_LATLON, WMO6_GAUSSIAN, WMO6_ROTLL, DWD6_ICOSAHEDRON)
    case default
       write(*,*) 'distribute_icon_gme: unsupported reference grid', &
            g_gme% gridtype
       call finish ('distribute_icon_gme','unsupported reference grid')
    end select
    !-------------------------------
    ! set up destination index array
    !-------------------------------
    lbi = g_icon% lb
    ubi = g_icon% ub
!   lbg = g_gme % lbg
!   ubg = g_gme % ubg - (/1,1,0,0/)  !upper diamod edges go to another PE

    allocate (ix (   lbi(1):ubi(1), lbi(2):ubi(2), lbi(4):ubi(4)))
    allocate (ir (4, lbi(1):ubi(1), lbi(2):ubi(2), lbi(4):ubi(4)))

!   !-------------
!   ! for testing:
!   !-------------
!   allocate (ic (   lbg(1):ubg(1), lbg(2):ubg(2), lbg(4):ubg(4)))
!   ic = 0

    id = 0
    ip = 0
    do d =     lbi(4),ubi(4)
      do j =   lbi(2),ubi(2)
        do i = lbi(1),ubi(1)
          !------------------------------------------------
          ! get interpolation indices to GME/reference grid
          !------------------------------------------------
          call grid_indices (r2d * g_icon% rlon(i,j,1,d),     &! + 36._wp,&
                             r2d * g_icon% rlat(i,j,1,d),     &
                                   g_gme,    ix(i,j,d)% ijdp, &
                                             ix(i,j,d)% w,    &
                                             ix(i,j,d)% n,    &
                             lunique=.false.                  )

          !-------------------------------------
          ! get reference gridpoint (lowest i,j)
          ! should not leave the diamond
          !-------------------------------------
          ii = ix(i,j,d)% ijdp (1,:)

          if (ix(i,j,d)% ijdp(2,3) /= ii(3) .or. &
              ix(i,j,d)% ijdp(3,3) /= ii(3)    ) &
            call finish('distribute_icon_gme','diamond mismatch')
          do k = 2, ix (i,j,d)% n
            if (ix(i,j,d)% ijdp(k,1)  < ii(1)) then
              ii(1) = ix(i,j,d)% ijdp(k,1)
              ii(3) = ix(i,j,d)% ijdp(k,3)
            endif
            if (ix(i,j,d)% ijdp(k,2)  < ii(2)) then
              ii(2) = ix(i,j,d)% ijdp(k,2)
              ii(3) = ix(i,j,d)% ijdp(k,3)
            endif
          end do
          !--------------------
          ! get processor index
          !--------------------
          select case (g_gme% gridtype)
          case default
            ii (4) = g_gme % marr(1, ii(1), ii(2), ii(3))
          case (DWD6_ICOSAHEDRON)
            i1 = (ii(1) - 0) * nproc1 / g_gme% ni
            i2 = (ii(2) - 1) * nproc2 / g_gme% ni
            ii (4) = i1 + i2 * nproc1
          end select
          ir (:, i, j, d) = ii
          !-------
          ! counts
          !-------
!         ic (ii(1), ii(2), ii(3)) = ic (ii(1), ii(2), ii(3)) + 1
          id (ii(3), ii(4))        = id (ii(3), ii(4))        + 1
          ip (ii(4))               = ip (ii(4))               + 1
        end do
      end do
    end do
    id = p_sum (id)
    ip = p_sum (ip)

!   !-----------------------------
!   ! print out counts for testing
!   !-----------------------------
!   ic = p_sum (ic)
!   if (dace% pe == 0) then
!     print *,'### distribute_icon_gme: ic =',minval(ic),maxval(ic)
!     print *,'### distribute_icon_gme: ic =',minloc(ic)
!     print *,'### distribute_icon_gme: ic =',maxloc(ic)
!     print *,'### distribute_icon_gme: ic =',count (ic/=2)
!     do d=1,1
!!    write(6,'(a,2i3,a,99i3)')'### distribute_icon_gme: ic =',0,0,' - ',(/(mod(i,10),i=lbg(1),ubg(1))/)
!     do j=lbg(2), ubg(2)
!       write(6,'(a,2i3,a,99i3)')'### distribute_icon_gme: ic =',d,j,' - ',ic(lbg(1):ubg(1),j,d)
!     end do
!     end do
!!    print *,'### distribute_icon_gme: id =',id
!     print *,'### distribute_icon_gme: ip =',ip
!   endif

    !---------
    ! printout
    !---------
    if(dace% lpio) then
      write(6,*)
      write(6,*) 'distribute_icon_gme: number of ICON gridpoints per PE'
      do i=1, nproc2
      k = (i-1)*nproc1
        write(6,*) i, k, '..', min(k+nproc1-1,dace% npe-1),' : ', &
                   (ip(j), j=k,min(k+nproc1-1,dace% npe-1))
      end do
      write(6,*)
      select case (g_gme% gridtype)
      case (DWD6_ICOSAHEDRON)
        write(6,*) 'distribute_icon_gme: number of GME gridpoints per PE'
      case default
        write(6,*) 'distribute_icon_gme: number of points / PE, reference grid'
      end select
      do i=1, nproc2
      k = (i-1)*nproc1
        write(6,*) i, k, '..', min(k+nproc1-1,dace% npe-1),' : ', &
                  ((g_gme% dc% ilim2(i)-g_gme% dc% ilim2(i-1))   &
                  *(g_gme% dc% ilim1(j)-g_gme% dc% ilim1(j-1)),j=1,nproc1)
      end do
      write(6,*)
    end if

    !-----------------------
    ! set up new index array
    !-----------------------
    n = g_icon% nxny
    if (associated (g_icon% marr)) deallocate (g_icon% marr)
    allocate (g_icon% marr (4,n,1,1))
    allocate (        ix_pe  (n)    )
    allocate (        ix_d   (n)    )

    do i=1, dace% npe-1
      ip (i)   = ip (i)   + ip (i-1)
    end do
    do i=2,10
      id (i,:) = id (i,:) + id (i-1,:)
    end do

    call p_allgather (ir (3,:,1,1), ix_d )
    call p_allgather (ir (4,:,1,1), ix_pe)

    !------------------------------------------
    ! adjust grid% d_gme (local diamond offset)
    !------------------------------------------
    g_icon% d_gme = id (:, dace% pe)

    !-----------------------------------
    ! adjust grid% lb, ub (local bounds)
    !-----------------------------------
    lb1               = g_icon% lb(1)
    ub1               = g_icon% ub(1)
    g_icon% lb   (1)  = ip (dace% pe - 1) + 1
    g_icon% ub   (1)  = ip (dace% pe    )
    g_icon% shape(1)  = ip (dace% pe    ) - ip (dace% pe - 1)
    g_icon% dc% ilim1 = ip + 1
    where (g_icon%m%i% lb(1) == lb1 .and. g_icon%m%i% ub(1) == ub1)
      g_icon%m%i% lb(1) = g_icon% lb(1)
      g_icon%m%i% ub(1) = g_icon% ub(1)
    endwhere

    !----------------------------------------
    ! adjust grid% marr (permutation indices)
    !----------------------------------------
    do i=1,n
      pe = ix_pe(i)
      d  = ix_d (i)
      j  = id (d-1, pe  ) + 1
      id      (d-1, pe  ) = j
      j  = ip (     pe-1) + j
      g_icon% marr (1,i,1,1) = pe
      g_icon% marr (2,i,1,1) = j
      g_icon% marr (3,i,1,1) = 1
      g_icon% marr (4,i,1,1) = 1
    end do

    !-------------------
    ! Consistency checks
    !-------------------
    select case (g_gme% gridtype)
    case default
       do pe = 0, dace% npe-1
          if (any (id(0,pe) /= id(1:,pe))) then
             write(0,*) dace% pe, ': distribute_icon_gme: pe, id(:,pe)=', &
                  pe, id(:,pe)
             call finish ("distribute_icon_gme", "internal error")
          end if
       end do
    case (DWD6_ICOSAHEDRON)
    end select

    !--------------------------------------------------
    ! adjust grid% rlon, rlat (coordinates already set)
    !--------------------------------------------------
    g_icon% rlat   (g_icon% marr (2,:,1,1),1,1,1) = g_icon% rlat   (:,1,1,1)
    g_icon% rlon   (g_icon% marr (2,:,1,1),1,1,1) = g_icon% rlon   (:,1,1,1)
    g_icon% xnglob (g_icon% marr (2,:,1,1),1,:,1) = g_icon% xnglob (:,1,:,1)
    g_icon% dlat   (g_icon% marr (2,:,1,1)      ) = g_icon% dlat   (:)
    g_icon% dlon   (g_icon% marr (2,:,1,1)      ) = g_icon% dlon   (:)

  end subroutine distribute_icon_gme
!==============================================================================
  function mode2text (mode) result (text)
    !---------------------------------------
    ! Describe horizontal interpolation mode
    !---------------------------------------
    integer, intent(in) :: mode
    integer, parameter  :: ltext = size (horint_modes)     &
                                 * len  (horint_modes% desc)
    character(ltext)    :: text

    integer :: i

    if (mode == 0) then
       text = horint_modes(1)% desc
    else
       text = ""
       do i = 1, size (horint_modes)
          if (iand (mode, horint_modes(i)% mode) /= 0) &
               text = trim (text) // " " // horint_modes(i)% desc
       end do
       text = adjustl (text)
    end if
  end function mode2text
  !----------------------------------------------------------------------------
  function text2mode (text) result (mode)
    !---------------------------------------------------
    ! Decode horizontal interpolation mode specification
    !---------------------------------------------------
    character(*), intent(in) :: text
    integer                  :: mode

    character(len(text)) :: s
    character(12)        :: desc(10)
    integer              :: i, j, n
    logical              :: ok

    ! Convert to lowercase and replace commas by blanks
    s = tolower (text)
    forall (i=1:len(s)) s(i:i) = merge (" ",s(i:i),s(i:i)==",")

    call split (desc, s, n)
    if (n < 0) call finish ("text2mode","too many criteria")

    mode = 0
    do i = 1, n
       if (desc(i) == "") cycle
       ok = .false.
       do j = 1, size (horint_modes)
          if (desc(i) == horint_modes(j)% desc) then
             mode = ior (horint_modes(j)% mode, mode)
             ok   = .true.
             exit
          end if
       end do
       if (.not. ok) call finish ("text2mode","bad string: " // trim (desc(i)))
    end do

    if (mode > 0 .and. iand (mode, HORINT_NEAREST + HORINT_ZDIST + &
                                   HORINT_FRLAND  + HORINT_FRSEA   ) == 0) &
         mode = ior (mode, HORINT_DEFAULT)

    if (.not. check_horint (mode)) then
       write(0,*) "text2mode: bad mode = ", trim (text)
       call finish ("text2mode","conflicting mode specifications")
    end if

  end function text2mode
  !----------------------------------------------------------------------------
  function check_horint (mode) result (res)
    !-----------------------------------
    ! Check for valid mode specification
    !-----------------------------------
    integer, intent(in)  :: mode
    logical              :: res

    res = .false.
    if (mode < 0 .or. mode >= HORINT_OTHER) return
    !-----------------------------------
    ! Cannot have both NEAREST and ZDIST
    ! or largest land or sea fraction.
    ! Cannot have conflicting choices
    !-----------------------------------
    if (count (iand (mode, [HORINT_DEFAULT,  &
                            HORINT_NEAREST,  &
                            HORINT_ZDIST,    &
                            HORINT_FRLAND,   &
                            HORINT_FRSEA  ]  ) /= 0) > 1) return
    !-----------------------------------
    ! Cannot have multiple surface types
    !-----------------------------------
    if (count (iand (mode, [HORINT_LAND,     &
                            HORINT_SEA,      &
                            HORINT_LANDONLY, &
                            HORINT_SEAONLY ] ) /= 0) > 1) return

    res = .true.
  end function check_horint
  !----------------------------------------------------------------------------
  subroutine filter_idx (hic, mc, grid, mode, href, np)
    !----------------------------------------------------
    ! Select gridpoints for horizontal interpolation mode
    !----------------------------------------------------
    type(t_hic)   ,intent(inout) :: hic      ! hor. interpolation coefficients
    type(t_mcols) ,intent(in)    :: mc       ! model column descriptor
    type(t_grid),  intent(in)    :: grid     ! model grid description
    integer,       intent(in)    :: mode     ! Gridpoint selection mode
    real(wp),      intent(in)    :: href     ! Reference height [m]
    integer,       intent(out)   :: np       ! Number of accepted gridpoints

    logical               :: land, landonly !
    logical               :: sea,  seaonly  !
    logical               :: nn, nz         ! nearest gridpoint method
    logical               :: fl, fs         ! maximum land/sea fraction method
    integer               :: ng             ! Number of provided gridpoints
    integer               :: i              ! Loop index
    integer               :: ic, ix, i0     ! Indices
    integer               :: inn            ! Aux. index
    integer               :: ijd (3)        ! (i,j,d)
    integer, parameter    :: mg = 16
    real(wp)              :: lsm (mg)       ! land-sea mask
    real(wp)              :: h   (mg)       ! model surface height
    logical               :: mask(mg)       ! gridpoint selection mask
    real(wp)              :: sumw           ! temporary (sum of weights)
    real(wp)              :: mf             ! temporary (min/max land fraction)
    integer               :: m(4)           ! Auxiliary bitmask
    type(t_mcol), pointer :: c
    logical               :: hs

    hs = associated (grid% hsurf)
    ng = 0
    do i = 1, size (hic% imc, dim=1)
       ic = hic% imc(i,1)
       if (ic == 0) exit
       ng =  i
       c  => mc% c(ic)
       ijd(:) = c% ijdtp(1:3)
       lsm(i) = grid% lsm  (ijd(1),ijd(2),1,ijd(3))
       if (hs) then
         h(i) = grid% hsurf(ijd(1),ijd(2),1,ijd(3))
       else
         h(i) = grid% geosp(ijd(1),ijd(2),1,ijd(3)) / gacc
       end if
    end do
    np = ng

    ! Nothing to process?
    if (np   == 0) return

    ! Default mode: horizontal interpolation, no surface type selection
    if (mode == 0) return

    m = iand (mode, [HORINT_LAND,     HORINT_SEA,   &
                     HORINT_LANDONLY, HORINT_SEAONLY])
    land     = m(1) /= 0
    sea      = m(2) /= 0
    landonly = m(3) /= 0
    seaonly  = m(4) /= 0

    if     (land .or. landonly) then
       mask(:ng) = lsm(:ng) >= 0.5
       if (land) then
          if (count (mask(:ng)) == 0) then
             mask(:ng) = .true.  ! Fallback: use all seapoints
             land      = .false.
          end if
       end if
    else if (sea .or. seaonly) then
       mask(:ng) = lsm(:ng) <  0.5
       if (sea) then
          if (count (mask(:ng)) == 0) then
             mask(:ng) = .true.  ! Fallback: use all landpoints
             sea       = .false.
          end if
       end if
    else
       mask(:ng) = .true.
    end if
!   mask(ng+1:)  = .false.

    ! Check if there are any points left to choose from
    np = count (mask(:ng))
    if (np == 0) then
       hic% imc(1,1) = 0
       return
    end if

    ! Remember eligible gridpoint with largest weight
    i0 = maxloc (hic% w(1:ng), mask=mask(:ng), dim=1)

    ix = -1
    nn = iand (mode, HORINT_NEAREST) /= 0
    nz = iand (mode, HORINT_ZDIST)   /= 0
    fl = iand (mode, HORINT_FRLAND)  /= 0
    fs = iand (mode, HORINT_FRSEA)   /= 0
    if (nn .or. nz .or. fl .or. fs) then
       !---------------------
       ! Use only 1 gridpoint
       !---------------------
       if (nn) then
          ic = i0
       else if (nz) then
          ic = minloc (abs (h(1:ng)-href), mask=mask(:ng), dim=1)
          !write(0,*) "### hsurf :", real (h(1:ng)), ":", ic, i0
       else if (fl) then
          ! Use nearest point with maximum land fraction
          mf = maxval (lsm   (1:ng), mask=mask(:ng))
          ic = maxloc (hic% w(1:ng), mask=mask(:ng) .and. (lsm(:ng)==mf), dim=1)
          !write(0,*) "### frland:", real (lsm(1:ng)), ":", ic, i0
       else !if (fs) then
          ! Use nearest point with minimum land fraction
          mf = minval (lsm   (1:ng), mask=mask(:ng))
          ic = maxloc (hic% w(1:ng), mask=mask(:ng) .and. (lsm(:ng)==mf), dim=1)
          !write(0,*) "### frsea :", real (lsm(1:ng)), ":", ic, i0
       end if
       hic% imc(1,:)  = hic% imc(ic,:)
       hic% imc(2:,:) = 0
       hic% w  (1)    = 1.
       inn            = 1
       ix             = hic% imc(inn,1)
       np             = 1
    else if (np == ng) then
       !-----------------------
       ! Retain all gridpoints?
       !-----------------------
       return
    else
       !----------------------------
       ! Use all eligible gridpoints
       !----------------------------
       sumw = sum (hic% w(1:ng), mask=mask(:ng))
       if (sumw > 0.) then
          ic = 0
          do i = 1, ng
             if (mask(i)) then
                ic = ic + 1
                if (ic < i) then
                   hic% imc(ic,:) = hic% imc(i,:)
                   hic% w  (ic)   = hic% w  (i)
                end if
             end if
          end do
          hic% imc(np+1:,:) = 0
          hic% w  (np+1:)   = 0.
          !------------
          ! Renormalize
          !------------
          hic% w(1:np) = hic% w(1:np) / sumw
          inn = maxloc (hic% w(1:np), dim=1)
          ix  = hic% imc (inn,1)
       else
          np = 0
       end if
    end if ! (nn .or. nz)

    if (ix > 0) then
       hic% ijdp(1:4) = mc% c(ix)% ijdtp([1,2,3,5])
    end if

  end subroutine filter_idx
!==============================================================================
end module mo_grid_intpol
