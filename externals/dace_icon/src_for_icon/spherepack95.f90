!
!+ Fortran95 Interfaces for SPHEREPACK 3
!
! $Id$
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NOTE: This file does require preprocessing!  See below for details.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE spherepack95
!
! Description:
!   Interface definitions and Fortran90/95 style wrapper routines
!   to selected routines from UCAR's SPHEREPACK for more convenient
!   and safer use.  The wrapper routines support only regular grids
!   (also Gaussian grid for u,v -> streamfct.,vel.pot. (sp95_sfvp_gg)
!   with the conventions for geophysical coordinates used at DWD.
!
!   Requires a prebuilt version of UCAR's SPHEREPACK 3
!   (http://www.scd.ucar.edu/css/software/spherepack)
!
! Current Code Owner: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Comment
! ------------ ---------- -------
!   1.0    2005-02-11  Initial version.                      (Harald Anlauf)
!          2006-04-19  sp95_sfvp_gg added for Gaussian grid  (Andreas Rhodin)
!          2007-07-03  sp95_grad_gg: gradient on Gaussian grid [ha]
!          2007-09-24  sp95_sha_gg:  scalar analysis,
!                      sp95_shs_gg:  scalar synthesis on Gaussian grid [ha]
!          2007-11-13  sp95_shp_gg:  scalar harmonic projection,
!                      sp95_slap_gg: scalar Laplacian,
!                      sp95_vrtdiv_sfvp_gg: sf,vp,u,v from vort.,div. [ha]
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_8         2009/12/09 Harald Anlauf
!  Add forgotten variable declaration
! V1_9         2010/04/20 Harald Anlauf
!  sp95_sha, sp95_shs: spectral analysis and synthesis of scalar fields on regular grid
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Author:              Harald Anlauf <harald.anlauf@dwd.de>
!------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Define here the actual working precision of the spherepack library
! (This kludge is required because the IMPORT statement for interfaces
! was not introduced before Fortran 2003.)
!
!!! Single precision version:
!#define _wp_  kind(1.0)
!!! Double precision version:
#define _wp_  kind(1.d0)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  use mo_kind,      only: wp
  use mo_exception, only: finish
  implicit none

!!! Global Declarations:

  public :: sp95_init, sp95_cleanup
  public :: sp95_vrtdiv_s, sp95_vrtdiv
  public :: sp95_sfvp_s !u,v->streamf.,vel.p.,regular grid,runtime optimization
  public :: sp95_sfvp   !u,v->streamf.,vel.p.,regular grid,memory  optimization
  public :: sp95_sfvp_gg!u,v->streamf.,vel.p.,Gaussiangrid,memory  optimization
  public :: sp95_grad_gg        ! gradient of scalar field, Gaussian grid (2d)
  public :: sp95_sha            ! scalar analysis,  regular grid (2d)
  public :: sp95_shs            ! scalar synthesis, regular grid (2d)
  public :: sp95_sha_gg         ! scalar analysis,  Gaussian grid (2d)
  public :: sp95_shs_gg         ! scalar synthesis, Gaussian grid (2d)
  public :: sp95_shp_gg         ! harmonic projection of scalar field, Gauss.
  public :: sp95_slap_gg        ! Laplacian of a scalar field, Gaussian grid
  public :: sp95_sfvp_vrtdiv_gg ! sf,vp,u,v from vorticity, divergence; Gauss.

  private

!  integer, parameter :: sp = kind (1.0)
!  integer, parameter :: dp = kind (1.d0)
!  integer, parameter :: wp = _wp_

  real(wp), parameter :: PI = 3.1415926535897932_wp

  ! Radius of sphere
  real(wp) :: radius = 0

  ! Grid properties
  integer  :: nlat = 0
  integer  :: nlon = 0

  ! Working storage for the spherepack routines
  real(wp), allocatable :: wshaec(:)
  real(wp), allocatable :: wshags(:)
  real(wp), allocatable :: wshagc(:)    ! Scalar Harmonic Analysis
  real(wp), allocatable :: wshses(:)
  real(wp), allocatable :: wshsec(:)
  real(wp), allocatable :: wshsgc(:)    ! Scalar Harmonic Synthesis
  real(wp), allocatable :: wvhaes(:)
  real(wp), allocatable :: wvhaec(:)
  real(wp), allocatable :: wvhagc(:)    ! Vector Harmonic Analysis
  real(wp), allocatable :: wvhsgs(:)
  real(wp), allocatable :: wvhsgc(:)    ! Vector Harmonic Synthesis

!!! Interfaces to "bare" spherepack routines

  interface
     subroutine shaec (nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                       wshaec, lshaec, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, isym, nt, idg, jdg, mdab, ndab, &
                               lshaec, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: g(idg,jdg,*), a(mdab,ndab,*), b(mdab,ndab,*), &
                               wshaec(lshaec), work(lwork)
     end subroutine shaec
  end interface

  interface
     subroutine shagc (nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                       wshagc, lshagc, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, isym, nt, idg, jdg, mdab, ndab, &
                               lshagc, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: g(idg,jdg,*), a(mdab,ndab,*), b(mdab,ndab,*), &
                               wshagc(lshagc), work(lwork)
     end subroutine shagc
  end interface

  interface
     subroutine shags (nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                       wshags, lshags, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, isym, nt, idg, jdg, mdab, ndab, &
                               lshags, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: g(idg,jdg,*), a(mdab,ndab,*), b(mdab,ndab,*), &
                               wshags(lshags), work(lwork)
     end subroutine shags
  end interface

  interface
     subroutine shaeci (nlat, nlon, wshaec, lshaec, dwork, ldwork, ierror)
       integer, intent(in)  :: nlat, nlon, lshaec, ldwork
       integer, intent(out) :: ierror
       real(_wp_)           :: wshaec(lshaec)
       double precision     :: dwork(ldwork)
     end subroutine shaeci
  end interface

  interface
     subroutine shagci (nlat, nlon, wshagc, lshagc, dwork, ldwork, ierror)
       integer, intent(in)  :: nlat, nlon, lshagc, ldwork
       integer, intent(out) :: ierror
       real(_wp_)           :: wshagc(lshagc)
       double precision     :: dwork(ldwork)
     end subroutine shagci
  end interface

  interface
     subroutine shagsi (nlat, nlon, wshags, lshags, work, lwork, &
                        dwork, ldwork, ierror)
       integer, intent(in)  :: nlat, nlon, lshags, lwork, ldwork
       integer, intent(out) :: ierror
       real(_wp_)           :: wshags(lshags), work(lwork)
       double precision     :: dwork(ldwork)
     end subroutine shagsi
  end interface

  interface
     subroutine shses (nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                       wshses, lshses, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, isym, nt, idg, jdg, mdab, ndab, &
                               lshses, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: g(idg,jdg,*), a(mdab,ndab,*), b(mdab,ndab,*), &
                               wshses(lshses), work(lwork)
     end subroutine shses
  end interface

  interface
     subroutine shsec (nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                       wshsec, lshsec, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, isym, nt, idg, jdg, mdab, ndab, &
                               lshsec, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: g(idg,jdg,*), a(mdab,ndab,*), b(mdab,ndab,*), &
                               wshsec(lshsec), work(lwork)
     end subroutine shsec
  end interface

  interface
     subroutine shsgc (nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                       wshsgc, lshsgc, work, lwork, ierror)
       implicit none
       integer, intent(in)  :: nlat, nlon, isym, nt, idg, jdg, mdab, ndab, &
                               lshsgc, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: g(idg,jdg,nt),                    &
                               a(mdab,ndab,nt), b(mdab,ndab,nt), &
                               wshsgc(lshsgc), work(lwork)
     end subroutine shsgc
  end interface

  interface
    subroutine shsesi (nlat, nlon, wshses, lshses, &
                       work, lwork, dwork, ldwork, ierror)
       integer, intent(in)  :: nlat, nlon, lshses, lwork, ldwork
       integer, intent(out) :: ierror
       real(_wp_)           :: wshses(*),work(*)
       double precision     :: dwork(*)
     end subroutine shsesi
  end interface

  interface
    subroutine shseci (nlat, nlon, wshsec, lshsec, dwork, ldwork, ierror)
       integer, intent(in)  :: nlat, nlon, lshsec, ldwork
       integer, intent(out) :: ierror
       real(_wp_)           :: wshsec(*)
       double precision     :: dwork(*)
     end subroutine shseci
  end interface

  interface
    subroutine shsgci (nlat, nlon, wshsgc, lshsgc, dwork, ldwork, ierror)
       implicit none
       integer, intent(in)  :: nlat, nlon, lshsgc, ldwork
       integer, intent(out) :: ierror
       real(_wp_)           :: wshsgc(lshsgc)
       double precision     :: dwork(ldwork)
     end subroutine shsgci
  end interface

  interface
     subroutine vhaes (nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                       br, bi, cr, ci, mdab, ndab,             &
                       wvhaes, lvhaes, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, ityp, nt, idvw, jdvw, &
                               mdab, ndab, lvhaes, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: v(idvw,jdvw,*), w(idvw,jdvw,*), &
                               br(mdab,ndab,*), bi(mdab,ndab,*), &
                               cr(mdab,ndab,*), ci(mdab,ndab,*), &
                               work(lwork), wvhaes(lvhaes)
     end subroutine vhaes
  end interface

  interface
     subroutine vhaec (nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                       br, bi, cr, ci, mdab, ndab,             &
                       wvhaec, lvhaec, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, ityp, nt, idvw, jdvw, &
                               mdab, ndab, lvhaec, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: v(idvw,jdvw,*), w(idvw,jdvw,*), &
                               br(mdab,ndab,*), bi(mdab,ndab,*), &
                               cr(mdab,ndab,*), ci(mdab,ndab,*), &
                               work(lwork), wvhaec(lvhaec)
     end subroutine vhaec
  end interface

  interface
     subroutine vhagc (nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                       br, bi, cr, ci, mdab, ndab,             &
                       wvhagc, lvhagc, work, lwork, ierror)
       implicit none
       integer, intent(in)  :: nlat, nlon, ityp, nt, idvw, jdvw, &
                               mdab, ndab, lvhagc, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: v(idvw,jdvw,nt), w(idvw,jdvw,nt), &
                               br(mdab,ndab,nt), bi(mdab,ndab,nt), &
                               cr(mdab,ndab,nt), ci(mdab,ndab,nt), &
                               work(lwork), wvhagc(lvhagc)
     end subroutine vhagc
  end interface

  interface
     subroutine vhaesi (nlat, nlon, wvhaes, lvhaes, work, lwork, dwork, &
          ldwork, ierror)
       integer, intent(in)  :: nlat, nlon, lvhaes, lwork, ldwork
       real(_wp_)           :: wvhaes(lvhaes), work(lwork)
       integer, intent(out) :: ierror
       double precision     :: dwork(ldwork)
     end subroutine vhaesi
  end interface

  interface
     subroutine vhaeci (nlat, nlon, wvhaec, lvhaec, dwork, ldwork, ierror)
       integer, intent(in)  :: nlat, nlon, lvhaec, ldwork
       real(_wp_)           :: wvhaec(lvhaec)
       integer, intent(out) :: ierror
       double precision     :: dwork(ldwork)
     end subroutine vhaeci
  end interface

  interface
     subroutine vhagci (nlat, nlon, wvhagc, lvhagc, dwork, ldwork, ierror)
       implicit none
       integer, intent(in)  :: nlat, nlon, lvhagc, ldwork
       real(_wp_)           :: wvhagc(lvhagc)
       integer, intent(out) :: ierror
       double precision     :: dwork(ldwork)
     end subroutine vhagci
  end interface

  interface
     subroutine vhsgsi (nlat, nlon, wvhsgs, lvhsgs, dwork, ldwork, ierror)
       integer, intent(in)  :: nlat, nlon, lvhsgs, ldwork
       integer, intent(out) :: ierror
       real(_wp_)           :: wvhsgs(lvhsgs)
       double precision     :: dwork(ldwork)
     end subroutine vhsgsi
  end interface

  interface
     subroutine vhsgc (nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                       br, bi, cr, ci, mdab, ndab,             &
                       wvhsgc, lvhsgc, work, lwork, ierror)
       implicit none
       integer, intent(in)  :: nlat, nlon, ityp, nt, idvw, jdvw, &
                               mdab, ndab, lvhsgc, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: v(idvw,jdvw,nt), w(idvw,jdvw,nt), &
                               br(mdab,ndab,nt), bi(mdab,ndab,nt), &
                               cr(mdab,ndab,nt), ci(mdab,ndab,nt), &
                               work(lwork), wvhsgc(lvhsgc)
     end subroutine vhsgc
  end interface

  interface
     subroutine vhsgci (nlat, nlon, wvhsgc, lvhsgc, dwork, ldwork, ierror)
       integer, intent(in)  :: nlat, nlon, lvhsgc, ldwork
       integer, intent(out) :: ierror
       real(_wp_)           :: wvhsgc(lvhsgc)
       double precision     :: dwork(ldwork)
     end subroutine vhsgci
  end interface

  interface
     subroutine gradgs (nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, &
                        mdab, ndab, wvhsgs, lvhsgs, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, isym, nt, idvw, jdvw, mdab, ndab, &
                               lvhsgs, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: v(idvw,jdvw,nt), w(idvw,jdvw,nt)
       real(_wp_)           :: a(mdab,ndab,nt), b(mdab,ndab,nt)
       real(_wp_)           :: wvhsgs(lvhsgs), work(lwork)
     end subroutine gradgs
  end interface

  interface
     subroutine dives (nlat, nlon, isym, nt, dv, idv, jdv, br, bi, &
                       mdb, ndb, wshses, lshses, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, isym, nt, idv, jdv, mdb, ndb, &
                               lshses, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: dv(idv,jdv,nt), br(mdb,ndb,nt), bi(mdb,ndb,nt)
       real(_wp_)           :: wshses(lshses), work(lwork)
     end subroutine dives
  end interface

  interface
     subroutine divec (nlat, nlon, isym, nt, dv, idv, jdv, br, bi, &
                       mdb, ndb, wshsec, lshsec, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, isym, nt, idv, jdv, mdb, ndb, &
                               lshsec, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: dv(idv,jdv,nt), br(mdb,ndb,nt), bi(mdb,ndb,nt)
       real(_wp_)           :: wshsec(lshsec), work(lwork)
     end subroutine divec
  end interface

  interface
     subroutine vrtes (nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, &
                       mdc, ndc, wshses, lshses, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, isym, nt, ivrt, jvrt, mdc, ndc, &
                               lshses, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: vort(ivrt,jvrt,nt), &
                               cr(mdc,ndc,nt), ci(mdc,ndc,nt)
       real(_wp_)           :: wshses(lshses), work(lwork)
     end subroutine vrtes
  end interface

  interface
     subroutine sfvpes (nlat, nlon, isym, nt, sf, vp, idv, jdv, &
                        br, bi, cr, ci, mdb, ndb,               &
                        wshses, lshses, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, isym, nt, idv, jdv, mdb, ndb, &
                               lshses, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: sf(idv,jdv,nt), vp(idv,jdv,nt)
       real(_wp_)           :: br(mdb,ndb,nt), bi(mdb,ndb,nt)
       real(_wp_)           :: cr(mdb,ndb,nt), ci(mdb,ndb,nt)
       real(_wp_)           :: wshses(lshses), work(lwork)
    end subroutine sfvpes
  end interface

  interface
     subroutine sfvpec (nlat, nlon, isym, nt, sf, vp, idv, jdv, &
                        br, bi, cr, ci, mdb, ndb,               &
                        wshsec, lshsec, work, lwork, ierror)
       integer, intent(in)  :: nlat, nlon, isym, nt, idv, jdv, mdb, ndb, &
                               lshsec, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: sf(idv,jdv,nt), vp(idv,jdv,nt)
       real(_wp_)           :: br(mdb,ndb,nt), bi(mdb,ndb,nt)
       real(_wp_)           :: cr(mdb,ndb,nt), ci(mdb,ndb,nt)
       real(_wp_)           :: wshsec(lshsec), work(lwork)
    end subroutine sfvpec
  end interface

  interface
     subroutine sfvpgc (nlat, nlon, isym, nt, sf, vp, idv, jdv, &
                        br, bi, cr, ci, mdb, ndb,               &
                        wshsgc, lshsgc, work, lwork, ierror)
       implicit none
       integer, intent(in)  :: nlat, nlon, isym, nt, idv, jdv, mdb, ndb, &
                               lshsgc, lwork
       integer, intent(out) :: ierror
       real(_wp_)           :: sf(idv,jdv,nt), vp(idv,jdv,nt)
       real(_wp_)           :: br(mdb,ndb,nt), bi(mdb,ndb,nt)
       real(_wp_)           :: cr(mdb,ndb,nt), ci(mdb,ndb,nt)
       real(_wp_)           :: wshsgc(lshsgc), work(lwork)
    end subroutine sfvpgc
  end interface

contains

  subroutine sp95_init ()
    ! Initialization
    ! Currently just a no-op...

!    real(wp), intent(in) :: rsphere
!    radius = rsphere
  end subroutine sp95_init

  ! --

  subroutine sp95_cleanup ()
    ! Deallocate auxiliary module variables

    if (allocated (wshaec)) deallocate (wshaec)
    if (allocated (wshags)) deallocate (wshags)
    if (allocated (wshagc)) deallocate (wshagc)
    if (allocated (wvhaes)) deallocate (wvhaes)
    if (allocated (wvhaec)) deallocate (wvhaec)
    if (allocated (wvhagc)) deallocate (wvhagc)

    if (allocated (wshses)) deallocate (wshses)
    if (allocated (wshsec)) deallocate (wshsec)
    if (allocated (wshsgc)) deallocate (wshsgc)
    if (allocated (wvhsgs)) deallocate (wvhsgs)
    if (allocated (wvhsgc)) deallocate (wvhsgc)

  end subroutine sp95_cleanup

  ! --

  subroutine sp95_vrtdiv_s (ugrid, vgrid, vrtgrid, divgrid, rsphere)
    !
    ! Given a wind field (u,v), compute vorticity and divergence.
    ! Uses intermediate storage to save calculations.
    !
    real(wp), dimension(:,:,:), intent(in)  :: ugrid, vgrid
    real(wp), dimension(:,:,:), intent(out) :: vrtgrid, divgrid
    real(wp),                   intent(in)  :: rsphere

#ifndef SPHEREPACK
call finish ('sp95_vrtdiv_s','libspherepack not linked')
#else

    ! Local variables
    integer :: nlat, nlon, isym, ityp, nt, &
               idv, jdv, idvw, jdvw, ivrt, jvrt, &
               mdab, ndab, mdb, ndb, mdc, ndc,   &
               lshses, lvhaes, lwork, ldwork, ierror
    integer :: l1, l2
    real(wp), allocatable, dimension(:,:,:) :: v, w, dv, vort
    real(wp), allocatable, dimension(:,:,:) :: br, bi, cr, ci

    real(wp),         allocatable :: work(:)
    double precision, allocatable :: dwork(:)

    ! Get shape of grid in geophysical conventions
    nlon = size (ugrid, dim=1)
    nlat = size (ugrid, dim=2)
    nt   = size (ugrid, dim=3)

    allocate (v(nlat,nlon,nt), w(nlat,nlon,nt))

    ! Initialize coefficients of vector spherical harmonic analysis
    l1 = min (nlat, (nlon+1)/2)
    l2 = (nlat+1)/2
    lvhaes = l1*l2*(2*nlat-l1+1)+nlon+15
    ! Size of lwork dominated by vhaes and/or dives
    lwork  = max ((2*nt+1)*nlat*nlon, &
                  nlat*((nt+1)*nlon+2*nt*(l1+1)+1))
    ldwork = 2*(nlat+1)

    allocate (work(lwork))
    allocate (dwork(ldwork))

    if (allocated (wvhaes)) then
       if (size (wvhaes) /= lvhaes) then
          print *, "size (wvhaes) /= lvhaes"
          stop
       end if
    else
       allocate (wvhaes(lvhaes))
    end if

    call vhaesi (nlat, nlon, wvhaes, lvhaes, &
                 work, lwork, dwork, ldwork, ierror)

    if (ierror /= 0) then
       print *, "vhaesi failed, ierror =", ierror
       stop
    end if

    ! Initialize coefficients of spherical harmonic synthesis
    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lshses = l1*l2*(2*nlat-l1+1)+nlon+15

    if (allocated (wshses)) then
       if (size (wshses) /= lshses) then
          print *, "size (wshses) /= lshses"
          stop
       end if
    else
       allocate (wshses(lshses))
    end if

    call shsesi (nlat, nlon, wshses, lshses, &
                 work, lwork, dwork, ldwork, ierror)

    if (ierror /= 0) then
       print *, "shsesi failed, ierror =", ierror
       stop
    end if

    ! Convert winds from geophysical to mathematical coordinates
    call sp95_geo2math_vector (ugrid, vgrid, v, w, 0)

    ! Vector spherical harmonic analysis
    ityp = 0
    idvw = nlat
    jdvw = nlon
    mdab = nlat
    ndab = nlon
    allocate (br(mdab,ndab,nt), bi(mdab,ndab,nt), &
              cr(mdab,ndab,nt), ci(mdab,ndab,nt))
    call vhaes (nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                br, bi, cr, ci, mdab, ndab,             &
                wvhaes, lvhaes, work, lwork, ierror)

    if (ierror /= 0) then
       print *, "vhaes failed, ierror =", ierror
       stop
    end if

    ! Calculate divergence ...
    isym = 0
    idv  = nlat
    jdv  = nlon
    mdb  = nlat
    ndb  = nlon
    allocate (dv(idv,jdv,nt))
    call dives (nlat, nlon, isym, nt, dv, idv, jdv, br, bi, &
                mdb, ndb, wshses, lshses, work, lwork, ierror)

    if (ierror /= 0) then
       print *, "dives failed, ierror =", ierror
       stop
    end if

    call sp95_math2geo_scaled (phi_math=dv, phi_geo=divgrid, scale=1/rsphere, ig=0)

    ! ... and vorticity
    isym = 0
    ivrt = nlat
    jvrt = nlon
    mdc  = nlat
    ndc  = nlon
    allocate (vort(ivrt,jvrt,nt))
    call vrtes (nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, &
                mdc, ndc, wshses, lshses, work, lwork, ierror)

    if (ierror /= 0) then
       print *, "vrtes failed, ierror =", ierror
       stop
    end if

    call sp95_math2geo_scaled (phi_math=vort, phi_geo=vrtgrid, scale=1/rsphere, ig=0)

    deallocate (br, bi, cr, ci)
    deallocate (dv, vort)
    deallocate (work, dwork)
    deallocate (v, w)

#endif

  end subroutine sp95_vrtdiv_s

  ! --

  subroutine sp95_vrtdiv (ugrid, vgrid, vrtgrid, divgrid, rsphere)
    !
    ! Given a wind field (u,v), compute vorticity and divergence.
    ! Saves intermediate storage by allowing coefficient recalculations.
    !
    real(wp), dimension(:,:,:), intent(in)  :: ugrid, vgrid
    real(wp), dimension(:,:,:), intent(out) :: vrtgrid, divgrid
    real(wp),                   intent(in)  :: rsphere

#ifndef SPHEREPACK
call finish ('sp95_vrtdiv','libspherepack not linked')
#else

    ! Local variables
    integer :: nlat, nlon, isym, ityp, nt, &
               idv, jdv, idvw, jdvw, ivrt, jvrt, &
               mdab, ndab, mdb, ndb, mdc, ndc,   &
               lshsec, lvhaec, lwork, ldwork, ierror
    integer :: l1, l2
    real(wp), allocatable, dimension(:,:,:) :: v, w, dv, vort
    real(wp), allocatable, dimension(:,:,:) :: br, bi, cr, ci

    real(wp),         allocatable :: work(:)
    double precision, allocatable :: dwork(:)

    ! Get shape of grid in geophysical conventions
    nlon = size (ugrid, dim=1)
    nlat = size (ugrid, dim=2)
    nt   = size (ugrid, dim=3)

    allocate (v(nlat,nlon,nt), w(nlat,nlon,nt))

    ! Initialize coefficients of vector spherical harmonic analysis
    l1 = min (nlat, (nlon+1)/2)
    l2 = (nlat+1)/2
    lvhaec = 4*nlat*l2+3*max (l1-2,0)*(2*nlat-l1-1)+nlon+15
    ! Size of lwork dominated by vhaec and/or divec
    lwork  = max (nlat*(2*nt*nlon+max (6*l2,nlon)), &
                  nlat*(nt*nlon+max (3*l2,nlon)+2*nt*(l1+1)+1))
    ldwork = 2*(nlat+2)

    allocate (work(lwork))
    allocate (dwork(ldwork))

    if (allocated (wvhaec)) then
       if (size (wvhaec) /= lvhaec) then
          print *, "size (wvhaec) /= lvhaec"
          stop
       end if
    else
       allocate (wvhaec(lvhaec))
    end if

    call vhaeci (nlat, nlon, wvhaec, lvhaec, dwork, ldwork, ierror)

    if (ierror /= 0) then
       print *, "vhaeci failed, ierror =", ierror
       stop
    end if

    ! Initialize coefficients of spherical harmonic synthesis
    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lshsec = 2*nlat*l2+3*((l1-2)*(2*nlat-l1-1))/2+nlon+15

    if (allocated (wshsec)) then
       if (size (wshsec) /= lshsec) then
          print *, "size (wshsec) /= lshsec"
          stop
       end if
    else
       allocate (wshsec(lshsec))
    end if

    call shseci (nlat, nlon, wshsec, lshsec, dwork, ldwork, ierror)

    if (ierror /= 0) then
       print *, "shseci failed, ierror =", ierror
       stop
    end if

    ! Convert winds from geophysical to mathematical coordinates
    call sp95_geo2math_vector (ugrid, vgrid, v, w, 0)

    ! Vector spherical harmonic analysis
    ityp = 0
    idvw = nlat
    jdvw = nlon
    mdab = nlat
    ndab = nlon
    allocate (br(mdab,ndab,nt), bi(mdab,ndab,nt), &
              cr(mdab,ndab,nt), ci(mdab,ndab,nt))
    call vhaec (nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                br, bi, cr, ci, mdab, ndab,             &
                wvhaec, lvhaec, work, lwork, ierror)

    if (ierror /= 0) then
       print *, "vhaec failed, ierror =", ierror
       stop
    end if

    ! Calculate divergence ...
    isym = 0
    idv  = nlat
    jdv  = nlon
    mdb  = nlat
    ndb  = nlon
    allocate (dv(idv,jdv,nt))
    call divec (nlat, nlon, isym, nt, dv, idv, jdv, br, bi, &
                mdb, ndb, wshsec, lshsec, work, lwork, ierror)

    if (ierror /= 0) then
       print *, "divec failed, ierror =", ierror
       stop
    end if

    call sp95_math2geo_scaled (phi_math=dv, phi_geo=divgrid, scale=1/rsphere, ig=0)

    ! ... and vorticity
    isym = 0
    ivrt = nlat
    jvrt = nlon
    mdc  = nlat
    ndc  = nlon
    allocate (vort(ivrt,jvrt,nt))
    call vrtec (nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, &
                mdc, ndc, wshsec, lshsec, work, lwork, ierror)

    if (ierror /= 0) then
       print *, "vrtec failed, ierror =", ierror
       stop
    end if

    call sp95_math2geo_scaled (phi_math=vort, phi_geo=vrtgrid, scale=1/rsphere, ig=0)

    deallocate (br, bi, cr, ci)
    deallocate (dv, vort)
    deallocate (work, dwork)
    deallocate (v, w)

#endif

  end subroutine sp95_vrtdiv

  ! --

  subroutine sp95_sfvp_s (ugrid, vgrid, psigrid, chigrid, rsphere)
    !
    ! Given a wind field (u,v), compute streamfunction psi and
    ! velocity potential chi.
    ! Uses intermediate storage to save calculations.
    !
    real(wp), dimension(:,:,:), intent(in)  :: ugrid, vgrid
    real(wp), dimension(:,:,:), intent(out) :: psigrid, chigrid
    real(wp),                   intent(in)  :: rsphere

#ifndef SPHEREPACK
call finish ('sfvp_s','libspherepack not linked')
#else


    ! Local variables
    integer :: nlat, nlon, isym, ityp, nt, &
               idv, jdv, idvw, jdvw, mdab, ndab, mdb, ndb, &
               lshses, lvhaes, lwork, ldwork, ierror
    integer :: l1, l2
    real(wp), allocatable, dimension(:,:,:) :: v, w, sf, vp
    real(wp), allocatable, dimension(:,:,:) :: br, bi, cr, ci

    real(wp),         allocatable :: work(:)
    double precision, allocatable :: dwork(:)

    ! Get shape of grid in geophysical conventions
    nlon = size (ugrid, dim=1)
    nlat = size (ugrid, dim=2)
    nt   = size (ugrid, dim=3)

    allocate (v(nlat,nlon,nt), w(nlat,nlon,nt))

    ! Initialize coefficients of vector spherical harmonic analysis
    l1 = min (nlat, (nlon+1)/2)
    l2 = (nlat+1)/2
    lvhaes = l1*l2*(2*nlat-l1+1)+nlon+15
    ! Size of lwork dominated by either vhaes or sfvpes
    lwork  = max ((2*nt+1)*nlat*(nlon+1), &
                  nlat*((nt+1)*nlon + 2*l2*nt+1))
    ldwork = 2*(nlat+1)

    allocate (work(lwork))
    allocate (dwork(ldwork))

    if (allocated (wvhaes)) then
       if (size (wvhaes) /= lvhaes) then
          print *, "size (wvhaes) /= lvhaes"
          stop
       end if
    else
       allocate (wvhaes(lvhaes))
    end if

    call vhaesi (nlat, nlon, wvhaes, lvhaes, &
                 work, lwork, dwork, ldwork, ierror)

    if (ierror /= 0) then
       print *, "vhaesi failed, ierror =", ierror
       stop
    end if

    ! Initialize coefficients of spherical harmonic synthesis
    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lshses = l1*l2*(2*nlat-l1+1)+nlon+15

    if (allocated (wshses)) then
       if (size (wshses) /= lshses) then
          print *, "size (wshses) /= lshses"
          stop
       end if
    else
       allocate (wshses(lshses))
    end if

    call shsesi (nlat, nlon, wshses, lshses, &
                 work, lwork, dwork, ldwork, ierror)

    if (ierror /= 0) then
       print *, "shsesi failed, ierror =", ierror
       stop
    end if

    ! Convert winds from geophysical to mathematical coordinates
    call sp95_geo2math_vector (ugrid, vgrid, v, w, 0)

    ! Vector spherical harmonic analysis
    ityp = 0
    idvw = nlat
    jdvw = nlon
    mdab = nlat
    ndab = nlon
    allocate (br(mdab,ndab,nt), bi(mdab,ndab,nt), &
              cr(mdab,ndab,nt), ci(mdab,ndab,nt))
    call vhaes (nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                br, bi, cr, ci, mdab, ndab,             &
                wvhaes, lvhaes, work, lwork, ierror)

    if (ierror /= 0) then
       print *, "vhaes failed, ierror =", ierror
       stop
    end if

    ! Calculate streamfunction and velocity potential
    isym = 0
    idv  = nlat
    jdv  = nlon
    mdb  = nlat
    ndb  = nlon
    allocate (sf(idv,jdv,nt), vp(idv,jdv,nt))
    call sfvpes (nlat, nlon, isym, nt, sf, vp, idv, jdv, &
                 br, bi, cr, ci, mdb, ndb,               &
                 wshses, lshses, work, lwork, ierror)

    if (ierror /= 0) then
       print *, "sfvpes failed, ierror =", ierror
       stop
    end if

    call sp95_math2geo_scaled (phi_math=sf, phi_geo=psigrid, scale=rsphere, ig=0)
    call sp95_math2geo_scaled (phi_math=vp, phi_geo=chigrid, scale=rsphere, ig=0)

    deallocate (sf, vp)
    deallocate (br, bi, cr, ci)
    deallocate (work, dwork)
    deallocate (v, w)

#endif

  end subroutine sp95_sfvp_s

  ! --

  subroutine sp95_sfvp (ugrid, vgrid, psigrid, chigrid, rsphere)
    !
    ! Given a wind field (u,v), compute streamfunction psi and
    ! velocity potential chi.
    ! Saves intermediate storage by allowing coefficient recalculations.
    !
    real(wp), dimension(:,:,:), intent(in)  :: ugrid, vgrid
    real(wp), dimension(:,:,:), intent(out) :: psigrid, chigrid
    real(wp),                   intent(in)  :: rsphere

#ifndef SPHEREPACK
call finish ('sp95_sfvp','libspherepack not linked')
#else

    ! Local variables
    integer :: nlat, nlon, isym, ityp, nt, &
               idv, jdv, idvw, jdvw, mdab, ndab, mdb, ndb, &
               lshsec, lvhaec, lwork, ldwork, ierror
    integer :: l1, l2
    real(wp), allocatable, dimension(:,:,:) :: v, w, sf, vp
    real(wp), allocatable, dimension(:,:,:) :: br, bi, cr, ci

    real(wp),         allocatable :: work(:)
    double precision, allocatable :: dwork(:)

    ! Get shape of grid in geophysical conventions
    nlon = size (ugrid, dim=1)
    nlat = size (ugrid, dim=2)
    ! Number of analyses or levels
    nt   = size (ugrid, dim=3)

    allocate (v(nlat,nlon,nt), w(nlat,nlon,nt))

    ! Initialize coefficients of vector spherical harmonic analysis
    l1 = min (nlat, (nlon+1)/2)
    l2 = (nlat+1)/2
    lvhaec = 4*nlat*l2+3*max (l1-2,0)*(2*nlat-l1-1)+nlon+15
    ! Size of lwork determined from vhaec and sfvpec
    lwork  = max (nlat*(2*nt*nlon+max (6*l2,nlon)), &
                  l2*(nt*nlon+max (3*nlat,nlon)) + nlat*(2*l1*nt+1))
    ldwork = 2*(nlat+2)

    allocate (work(lwork))
    allocate (dwork(ldwork))

    if (allocated (wvhaec)) then
       if (size (wvhaec) /= lvhaec) then
          print *, "size (wvhaec) /= lvhaec"
          stop
       end if
    else
       allocate (wvhaec(lvhaec))
    end if

    call vhaeci (nlat, nlon, wvhaec, lvhaec, dwork, ldwork, ierror)

    if (ierror /= 0) then
       print *, "vhaeci failed, ierror =", ierror
       stop
    end if

    ! Initialize coefficients of spherical harmonic synthesis
    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lshsec = 2*nlat*l2+3*((l1-2)*(2*nlat-l1-1))/2+nlon+15

    if (allocated (wshsec)) then
       if (size (wshsec) /= lshsec) then
          print *, "size (wshsec) /= lshsec"
          stop
       end if
    else
       allocate (wshsec(lshsec))
    end if

    call shseci (nlat, nlon, wshsec, lshsec, dwork, ldwork, ierror)

    if (ierror /= 0) then
       print *, "shseci failed, ierror =", ierror
       stop
    end if

    ! Convert winds from geophysical to mathematical coordinates
    call sp95_geo2math_vector (ugrid, vgrid, v, w, 0)

    ! Vector spherical harmonic analysis
    ityp = 0
    idvw = nlat
    jdvw = nlon
    mdab = nlat
    ndab = nlon
    allocate (br(mdab,ndab,nt), bi(mdab,ndab,nt), &
              cr(mdab,ndab,nt), ci(mdab,ndab,nt))
    call vhaec (nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                br, bi, cr, ci, mdab, ndab,             &
                wvhaec, lvhaec, work, lwork, ierror)

    if (ierror /= 0) then
       print *, "vhaec failed, ierror =", ierror
       stop
    end if

    ! Calculate streamfunction and velocity potential
    isym = 0
    idv  = nlat
    jdv  = nlon
    mdb  = nlat
    ndb  = nlon
    allocate (sf(idv,jdv,nt), vp(idv,jdv,nt))
    call sfvpec (nlat, nlon, isym, nt, sf, vp, idv, jdv, &
                 br, bi, cr, ci, mdb, ndb,               &
                 wshsec, lshsec, work, lwork, ierror)

    if (ierror /= 0) then
       print *, "sfvpec failed, ierror =", ierror
       stop
    end if

    call sp95_math2geo_scaled (phi_math=sf, phi_geo=psigrid, scale=rsphere, ig=0)
    call sp95_math2geo_scaled (phi_math=vp, phi_geo=chigrid, scale=rsphere, ig=0)

    deallocate (sf, vp)
    deallocate (br, bi, cr, ci)
    deallocate (work, dwork)
    deallocate (v, w)

#endif

  end subroutine sp95_sfvp
  !============================================================================
  subroutine sp95_sfvp_gg (ugrid, vgrid, psigrid, chigrid, rsphere, ig)
    !
    ! Given a wind field (u,v), compute streamfunction psi and
    ! velocity potential chi for a gaussian grid.
    ! Saves intermediate storage by allowing coefficient recalculations.
    !
    real(wp), intent(in)  :: ugrid   (:,:,:)
    real(wp), intent(in)  :: vgrid   (:,:,:)
    real(wp), intent(out) :: psigrid (:,:,:)
    real(wp), intent(out) :: chigrid (:,:,:)
    real(wp), intent(in)  :: rsphere         ! radius of the sphere
    integer,  intent(in)  :: ig              ! ordering 0: S->N;  1-> N->S

#ifndef SPHEREPACK
call finish ('sp95_sfvp_gg','libspherepack not linked')
#else

    ! Local variables
    integer :: nlat, nlon, isym, ityp, nt, idvw, jdvw, mdab, ndab, &
               lshsgc, lvhagc, lwork, ierror
    integer :: l1, l2
    real(wp), allocatable :: work(:)
    real(wp), allocatable, dimension(:,:,:) :: v, w, sf, vp
    real(wp), allocatable, dimension(:,:,:) :: br, bi, cr, ci

    ! Get shape of grid in geophysical conventions
    nlon = size (ugrid, dim=1)
    nlat = size (ugrid, dim=2)
    ! Number of analyses or levels
    nt   = size (ugrid, dim=3)

    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    ! Size of lwork determined from vhagc and sfvpgc
    lwork  = max (2*nlat*(2*nlon*nt+3*l2)                      &! vhagc
                 ,nlat*((nt*nlon+max0(3*l2,nlon))+2*l1*nt+1))   ! sfvpgc
    allocate (work(lwork))

    ! Initialize coefficients of vector spherical harmonic analysis
    call set_wvhagc (nlat, nlon)
    lvhagc = size (wvhagc)

    ! Initialize coefficients of scalar spherical harmonic synthesis
    call set_wshsgc (nlat, nlon)
    lshsgc = size (wshsgc)

    ! Convert winds from geophysical to mathematical coordinates
    idvw = nlat
    jdvw = nlon
    allocate (v(nlat,nlon,nt), w(nlat,nlon,nt))
    call sp95_geo2math_vector (ugrid, vgrid, v, w, ig)

    ! Vector spherical harmonic analysis
    ityp = 0
    mdab = l1
    ndab = nlon
    allocate (br(mdab,ndab,nt), bi(mdab,ndab,nt), &
              cr(mdab,ndab,nt), ci(mdab,ndab,nt))
    call vhagc (nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                br, bi, cr, ci, mdab, ndab,             &
                wvhagc, lvhagc, work, lwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_sfvp_gg","vhagc")
    end if

    ! Calculate streamfunction and velocity potential
    isym = 0
    allocate (sf(idvw,jdvw,nt), vp(idvw,jdvw,nt))
    call sfvpgc (nlat, nlon, isym, nt, sf, vp, idvw, jdvw, &
                 br, bi, cr, ci, mdab, ndab,               &
                 wshsgc, lshsgc, work, lwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_sfvp_gg","sfvpgc")
    end if

    call sp95_math2geo_scaled (phi_math=sf, phi_geo=psigrid, scale=rsphere, ig=ig)
    call sp95_math2geo_scaled (phi_math=vp, phi_geo=chigrid, scale=rsphere, ig=ig)

    deallocate (sf, vp)
    deallocate (br, bi, cr, ci)
    deallocate (work)
    deallocate (v, w)

#endif

  end subroutine sp95_sfvp_gg
  !============================================================================
  subroutine sp95_sfvp_vrtdiv_gg (vrt, div, psi, chi, &
                                  ugrid, vgrid, rsphere, ig)
    !
    ! Given the vorticity and divergence of a vector field,
    ! compute streamfunction psi, velocity potential chi,
    ! and wind field (u,v) for a Gaussian grid.
    ! Saves intermediate storage by allowing coefficient recalculations.
    !
    real(wp), intent(in)  :: vrt     (:,:,:)    ! Vorticity
    real(wp), intent(in)  :: div     (:,:,:)    ! Divergence
    real(wp), intent(out) :: psi     (:,:,:)    ! Streamfunction
    real(wp), intent(out) :: chi     (:,:,:)    ! Velocity potential
    real(wp), intent(out) :: ugrid   (:,:,:)    ! u-wind
    real(wp), intent(out) :: vgrid   (:,:,:)    ! v-wind
    real(wp), intent(in)  :: rsphere            ! Radius of the sphere
    integer,  intent(in)  :: ig                 ! ordering 0: S->N;  1-> N->S

#ifndef SPHEREPACK
call finish ('sp95_sfvp_vrtdiv_gg','libspherepack not linked')
#else

    ! Local variables
    integer  :: nlat, nlon, isym, ityp, nt, idvw, jdvw, mdab, ndab, &
                lshagc, lshsgc, lvhsgc, lwork, ierror, n, mm
    integer  :: l1, l2, mmax
    real(wp) :: fac
    real(wp), allocatable :: work(:)
    real(wp), allocatable, dimension(:,:,:) :: v, w, sf, vp
    real(wp), allocatable, dimension(:,:,:) :: br, bi, cr, ci

    ! Get shape of grid in geophysical conventions
    nlon = size (vrt, dim=1)
    nlat = size (vrt, dim=2)
    ! Number of analyses or levels
    nt   = size (vrt, dim=3)

    if ( any (shape (vrt) /= shape (div))     .or. &
         any (shape (vrt) /= shape (psi)) .or. &
         any (shape (vrt) /= shape (chi)) .or. &
         any (shape (vrt) /= shape (ugrid))   .or. &
         any (shape (vrt) /= shape (vgrid)) ) then
       write (0,*) shape (vrt), shape (div)
       write (0,*) shape (psi), shape (chi)
       write (0,*) shape (ugrid), shape (vgrid)
       call finish ('sp95_sfvp_vrtdiv_gg','non-conforming shapes!')
    end if
    if (ig /= 0) then
       call finish ('sp95_sfvp_vrtdiv_gg','unsupported grid orientation N->S!')
    end if

    ! Initialize coefficients of scalar spherical harmonic analysis
    call set_wshagc (nlat, nlon)
    lshagc = size (wshagc)

    ! Initialize coefficients of scalar spherical harmonic synthesis
    call set_wshsgc (nlat, nlon)
    lshsgc = size (wshsgc)

    ! Initialize coefficients of vector spherical harmonic synthesis
    call set_wvhsgc (nlat, nlon)
    lvhsgc = size (wvhsgc)

    ! Size of lwork determined from shagc, shsgc and vhsgc
    l1 = min (nlat, (nlon+2)/2)                         ! sh*
    l2 = (nlat+1)/2
    lwork = max (nlat*(  nt*nlon+max (3*l2,nlon))      &! shagc, shsgc
                ,nlat*(2*nt*nlon+max (6*l2,nlon)))      ! vhsgc
    allocate (work(lwork))

    mmax = min (nlat, (nlon+1)/2)
    idvw = nlat
    jdvw = nlon
    allocate (v(nlat,nlon,nt), w(nlat,nlon,nt))
    allocate (sf(idvw,jdvw,nt), vp(idvw,jdvw,nt))
    mdab = l1                                           ! sh*
    ndab = nlat
    allocate (br(mdab,ndab,nt), bi(mdab,ndab,nt), &
              cr(mdab,ndab,nt), ci(mdab,ndab,nt))

    ! Convert vorticity, divergence from geophysical to mathematical
    ! coordinates (Spherepack's internal representation)
    call sp95_geo2math_scalar (vrt, v)
    call sp95_geo2math_scalar (div, w)

    ! Scalar harmonic analysis
    isym = 0
    call shagc (nlat, nlon, isym, nt, w, idvw, jdvw, br, bi, mdab, ndab, &
                wshagc, lshagc, work, lwork, ierror)
    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_sfvp_vrtdiv_gg", "shagc failed to convert div")
    end if

    call shagc (nlat, nlon, isym, nt, v, idvw, jdvw, cr, ci, mdab, ndab, &
                wshagc, lshagc, work, lwork, ierror)
    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_sfvp_vrtdiv_gg", "shagc failed to convert vrt")
    end if

    ! Check/set mode 0 to zero and apply inverse Laplacian
    br(:,1,:) = 0
    bi(:,1,:) = 0
    cr(:,1,:) = 0
    ci(:,1,:) = 0
    do n = 2, nlat
       fac = - Rsphere**2 / (n*(n-1))
       mm = min (n, mmax)
       br(:mm,n,:) = br(:mm,n,:) * fac
       bi(:mm,n,:) = bi(:mm,n,:) * fac
       cr(:mm,n,:) = cr(:mm,n,:) * fac
       ci(:mm,n,:) = ci(:mm,n,:) * fac
    end do
    br(mmax+1:,:,:) = 0
    bi(mmax+1:,:,:) = 0
    cr(mmax+1:,:,:) = 0
    bi(mmax+1:,:,:) = 0

    ! Synthesize streamfunction and velocity potential
    call shsgc (nlat, nlon, isym, nt, vp, idvw, jdvw, br, bi, mdab, ndab, &
                wshsgc, lshsgc, work, lwork, ierror)
    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_sfvp_vrtdiv_gg", "shsgc failed to convert vp")
    end if

    call shsgc (nlat, nlon, isym, nt, sf, idvw, jdvw, cr, ci, mdab, ndab, &
                wshsgc, lshsgc, work, lwork, ierror)
    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_sfvp_vrtdiv_gg", "shsgc failed to convert sf")
    end if

    ! Copy fields from Spherepack's internal representation
    call sp95_math2geo_scalar (phi_math=sf, phi_geo=psi)
    call sp95_math2geo_scalar (phi_math=vp, phi_geo=chi)

    ! Calculate winds (u,v) from sf, vp
    do n = 2, nlat
       fac = sqrt (real (n*(n-1), wp)) / Rsphere
       mm = min (n, mmax)
       br(:mm,n,:) = br(:mm,n,:) *   fac
       bi(:mm,n,:) = bi(:mm,n,:) *   fac
       cr(:mm,n,:) = cr(:mm,n,:) * (-fac)
       ci(:mm,n,:) = ci(:mm,n,:) * (-fac)
    end do

    ! Vector spherical harmonic synthesis
    ityp = 0
    call vhsgc (nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                mdab, ndab, wvhsgc, lvhsgc, work, lwork, ierror)
    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_sfvp_vrtdiv_gg", "vhsgc failed to synthesize winds")
    end if

    ! Convert winds from Spherepack's internal representation
    call sp95_math2geo_vector (v_math=v, w_math=w, u_geo=ugrid, v_geo=vgrid)

    deallocate (br, bi, cr, ci)
    deallocate (sf, vp)
    deallocate (v, w)
    deallocate (work)

#endif

  end subroutine sp95_sfvp_vrtdiv_gg
  !============================================================================
  subroutine sp95_grad_gg (f, f_x, f_y)
    !--------------------------------------------------------------
    ! Calculate gradient of a scalar function f on a Gaussian grid.
    ! f_x = eastern (zonal), f_y = northern (meridional) component.
    !--------------------------------------------------------------
    real(wp), intent(in)  :: f(:,:)
    real(wp), intent(out) :: f_x(:,:), f_y(:,:)

#ifndef SPHEREPACK
call finish ('sp95_grad_gg','libspherepack not linked')
#else
    !----------------
    ! Local variables
    !----------------
    integer :: nx, ny, nt, isym, ierror, l1, l2, lshags, lvhsgs
    integer :: lwork, ldwork, idg, jdg, idvw, jdvw, mdab, ndab
    logical :: need_init
    real(wp), allocatable, dimension(:)     :: work
    real(wp), allocatable, dimension(:,:,:) :: g, a, b, v, w
    double precision, allocatable           :: dwork(:)

    nx = size (f, dim=1)
    ny = size (f, dim=2)
    nt = 1
    isym = 0

    need_init = .false.
    if (nx /= nlon .or. ny /= nlat) then
       need_init = .true.
       nlon = nx
       nlat = ny
    end if

    if (allocated (wshags)) then
       if (need_init) deallocate (wshags)
    end if
    if (allocated (wvhsgs)) then
       if (need_init) deallocate (wvhsgs)
    end if
    !--------------------------------------------
    ! Spherical harmonic analysis of scalar field
    !--------------------------------------------
    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lshags = nlat*(3*(l1+l2)-2) + (l1-1)*(l2*(2*nlat-l1)-3*l1)/2 + nlon + 15
    ! Size of lwork for shags and/or gradgs
    lwork = max (4*nlat*(nlat+2) + 2, &
                 nlat*((2*nt+1)*nlon + 2*l1*nt + 1))
    allocate (work(lwork))
    !------------------------------
    ! Set up coefficients if needed
    !------------------------------
    if (.not. allocated (wshags)) then
       allocate (wshags(lshags))

       ldwork = nlat * (nlat+4)
       allocate (dwork(ldwork))

       call shagsi (nlat, nlon, wshags, lshags, &
                    work, lwork, dwork, ldwork, ierror)

       if (ierror /= 0) then
          write (0,*) "ierror =", ierror
          call finish ("sp95_grad_gg", "failure in shagsi")
       end if
       deallocate (dwork)
    end if

    ! Copy field to Spherepack's internal representation (mathematical coord.)
    idg = nlat
    jdg = nlon
    allocate (g(idg,jdg,1))
    call sp95_geo2math_scalar_2d (phi_geo=f, phi_math=g(:,:,1))

    mdab = l1
    ndab = nlat
    allocate (a(mdab,ndab,1))
    allocate (b(mdab,ndab,1))

    ! Spherical harmonic analysis of a scalar field
    call shags (nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                wshags, lshags, work, lwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_grad_gg", "failure in shags")
    end if

    deallocate (g)
    !------------------------
    ! Calculation of gradient
    !------------------------
    l1 = min (nlat, (nlon+1)/2)
    l2 = (nlat+1)/2
    lvhsgs = l1*l2*(2*nlat-l1+1) + nlon + 15 + 2*nlat
    !-------------------------------------------------------------
    ! Set up coefficients for vector spherical harmonics synthesis
    !-------------------------------------------------------------
    if (.not. allocated (wvhsgs)) then
       allocate (wvhsgs(lvhsgs))

       ldwork = (3*nlat*(nlat+3)+2) / 2
       allocate (dwork(ldwork))

       call vhsgsi (nlat, nlon, wvhsgs, lvhsgs, dwork, ldwork, ierror)

       if (ierror /= 0) then
          write (0,*) "ierror =", ierror
          call finish ("sp95_grad_gg", "failure in vhsgsi")
       end if
       deallocate (dwork)
    end if

    idvw = nlat
    jdvw = nlon
    allocate (v(idvw,jdvw,1))
    allocate (w(idvw,jdvw,1))

    ! Calculate gradient
    call gradgs (nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, &
                 mdab, ndab, wvhsgs, lvhsgs, work, lwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_grad_gg", "failure in gradgs")
    end if

    ! Copy vector fields from Spherepack's internal representation
    call sp95_math2geo_vector_2d (v_math=v(:,:,1), w_math=w(:,:,1), &
                                  u_geo =f_x(:,:), v_geo =f_y(:,:))

    ! Final cleanups
    deallocate (a, b, v, w)
    deallocate (work)
#endif
  end subroutine sp95_grad_gg
  !============================================================================
  subroutine sp95_sha_gg (f, a, b)
    !--------------------------------------------------------------
    ! Calculate spherical harmonic analysis of a scalar function f
    ! on a Gaussian grid.
    ! Output: a,b = real & imaginary parts of spectral coefficients
    !         Coefficients(m,n) of a, b: m=1...mmax, n=m..nlat,
    !         where m is the maximum(+1) longitudinal wavenumber.
    ! (For the normalization used c.f. sections 4 and 5 of the
    ! SPHEREPACK manual.)
    !--------------------------------------------------------------
    real(wp), intent(in)  :: f(:,:)
    real(wp), intent(out) :: a(:,:), b(:,:)

#ifndef SPHEREPACK
call finish ('sp95_sha_gg','libspherepack not linked')
#else
    !----------------
    ! Local variables
    !----------------
    integer :: nx, ny, nt, isym, ierror, l1, l2, lshagc
    integer :: lwork, idg, jdg, mdab, ndab, mmax
    logical :: need_init
    real(wp), allocatable, dimension(:)     :: work
    real(wp), allocatable, dimension(:,:,:) :: g, aa, bb

    nx = size (f, dim=1)
    ny = size (f, dim=2)
    nt = 1
    isym = 0
    if (size (a,dim=2) /= ny .or. size (b,dim=2) /= ny .or. &
         any (shape (a) /= shape (b))) then
       write (0,*) "shape (f) =", shape (f)
       write (0,*) "shape (a) =", shape (a)
       write (0,*) "shape (b) =", shape (b)
       call finish ("sp95_sha_gg", "incompatible shapes")
    end if

    need_init = .false.
    if (nx /= nlon .or. ny /= nlat) then
       need_init = .true.
       nlon = nx
       nlat = ny
    end if

    if (allocated (wshagc)) then
       if (need_init) deallocate (wshagc)
    end if
    call set_wshagc (nlat, nlon)
    lshagc = size (wshagc)

    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lwork = nlat*(nt*nlon + max (3*l2,nlon))    ! Size of lwork for shagc
    allocate (work(lwork))

    !------------------------------------
    ! Maximum(+1) longitudinal wavenumber
    !------------------------------------
    mmax = min (l1, size (a,dim=1), size(b,dim=1))

    ! Copy field to Spherepack's internal representation (mathematical coord.)
    idg = nlat
    jdg = nlon
    allocate (g(idg,jdg,1))
    call sp95_geo2math_scalar_2d (phi_geo=f, phi_math=g(:,:,1))

    mdab = l1
    ndab = nlat
    allocate (aa(mdab,ndab,1))
    allocate (bb(mdab,ndab,1))

    ! Spherical harmonic analysis of a scalar field
    call shagc (nlat, nlon, isym, nt, g, idg, jdg, aa, bb, mdab, ndab, &
                wshagc, lshagc, work, lwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_sha_gg", "failure in shagc")
    end if

    a(1:mmax, :) = aa(1:mmax,:,1)
    a(mmax+1:,:) = 0
    b(1:mmax, :) = bb(1:mmax,:,1)
    b(mmax+1:,:) = 0

    ! Final cleanups
    deallocate (aa, bb)
    deallocate (g)
    deallocate (work)
#endif
  end subroutine sp95_sha_gg
  !============================================================================
  subroutine sp95_sha (f, a, b)
    !--------------------------------------------------------------
    ! Calculate spherical harmonic analysis of a scalar function f
    ! on a regular grid.
    ! Output: a,b = real & imaginary parts of spectral coefficients
    !         Coefficients(m,n) of a, b: m=1...mmax, n=m..nlat,
    !         where m is the maximum(+1) longitudinal wavenumber.
    ! (For the normalization used c.f. sections 4 and 5 of the
    ! SPHEREPACK manual.)
    !--------------------------------------------------------------
    real(wp), intent(in)  :: f(:,:)
    real(wp), intent(out) :: a(:,:), b(:,:)

#ifndef SPHEREPACK
call finish ('sp95_sha','libspherepack not linked')
#else
    !----------------
    ! Local variables
    !----------------
    integer :: nx, ny, nt, isym, ierror, l1, l2, lshaec
    integer :: lwork, idg, jdg, mdab, ndab, mmax
    logical :: need_init
    real(wp), allocatable, dimension(:)     :: work
    real(wp), allocatable, dimension(:,:,:) :: g, aa, bb

    nx = size (f, dim=1)
    ny = size (f, dim=2)
    nt = 1
    isym = 0
    if (size (a,dim=2) /= ny .or. size (b,dim=2) /= ny .or. &
         any (shape (a) /= shape (b))) then
       write (0,*) "shape (f) =", shape (f)
       write (0,*) "shape (a) =", shape (a)
       write (0,*) "shape (b) =", shape (b)
       call finish ("sp95_sha", "incompatible shapes")
    end if

    need_init = .false.
    if (nx /= nlon .or. ny /= nlat) then
       need_init = .true.
       nlon = nx
       nlat = ny
    end if

    if (allocated (wshaec)) then
       if (need_init) deallocate (wshaec)
    end if
    call set_wshaec (nlat, nlon)
    lshaec = size (wshaec)

    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lwork = nlat*(nt*nlon + max (3*l2,nlon))    ! Size of lwork for shaec
    allocate (work(lwork))

    !------------------------------------
    ! Maximum(+1) longitudinal wavenumber
    !------------------------------------
    mmax = min (l1, size (a,dim=1), size(b,dim=1))

    ! Copy field to Spherepack's internal representation (mathematical coord.)
    idg = nlat
    jdg = nlon
    allocate (g(idg,jdg,1))
    call sp95_geo2math_scalar_2d (phi_geo=f, phi_math=g(:,:,1))

    mdab = l1
    ndab = nlat
    allocate (aa(mdab,ndab,1))
    allocate (bb(mdab,ndab,1))

    ! Spherical harmonic analysis of a scalar field
    call shaec (nlat, nlon, isym, nt, g, idg, jdg, aa, bb, mdab, ndab, &
                wshaec, lshaec, work, lwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_sha", "failure in shagc")
    end if

    a(1:mmax, :) = aa(1:mmax,:,1)
    a(mmax+1:,:) = 0
    b(1:mmax, :) = bb(1:mmax,:,1)
    b(mmax+1:,:) = 0

    ! Final cleanups
    deallocate (aa, bb)
    deallocate (g)
    deallocate (work)
#endif
  end subroutine sp95_sha
  !============================================================================
  subroutine sp95_shs_gg (f, a, b)
    !--------------------------------------------------------------
    ! Calculate spherical harmonic synthesis of a scalar function f
    ! on a Gaussian grid.
    ! Input: a,b = real & imaginary parts of spectral coefficients.
    !        Coefficients(m,n) of a, b: m=1...mmax, n=m..nlat,
    !        where m is the maximum(+1) longitudinal wavenumber.
    !--------------------------------------------------------------
    real(wp), intent(out) :: f(:,:)
    real(wp), intent(in)  :: a(:,:), b(:,:)

#ifndef SPHEREPACK
call finish ('sp95_shs_gg','libspherepack not linked')
#else
    !----------------
    ! Local variables
    !----------------
    integer :: nx, ny, nt, isym, ierror, l1, l2, lshsgc
    integer :: lwork, idg, jdg, mdab, ndab, mmax
    logical :: need_init
    real(wp), allocatable, dimension(:)     :: work
    real(wp), allocatable, dimension(:,:,:) :: g, aa, bb

    nx = size (f, dim=1)
    ny = size (f, dim=2)
    nt = 1
    isym = 0
    if (size (a,dim=2) /= ny .or. size (b,dim=2) /= ny .or. &
         any (shape (a) /= shape (b))) then
       write (0,*) "shape (f) =", shape (f)
       write (0,*) "shape (a) =", shape (a)
       write (0,*) "shape (b) =", shape (b)
       call finish ("sp95_shs_gg", "incompatible shapes")
    end if

    need_init = .false.
    if (nx /= nlon .or. ny /= nlat) then
       need_init = .true.
       nlon = nx
       nlat = ny
    end if

    if (allocated (wshsgc)) then
       if (need_init) deallocate (wshsgc)
    end if
    call set_wshsgc (nlat, nlon)
    lshsgc = size (wshsgc)

    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lwork = nlat*(nt*nlon + max (3*l2,nlon))    ! Size of lwork for shsgc
    allocate (work(lwork))

    !------------------------------------
    ! Maximum(+1) longitudinal wavenumber
    !------------------------------------
    mmax = min (l1, size (a,dim=1), size(b,dim=1))

    ! Transfer fields to Spherepack's representation, pad with zeros
    mdab = l1
    ndab = nlat
    allocate (aa(mdab,ndab,1))
    allocate (bb(mdab,ndab,1))
    aa(1:mmax, :,1) = a(1:mmax,:)
    aa(mmax+1:,:,1) = 0
    bb(1:mmax, :,1) = b(1:mmax,:)
    bb(mmax+1:,:,1) = 0
    idg = nlat
    jdg = nlon
    allocate (g(idg,jdg,1))

    ! Spherical harmonic synthesis of a scalar field
    call shsgc (nlat, nlon, isym, nt, g, idg, jdg, aa, bb, mdab, ndab, &
                wshsgc, lshsgc, work, lwork, ierror)
    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_shs_gg", "failure in shsgc")
    end if

    ! Copy field from Spherepack's internal representation
    call sp95_math2geo_scalar_2d (phi_math=g(:,:,1), phi_geo=f)

    ! Final cleanups
    deallocate (aa, bb, g)
    deallocate (work)
#endif
  end subroutine sp95_shs_gg
  !============================================================================
  subroutine sp95_shs (f, a, b)
    !--------------------------------------------------------------
    ! Calculate spherical harmonic synthesis of a scalar function f
    ! on a regular grid.
    ! Input: a,b = real & imaginary parts of spectral coefficients.
    !        Coefficients(m,n) of a, b: m=1...mmax, n=m..nlat,
    !        where m is the maximum(+1) longitudinal wavenumber.
    !--------------------------------------------------------------
    real(wp), intent(out) :: f(:,:)
    real(wp), intent(in)  :: a(:,:), b(:,:)

#ifndef SPHEREPACK
call finish ('sp95_shs','libspherepack not linked')
#else
    !----------------
    ! Local variables
    !----------------
    integer :: nx, ny, nt, isym, ierror, l1, l2, lshsec
    integer :: lwork, idg, jdg, mdab, ndab, mmax
    logical :: need_init
    real(wp), allocatable, dimension(:)     :: work
    real(wp), allocatable, dimension(:,:,:) :: g, aa, bb

    nx = size (f, dim=1)
    ny = size (f, dim=2)
    nt = 1
    isym = 0
    if (size (a,dim=2) /= ny .or. size (b,dim=2) /= ny .or. &
         any (shape (a) /= shape (b))) then
       write (0,*) "shape (f) =", shape (f)
       write (0,*) "shape (a) =", shape (a)
       write (0,*) "shape (b) =", shape (b)
       call finish ("sp95_shs", "incompatible shapes")
    end if

    need_init = .false.
    if (nx /= nlon .or. ny /= nlat) then
       need_init = .true.
       nlon = nx
       nlat = ny
    end if

    if (allocated (wshsec)) then
       if (need_init) deallocate (wshsec)
    end if
    call set_wshsec (nlat, nlon)
    lshsec = size (wshsec)

    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lwork = nlat*(nt*nlon + max (3*l2,nlon))    ! Size of lwork for shsec
    allocate (work(lwork))

    !------------------------------------
    ! Maximum(+1) longitudinal wavenumber
    !------------------------------------
    mmax = min (l1, size (a,dim=1), size(b,dim=1))

    ! Transfer fields to Spherepack's representation, pad with zeros
    mdab = l1
    ndab = nlat
    allocate (aa(mdab,ndab,1))
    allocate (bb(mdab,ndab,1))
    aa(1:mmax, :,1) = a(1:mmax,:)
    aa(mmax+1:,:,1) = 0
    bb(1:mmax, :,1) = b(1:mmax,:)
    bb(mmax+1:,:,1) = 0
    idg = nlat
    jdg = nlon
    allocate (g(idg,jdg,1))

    ! Spherical harmonic synthesis of a scalar field
    call shsec (nlat, nlon, isym, nt, g, idg, jdg, aa, bb, mdab, ndab, &
                wshsec, lshsec, work, lwork, ierror)
    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_shs", "failure in shsec")
    end if

    ! Copy field from Spherepack's internal representation
    call sp95_math2geo_scalar_2d (phi_math=g(:,:,1), phi_geo=f)

    ! Final cleanups
    deallocate (aa, bb, g)
    deallocate (work)
#endif
  end subroutine sp95_shs
  !============================================================================
  subroutine sp95_shp_gg (x, mtrunc, ntrunc, zeromean)
    !------------------------------------------------------------
    ! Calculate spherical harmonic projection of a scalar field x
    ! on a Gaussian grid, performed as a spectral analysis, an
    ! (optional) truncation, and a subsequent synthesis.
    !------------------------------------------------------------
    real(wp), intent(inout)        :: x(:,:,:)  ! Field to project
    integer,  intent(in), optional :: mtrunc    ! Truncation (zonal)
    integer,  intent(in), optional :: ntrunc    ! Truncation (meridional)
    logical,  intent(in), optional :: zeromean  ! Remove mean

#ifndef SPHEREPACK
call finish ('sp95_shp_gg','libspherepack not linked')
#else
    !----------------
    ! Local variables
    !----------------
    integer :: nx, ny, nt, isym, ierror, l1, l2, lshagc, lshsgc
    integer :: lwork, idg, jdg, mdab, ndab, mmax
    logical :: need_init
    real(wp), allocatable, dimension(:)     :: work
    real(wp), allocatable, dimension(:,:,:) :: g, aa, bb

    nx = size (x, dim=1)
    ny = size (x, dim=2)
    nt = size (x, dim=3)
    isym = 0

    need_init = .false.
    if (nx /= nlon .or. ny /= nlat) then
       need_init = .true.
       nlon = nx
       nlat = ny
    end if

    if (allocated (wshagc)) then
       if (need_init) deallocate (wshagc)
    end if
    call set_wshagc (nlat, nlon)
    lshagc = size (wshagc)

    if (allocated (wshsgc)) then
       if (need_init) deallocate (wshsgc)
    end if
    call set_wshsgc (nlat, nlon)
    lshsgc = size (wshsgc)

    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lwork = nlat*(nt*nlon + max (3*l2,nlon))    ! Size of lwork for shsgc
    allocate (work(lwork))

    !------------------------------------
    ! Maximum(+1) longitudinal wavenumber
    !------------------------------------
    mmax = min (nlat, (nlon+1)/2)
    if (present (mtrunc)) mmax = min (mmax, mtrunc)
    if (present (ntrunc)) mmax = min (mmax, ntrunc)

    ! Copy field to Spherepack's internal representation (mathematical coord.)
    idg = nlat
    jdg = nlon
    allocate (g(idg,jdg,nt))
    call sp95_geo2math_scalar (phi_geo=x, phi_math=g)

    mdab = l1
    ndab = nlat
    allocate (aa(mdab,ndab,nt))
    allocate (bb(mdab,ndab,nt))

    ! Spherical harmonic analysis of a scalar field
    call shagc (nlat, nlon, isym, nt, g, idg, jdg, aa, bb, mdab, ndab, &
                wshagc, lshagc, work, lwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_shp_gg", "failure in shagc")
    end if

    if (present (zeromean)) then
       if (zeromean) then
          aa(:,1,:) = 0
          bb(:,1,:) = 0
       end if
    end if
    ! Truncation
    aa(mmax+1:,:,:) = 0
    bb(mmax+1:,:,:) = 0
    if (present (ntrunc)) then
       aa(:,ntrunc+1:,:) = 0
       bb(:,ntrunc+1:,:) = 0
    end if

    ! Spherical harmonic synthesis of a scalar field
    call shsgc (nlat, nlon, isym, nt, g, idg, jdg, aa, bb, mdab, ndab, &
                wshsgc, lshsgc, work, lwork, ierror)
    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_shp_gg", "failure in shsgc")
    end if

    ! Copy field from Spherepack's internal representation
    call sp95_math2geo_scalar (phi_math=g, phi_geo=x)

    ! Final cleanups
    deallocate (aa, bb, g)
    deallocate (work)
#endif

  end subroutine sp95_shp_gg
  !============================================================================
  subroutine sp95_slap_gg (x, Rsphere)
    !------------------------------------------------------------
    ! Calculate Laplacian of a scalar field x on a Gaussian grid.
    !------------------------------------------------------------
    real(wp), intent(inout) :: x(:,:,:)         ! Scalar field
    real(wp), intent(in)    :: Rsphere          ! Radius of sphere

#ifndef SPHEREPACK
call finish ('sp95_slap_gg','libspherepack not linked')
#else
    !----------------
    ! Local variables
    !----------------
    integer  :: nx, ny, nt, isym, ierror, l1, l2, lshagc, lshsgc
    integer  :: lwork, idg, jdg, mdab, ndab, n
    logical  :: need_init
    real(wp) :: fac
    real(wp), allocatable, dimension(:)     :: work
    real(wp), allocatable, dimension(:,:,:) :: g, aa, bb

    nx = size (x, dim=1)
    ny = size (x, dim=2)
    nt = size (x, dim=3)
    isym = 0

    need_init = .false.
    if (nx /= nlon .or. ny /= nlat) then
       need_init = .true.
       nlon = nx
       nlat = ny
    end if

    if (allocated (wshagc)) then
       if (need_init) deallocate (wshagc)
    end if
    call set_wshagc (nlat, nlon)
    lshagc = size (wshagc)

    if (allocated (wshsgc)) then
       if (need_init) deallocate (wshsgc)
    end if
    call set_wshsgc (nlat, nlon)
    lshsgc = size (wshsgc)

    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lwork = nlat*(nt*nlon + max (3*l2,nlon))    ! Size of lwork for shagc,shsgc
    allocate (work(lwork))

    ! Copy field to Spherepack's internal representation (mathematical coord.)
    idg = nlat
    jdg = nlon
    allocate (g(idg,jdg,nt))
    call sp95_geo2math_scalar (phi_geo=x, phi_math=g)

    mdab = l1
    ndab = nlat
    allocate (aa(mdab,ndab,nt))
    allocate (bb(mdab,ndab,nt))

    ! Spherical harmonic analysis of a scalar field
    call shagc (nlat, nlon, isym, nt, g, idg, jdg, aa, bb, mdab, ndab, &
                wshagc, lshagc, work, lwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_slap_gg", "failure in shagc")
    end if

    ! Apply Laplacian
    do n = 1, nlat
       fac = - (n*(n-1)) / Rsphere**2
       aa(:n,n,:) = aa(:n,n,:) * fac
       bb(:n,n,:) = bb(:n,n,:) * fac
    end do

    ! Spherical harmonic synthesis of a scalar field
    call shsgc (nlat, nlon, isym, nt, g, idg, jdg, aa, bb, mdab, ndab, &
                wshsgc, lshsgc, work, lwork, ierror)
    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("sp95_slap_gg", "failure in shsgc")
    end if

    ! Copy field from Spherepack's internal representation
    call sp95_math2geo_scalar (phi_math=g, phi_geo=x)

    ! Final cleanups
    deallocate (aa, bb, g)
    deallocate (work)
#endif

  end subroutine sp95_slap_gg
  !============================================================================

#ifdef SPHEREPACK
  !============================================================================
  ! Internal routines for initializations that call Spherepack
  !============================================================================
  subroutine set_wshagc (nlat, nlon)
    !------------------------------------------------------------
    ! Initialize coefficients wshagc for scalar harmonic analysis
    !------------------------------------------------------------
    integer, intent(in)  :: nlat, nlon
    !----------------
    ! Local variables
    !----------------
    integer                       :: l1, l2, ldwork, ierror, lshagc
    double precision, allocatable :: dwork(:)

    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lshagc = nlat*(3*l1+2*l2-2) + 3*l1*(1-l1)/2 + nlon + 15

    if (allocated (wshagc)) then
       if (size (wshagc) == lshagc) return      ! Assume properly initialized
       deallocate (wshagc)                      ! Wrong size, needs realloc.
    end if
    allocate (wshagc(lshagc))

    ldwork = nlat * (nlat+4)
    allocate (dwork(ldwork))

    call shagci (nlat, nlon, wshagc, lshagc, dwork, ldwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("set_wshagc", "failure in shagci")
    end if
    deallocate (dwork)
  end subroutine set_wshagc
  !============================================================================
  subroutine set_wshaec (nlat, nlon)
    !------------------------------------------------------------
    ! Initialize coefficients wshaec for scalar harmonic analysis
    !------------------------------------------------------------
    integer, intent(in)  :: nlat, nlon
    !----------------
    ! Local variables
    !----------------
    integer                       :: l1, l2, ldwork, ierror, lshaec
    double precision, allocatable :: dwork(:)

    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lshaec = 2*nlat*l2 + 3*(l1-2)*(2*nlat-l1-1)/2 + nlon + 15

    if (allocated (wshaec)) then
       if (size (wshaec) == lshaec) return      ! Assume properly initialized
       deallocate (wshaec)                      ! Wrong size, needs realloc.
    end if
    allocate (wshaec(lshaec))

    ldwork = nlat+1
    allocate (dwork(ldwork))

    call shaeci (nlat, nlon, wshaec, lshaec, dwork, ldwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("set_wshaec", "failure in shaeci")
    end if
    deallocate (dwork)
  end subroutine set_wshaec
  !============================================================================
  subroutine set_wshsgc (nlat, nlon)
    !-------------------------------------------------------------
    ! Initialize coefficients wshsgc for scalar harmonic synthesis
    !-------------------------------------------------------------
    integer, intent(in)  :: nlat, nlon
    !----------------
    ! Local variables
    !----------------
    integer                       :: l1, l2, ldwork, ierror, lshsgc
    double precision, allocatable :: dwork(:)

    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lshsgc = nlat*(3*l1+2*l2-2) + 3*l1*(1-l1)/2 + nlon + 15

    if (allocated (wshsgc)) then
       if (size (wshsgc) == lshsgc) return      ! Assume properly initialized
       deallocate (wshsgc)                      ! Wrong size, needs realloc.
    end if
    allocate (wshsgc(lshsgc))

    ldwork = nlat * (nlat+4)
    allocate (dwork(ldwork))

    call shsgci (nlat, nlon, wshsgc, lshsgc, dwork, ldwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("set_wshsgc", "failure in shsgci")
    end if
    deallocate (dwork)
  end subroutine set_wshsgc
  !============================================================================
  subroutine set_wshsec (nlat, nlon)
    !-------------------------------------------------------------
    ! Initialize coefficients wshsec for scalar harmonic synthesis
    !-------------------------------------------------------------
    integer, intent(in)  :: nlat, nlon
    !----------------
    ! Local variables
    !----------------
    integer                       :: l1, l2, ldwork, ierror, lshsec
    double precision, allocatable :: dwork(:)

    l1 = min (nlat, (nlon+2)/2)
    l2 = (nlat+1)/2
    lshsec = 2*nlat*l2 + 3*(l1-2)*(nlat+nlat-l1-1)/2 + nlon + 15

    if (allocated (wshsec)) then
       if (size (wshsec) == lshsec) return      ! Assume properly initialized
       deallocate (wshsec)                      ! Wrong size, needs realloc.
    end if
    allocate (wshsec(lshsec))

    ldwork = nlat+1
    allocate (dwork(ldwork))

    call shseci (nlat, nlon, wshsec, lshsec, dwork, ldwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("set_wshsec", "failure in shseci")
    end if
    deallocate (dwork)
  end subroutine set_wshsec
  !============================================================================
  subroutine set_wvhagc (nlat, nlon)
    !------------------------------------------------------------
    ! Initialize coefficients wvhagc for vector harmonic analysis
    !------------------------------------------------------------
    integer, intent(in)  :: nlat, nlon
    !----------------
    ! Local variables
    !----------------
    integer                       :: l1, l2, ldwork, ierror, lvhagc
    double precision, allocatable :: dwork(:)

    l1 = min (nlat, (nlon+1)/2)
    l2 = (nlat+1)/2
    lvhagc = 4*nlat*l2+3*max (l1-2,0)*(2*nlat-l1-1)+nlon+l2+15

    if (allocated (wvhagc)) then
       if (size (wvhagc) == lvhagc) return      ! Assume properly initialized
       deallocate (wvhagc)                      ! Wrong size, needs realloc.
    end if
    allocate (wvhagc(lvhagc))

    ldwork = 2*nlat*(nlat+1)+1
    allocate (dwork(ldwork))

    call vhagci (nlat, nlon, wvhagc, lvhagc, dwork, ldwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("set_wvhagc", "failure in vhagci")
    end if
    deallocate (dwork)
  end subroutine set_wvhagc
  !============================================================================
  subroutine set_wvhsgc (nlat, nlon)
    !-------------------------------------------------------------
    ! Initialize coefficients wvhsgc for vector harmonic synthesis
    !-------------------------------------------------------------
    integer, intent(in)  :: nlat, nlon
    !----------------
    ! Local variables
    !----------------
    integer                       :: l1, l2, ldwork, ierror, lvhsgc
    double precision, allocatable :: dwork(:)

    l1 = min (nlat, (nlon+1)/2)
    l2 = (nlat+1)/2
    lvhsgc = 4*nlat*l2+3*max (l1-2,0)*(2*nlat-l1-1)+nlon+l2+15

    if (allocated (wvhsgc)) then
       if (size (wvhsgc) == lvhsgc) return      ! Assume properly initialized
       deallocate (wvhsgc)                      ! Wrong size, needs realloc.
    end if
    allocate (wvhsgc(lvhsgc))

    ldwork = 2*nlat*(nlat+1)+1
    allocate (dwork(ldwork))

    call vhsgci (nlat, nlon, wvhsgc, lvhsgc, dwork, ldwork, ierror)

    if (ierror /= 0) then
       write (0,*) "ierror =", ierror
       call finish ("set_wvhsgc", "failure in vhsgci")
    end if
    deallocate (dwork)
  end subroutine set_wvhsgc
  !============================================================================
#endif

  !=====================================================================
  !
  ! Conversion of grids between the geophysical convention used at DWD
  ! (south pole at (1:nlon,1), north pole at (1:nlon,nlat))
  ! and the mathematical convention applied in spherepack
  ! (north pole at (1,1:nlon), south pole at (nlat,1:nlon)).
  !
  ! Mapping of a scalar field:        g_math(i,j) <-> g_geo(j,nlat+1-i)
  ! Vector field (additionally): (v_math,w_math)) <-> (-v_geo,u_geo)
  !
  ! See the description of the routines geo2math[sv] and math2geo[sv]
  ! for more information on the grid types and conventions.
  ! The implementation below does not need working storage but requires
  ! that input and output arrays do not overlap.  It also allows a
  ! rescaling transformation for scalar fields.
  !
  !=====================================================================

  subroutine sp95_geo2math_scalar (phi_geo, phi_math)
    real(wp), intent(in)  :: phi_geo(:,:,:)
    real(wp), intent(out) :: phi_math(:,:,:)

    integer :: nlon, nlat, nt, i

    nlon = size (phi_geo, dim=1)
    nlat = size (phi_geo, dim=2)
    nt   = size (phi_geo, dim=3)
    if ( size (phi_math, dim=1) /= nlat .or. &
         size (phi_math, dim=2) /= nlon .or. &
         size (phi_math, dim=3) /= nt) then
       write (0,*) shape (phi_geo), shape (phi_math)
       call finish ("sp95_geo2math_scalar", "shapes do not match!")
    end if

    do i = 1, nlat
       phi_math(i,1:nlon,1:nt) = phi_geo(1:nlon,(nlat+1-i),1:nt)
    end do
  end subroutine sp95_geo2math_scalar

  ! --

  subroutine sp95_math2geo_scalar (phi_math, phi_geo)
    real(wp), intent(in)  :: phi_math(:,:,:)
    real(wp), intent(out) :: phi_geo(:,:,:)

    integer :: nlon, nlat, nt, i

    nlon = size (phi_geo, dim=1)
    nlat = size (phi_geo, dim=2)
    nt   = size (phi_geo, dim=3)
    if ( size (phi_math, dim=1) /= nlat .or. &
         size (phi_math, dim=2) /= nlon .or. &
         size (phi_math, dim=3) /= nt) then
       write (0,*) shape (phi_geo), shape (phi_math)
       call finish ("sp95_math2geo_scalar", "shapes do not match!")
    end if

    do i = 1, nlat
       phi_geo(1:nlon,(nlat+1-i),1:nt) = phi_math(i,1:nlon,1:nt)
    end do
  end subroutine sp95_math2geo_scalar

  ! --

  subroutine sp95_geo2math_scaled (phi_geo, phi_math, scale)
    real(wp), intent(in)  :: phi_geo(:,:,:)
    real(wp), intent(out) :: phi_math(:,:,:)
    real(wp), intent(in)  :: scale

    integer :: nlon, nlat, nt, i

    nlon = size (phi_geo, dim=1)
    nlat = size (phi_geo, dim=2)
    nt   = size (phi_geo, dim=3)
    if ( size (phi_math, dim=1) /= nlat .or. &
         size (phi_math, dim=2) /= nlon .or. &
         size (phi_math, dim=3) /= nt) then
       write (0,*) shape (phi_geo), shape (phi_math)
       call finish ("sp95_geo2math_scaled", "shapes do not match!")
    end if

    do i = 1, nlat
       phi_math(i,1:nlon,1:nt) = scale * phi_geo(1:nlon,(nlat+1-i),1:nt)
    end do
  end subroutine sp95_geo2math_scaled

  ! --

  subroutine sp95_math2geo_scaled (phi_math, phi_geo, scale, ig)
    real(wp), intent(in)  :: phi_math(:,:,:)
    real(wp), intent(out) :: phi_geo(:,:,:)
    real(wp), intent(in)  :: scale
    integer,  intent(in)  :: ig

    integer :: nlon, nlat, nt, i

    nlon = size (phi_geo, dim=1)
    nlat = size (phi_geo, dim=2)
    nt   = size (phi_geo, dim=3)
    if ( size (phi_math, dim=1) /= nlat .or. &
         size (phi_math, dim=2) /= nlon .or. &
         size (phi_math, dim=3) /= nt) then
       write (0,*) shape (phi_geo), shape (phi_math)
       call finish ("sp95_math2geo_scaled", "shapes do not match!")
    end if

    if (ig==0) then
      do i = 1, nlat
         phi_geo(1:nlon,(nlat+1-i),1:nt) = scale * phi_math(i,1:nlon,1:nt)
      end do
    else
      do i = 1, nlat
         phi_geo(1:nlon,i,1:nt) = scale * phi_math(i,1:nlon,1:nt)
      end do
    endif

  end subroutine sp95_math2geo_scaled

  ! --

  subroutine sp95_geo2math_vector (u_geo, v_geo, v_math, w_math, ig)
    real(wp), dimension(:,:,:), intent(in)  :: u_geo, v_geo
    real(wp), dimension(:,:,:), intent(out) :: v_math, w_math
    integer,                    intent(in)  :: ig
    integer :: nlon, nlat, nt, i

    nlon = size (v_geo, dim=1)
    nlat = size (v_geo, dim=2)
    nt   = size (v_geo, dim=3)
    if ( size (v_math, dim=1) /= nlat .or. &
         size (v_math, dim=2) /= nlon .or. &
         size (v_math, dim=3) /= nt   .or. &
         any (shape (u_geo)  /= shape (v_geo)) .or. &
         any (shape (v_math) /= shape (w_math)) ) then
       write (0,*) shape (u_geo), shape (v_geo), shape (v_math), shape (w_math)
       call finish ("sp95_math2geo_scaled", "shapes do not match!")
    end if

    if (ig==0) then
      do i = 1, nlat
         w_math(i,1:nlon,1:nt) =  u_geo(1:nlon,(nlat+1-i),1:nt)
         v_math(i,1:nlon,1:nt) = -v_geo(1:nlon,(nlat+1-i),1:nt)
      end do
    else
      do i = 1, nlat
         w_math(i,1:nlon,1:nt) =  u_geo(1:nlon,i,1:nt)
         v_math(i,1:nlon,1:nt) = -v_geo(1:nlon,i,1:nt)
      end do
    endif

  end subroutine sp95_geo2math_vector

  ! --

  subroutine sp95_math2geo_vector (v_math, w_math, u_geo, v_geo)
    real(wp), dimension(:,:,:), intent(in)  :: v_math, w_math
    real(wp), dimension(:,:,:), intent(out) :: u_geo, v_geo

    integer :: nlon, nlat, nt, i

    nlon = size (v_geo, dim=1)
    nlat = size (v_geo, dim=2)
    nt   = size (v_geo, dim=3)
    if ( size (v_math, dim=1) /= nlat .or. &
         size (v_math, dim=2) /= nlon .or. &
         size (v_math, dim=3) /= nt   .or. &
         any (shape (u_geo)  /= shape (v_geo)) .or. &
         any (shape (v_math) /= shape (w_math)) ) then
       write (0,*) shape (u_geo), shape (v_geo), shape (v_math), shape (w_math)
       call finish ("sp95_math2geo_vector", "shapes do not match!")
    end if

    do i = 1, nlat
       u_geo(1:nlon,(nlat+1-i),1:nt) =  w_math(i,1:nlon,1:nt)
       v_geo(1:nlon,(nlat+1-i),1:nt) = -v_math(i,1:nlon,1:nt)
    end do
  end subroutine sp95_math2geo_vector

  ! --

  subroutine sp95_geo2math_scalar_2d (phi_geo, phi_math)
    real(wp), intent(in)  :: phi_geo(:,:)
    real(wp), intent(out) :: phi_math(:,:)

    integer :: nlon, nlat, i

    nlon = size (phi_geo, dim=1)
    nlat = size (phi_geo, dim=2)
    if ( size (phi_math, dim=1) /= nlat .or. &
         size (phi_math, dim=2) /= nlon) then
       write (0,*) shape (phi_geo), "/=", shape (phi_math)
       call finish ("sp95_geo2math_scalar_2d", "shapes do not match!")
    end if

    do i = 1, nlat
       phi_math(i,1:nlon) = phi_geo(1:nlon,(nlat+1-i))
    end do
  end subroutine sp95_geo2math_scalar_2d

  ! --

  subroutine sp95_math2geo_scalar_2d (phi_math, phi_geo)
    real(wp), intent(in)  :: phi_math(:,:)
    real(wp), intent(out) :: phi_geo(:,:)

    integer :: nlon, nlat, i

    nlon = size (phi_geo, dim=1)
    nlat = size (phi_geo, dim=2)
    if ( size (phi_math, dim=1) /= nlat .or. &
         size (phi_math, dim=2) /= nlon) then
       write (0,*) shape (phi_geo), shape (phi_math)
       call finish ("sp95_math2geo_scalar_2d", "shapes do not match!")
    end if

    do i = 1, nlat
       phi_geo(1:nlon,(nlat+1-i)) = phi_math(i,1:nlon)
    end do
  end subroutine sp95_math2geo_scalar_2d

  ! --

  subroutine sp95_math2geo_vector_2d (v_math, w_math, u_geo, v_geo)
    real(wp), dimension(:,:), intent(in)  :: v_math, w_math
    real(wp), dimension(:,:), intent(out) :: u_geo, v_geo

    integer :: nlon, nlat, i

    nlon = size (v_geo, dim=1)
    nlat = size (v_geo, dim=2)
    if ( size (v_math, dim=1) /= nlat .or. &
         size (v_math, dim=2) /= nlon .or. &
         any (shape (u_geo)  /= shape (v_geo)) .or. &
         any (shape (v_math) /= shape (w_math)) ) then
       write (0,*) shape (u_geo), shape (v_geo), shape (v_math), shape (w_math)
       call finish ("sp95_math2geo_vector_2d", "shapes do not match!")
    end if

    do i = 1, nlat
       u_geo(1:nlon,(nlat+1-i)) =  w_math(i,1:nlon)
       v_geo(1:nlon,(nlat+1-i)) = -v_math(i,1:nlon)
    end do
  end subroutine sp95_math2geo_vector_2d

end module spherepack95

!- End of module spherepack95
