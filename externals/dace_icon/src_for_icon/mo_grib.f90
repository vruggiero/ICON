!
!+ Routines to read and write GRIB files
!
MODULE mo_grib
!
! Description:
!   Routines to read and write GRIB files
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
! V1_2         2008/12/04 Andreas Rhodin
!  read_atm_state: new optional parameters: optionals, unsp_type
! V1_4         2009/03/26 Andreas Rhodin
!  changes for COSMO; optimisations in read_multi (bcast/scatter)
! V1_5         2009/05/25 Andreas Rhodin
!  read meta data from first valid GRIB field
! V1_7         2009/08/24 Andreas Rhodin
!  Changes for COSMO/LETKF
! V1_8         2009/12/09 Andreas Rhodin
!  read_atm_grid: changed error handling for subroutine get_inventory
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changes for GRIB2 API; ensemble IO; COSMO vertical coordinates
! V1_15        2011/12/06 Harald Anlauf
!  formatted printing of some diagnostic output
! V1_19        2012-04-16 Andreas Rhodin
!  changes for GME LETKF (set Zusatzelementnummer for GME W_SO to 255)
!  enhanced communication for LETKF parallel I/O
! V1_20        2012-06-18 Andreas Rhodin
!  fix parallal ensemble input; initialize variables before first use
! V1_22        2013-02-13 Harald Anlauf
!  adapt to ICON
! V1_23        2013-03-26 Harald Anlauf
!  write_atm_gen: set GRIB1/2 codes for all fields to be written
! V1_26        2013/06/27 Harald Anlauf
!  Changes for GRIB2/GRIB_API and ICON
! V1_27        2013-11-08 Harald Anlauf
!  changes for GRIB2 API and ICON
! V1_28        2014/02/26 Harald Anlauf
!  Fixes for LETKF; GRIB_API, GRIB1 and GRIB2, IFS
! V1_31        2014-08-21 Andreas Rhodin
!  new parameter 'leveltypes' to subroutine read (atmospheric state from GRIB)
! V1_37        2014-12-23 Harald Anlauf
!  read_atm_grid: enforce unique and consistent numberOfGridUsed
!  write_var_gen: implement blocking of gather operation to limit memory usage
! V1_42        2015-06-08 Andreas Rhodin
!  LETKF coarse grid; ICON local patch
! V1_43        2015-08-19 Andreas Rhodin
!  handle 'fr_lake','depth_lk'.
! V1_44        2015-09-30 Andreas Rhodin
!  read_atm_grid: optional parameter 'hhl'
!  LETKF: option to handle missing members
! V1_45        2015-12-15 Harald Anlauf
!  read_atm_ens/write_atm_ens: optimize I/O processor pattern for ensembles;
!  Cleanup towards TR15581/F2003 compatibility; use OpenMP workshare;
!  Skip isobaric layers when deriving vertical grid
! V1_46        2016-02-05 Andreas Rhodin
!  revise setup of ivctype, refatm; read_atm_grid: new parameter 'optionals'
! V1_47        2016-06-06 Harald Anlauf
!  Enable GRIB encoding of selected analysis fields with 24 bits.
!  subroutine read_atm_ens, read_atm_state: new optional parameter 'ierr'.
!  new optional parameter 'pos' to subroutine read_single.
!  write COSMO vcp to GRIB, changes for refatm=2,add 'u','v','w' (COMET).
! V1_48        2016-10-06 Harald Anlauf
!  estimate unperturbed model level heights for gen.vert.coordinate (ICON)
!  support COSMO ivctype=3,4 (sleeve coordinates), reference atmosphere 2
! V1_49        2016-10-25 Harald Anlauf
!  Generalized vertical coordinate: use nlev from GRIB2
! V1_50        2017-01-09 Harald Anlauf
!  Scaling factor for ECMWF grib edition 1 fields from mars
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Author:
! Andreas Rhodin  DWD  2002-2008  original source
! Harald Anlauf   DWD  2007       extend to (ECMWF) ensemble members
!==============================================================================
#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION) && 0
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!==============================================================================
#include "tr15581.incf"
  !-------------
  ! Modules used
  !-------------
  use mo_kind,          only : wp, dp, i8         ! kind parameters
  use mo_exception,     only : finish, message    ! abort routine
  use mo_dace_string,   only : char3,            &! integer -> char(3)
                               split,            &! put words into array
                               concat,           &! concat array to string
                               eval_string,      &! add,remove,... from string
                               byte2hex,         &! CHAR to 'hexadecimal'
                               toupper            ! convert to uppercase
  use mo_memory,        only : t_m,              &! type to hold 3d fields
                               t_mi               ! field metadata type
  use mo_atm_grid,      only : t_grid,           &! grid derived type
                               t_vcoord,         &! vertical coordinates type
                               t_refatm,         &! reference atmosphere type
                               construct,        &! constructor routine
                               destruct,         &! destructor  routine
                               allocate,         &! alloc. pointer components
                               deallocate,       &! dealloc. pointer components
                               set_ptopf,        &! set topmost level pressure
                               set_plev_indices, &! set some level indices
                               set_zlev_ref,     &! estimate unperturbed levels
                               cosmo_ref_atm,    &! set COSMO reference atm.
                               cosmo_vcp,        &! compute GDS vert.coord.par.
                               get_vertcoord_new,&! COSMO vert.coord.parameters
                               IVCTYPE_ICON,     &! ICON vertical coord.type
                               IVCTYPE_GME,      &! GME  vertical coord.type
                               p_bcast,          &! broadcast t_vcoord, t_refatm
                               reduced_grid,     &! derive reduced grid (weights)
                               read_icon_metadata ! setup ICON grid metadata
  use mo_icon_grid,     only : check_neighbours   ! check for correct order
  use mo_atm_state,     only : t_atm,            &! atmosphere derived type
                               allocate,         &! alloc. pointer components
                               set_pointers,     &! update pointers
                               deallocate         ! deallocate pointers
  use mo_grid_intpol,   only : distribute_icon_gme! re-distributed ICON points
  use mo_emos_grib1,    only : t_grib1,          &! GRIB record data type
!                              pbopen,           &! open     GRIB file
!                              pbclose,          &! close    GRIB file
                               pbseek,           &! position GRIB file
!                              pbgrib,           &! read     GRIB record
!                              gribex,           &! decode   GRIB record
                               local_dwd,        &! DWD   extension present?
                               local_ecmwf        ! ECMWF extension present?

  use mo_grib_handling, only : t_inventory,      &! GRIB inventory datatype
                               get_inventory,    &! derive inventory from GRIB
                               print_inventory,  &! print  inventory
                               p_bcast,          &! broadcast inventory via MPI
                               set_defaults,     &! set center, subcenter, process
                               set_ver_time,     &! specify verification time
                               set_ref_time,     &! specify reference time
                               set_dwd_ens,      &! specify DWD   local extens.
                               set_ecmwf_ens,    &! specify ECMWF local extens.
                               set_grid_latlon,  &! specify lat/lon grid
                               set_grid_gauss,   &! specify gaussian grid
                               set_grid_tri,     &! specify triangular grid
                               set_grid_icon,    &! specify ICON grid
                               set_grid_vert,    &! specify leveltype, vcp
                               set_code,         &! specify code, table
                               set_level,        &! specify level
                               set_data,         &! transfer data
                               open_gribfile,    &! GRIB open wrapper
                               write_gribrecord, &!
                               close_gribfile     ! GRIB close wrapper

  use mo_grib_invt,     only : operator(==)       ! compare inventory

  use mo_grib12,        only : read_gribex        ! read + decode GRIB 1 or 2

  use mo_wmo_tables,    only : WMO3_HYBRID,      &!
                               WMO3_HYBRIDB,     &!
                               WMO3_HHYBRID,     &! height-based hybrid levels
                               WMO3_GENV,        &! generalized vertical coord.
                               WMO3_ISOBARIC,    &!
                               WMO3_SURFACE,     &!
                               WMO3_SEALEVEL,    &!
                               WMO3_ABOVESUR,    &!
                               WMO3_BELOWSUR,    &! Levels below surface
                               WMO3_TILE_LAND,   &! grid-tile land fraction
                               WMO6_GAUSSIAN,    &! Gaussian lat/lon grid
                               DWD6_ICOSAHEDRON, &! Icosahedral triangular grid
                               DWD6_ICON,        &! ICON unstructured grid
                               WMO6_LATLON,      &! Latitude/Longitude Grid
                               WMO6_ROTLL,       &! Rotated lat/lon grid
                               WMO8_J_POSITIVE    !
!                              WMO5_FORECAST      !
  use mo_gribtables,    only : ar_des,           &! entry descriptor
                               search             ! search GRIB code
  use mo_mpi_dace,      only : dace,             &! MPI communication info
                               p_parallel,       &! true for parallel run
                               p_bcast,          &! generic broadcast routine
                               p_bcast_ptr,      &! pointer broadcast routine
                               p_min,            &! minimum over PEs
                               self               ! communicator group
  use mo_atm_transp,    only : scatter_level,    &! scatter model levels
                               scatter_multi,    &!++ specific routine for NEC
                               gather_level,     &! gather model level
!                              gather_multi,     &! gather multi level field
                               alltoallv_multi,  &! replace scatter_multi
                               alltoallv_multi_w  ! replace gather_multi
  use mo_time,          only : t_time,           &! time data type
                               init_time,        &! initialisation routine
                               operator(==),     &! time == time
                               operator(/=),     &! time /= time
                               p_bcast,          &! generic broadcast routine
                               cyyyymmddhhmm,    &! derive string from time
                               imm
  use mo_run_params,    only : p_readgrib,       &! PE used to read GRIB
                               model,            &! 'COSMO' or 'GME'
                               path_file,        &! derive path/file.suffix
                               io_max_gather,    &! gather buffer size for I/O
                               io_max_scatter,   &! scatter buffer size for I/O
                               chk_icon_grid      ! check neighbour references
  use mo_physics,       only : gacc,             &! gravity acceleration
                               d2r                ! conversion degree->radian
  use mo_grib12_dwd,    only : t_par_grib1,      &! Derived type: GRIB1 params.
                               t_par_grib2,      &! Derived type: GRIB2 params.
                               empty2,           &! Empty instance: "     "
                               search_grib1,     &! Search param. in GRIB1 tab.
                               search_grib2       ! Search param. in GRIB2 tab.
  use mo_algorithms,    only : index              ! Index an array for sorting
  implicit none

!==============================================================================
  !----------------
  ! public entities
  !----------------
  private
  public :: read                ! read model state or grid information
  public :: read_multi          ! read a multi-level field
  public :: read_single         ! read a single-level field
  public :: SCATTER, BCAST      ! flag values

  public :: write_grib          ! write atmosphere to GRIB file
  public :: set_gribcodes       ! set grib1/grib2 codes of atmospheric fields
  public :: set_iope_ens        ! Set I/O PE pattern for ensemble members
  public :: check_iope_ens      ! Compare I/O PE patterns
  public :: set_bits            ! Set grib bits/value for specified fields

  public :: IE_OK               ! error return values
  public :: IE_MISSING          ! requested field is missing
  public :: IE_AMBIGUOUS        ! requested field is ambiguous
  public :: IE_FILE             ! file not present or empty
  public :: IE_MEMBER           ! requested member is missing
  public :: IE_LEVELS           ! some levels of field are missing

  public :: t_ctr_inv           ! inventory container
  public :: t_ctr_slot          ! slot      container
  public :: construct           ! construct inventory container
  public :: destruct            ! destruct  inventory container
  public :: read_invt           ! fill inventory container
!==============================================================================
  !--------------------
  ! Inventory Container
  !--------------------
  type t_ctr_inv
    character(len=256)              :: name    = ''      ! file name
    logical                         :: global  = .false. ! inventory broadcasted?
    logical                         :: det     = .false. ! deterministic?
    integer                         :: pe      = -1      ! reading pe
    integer                         :: size    = -1      ! inventory size
    type(t_inventory) ,pointer      :: invt (:) => NULL ()
  end type t_ctr_inv

  type t_ctr_slot
    integer                         :: slot    = -1      ! slot
    integer                         :: nslots  = -1      ! number of slots
    integer                         :: nens    = -1      ! ensemble size
    type(t_ctr_inv)   ,pointer      :: c(:) => NULL ()   ! container for member
  end type t_ctr_slot
!==============================================================================
  !-----------
  ! interfaces
  !-----------
  interface read
  !--------------------------------------------------
  ! generic interface to read model or ensemble state
  ! or derive grid description from a GRIB file
  !--------------------------------------------------
    module procedure read_atm_state   ! read atmospheric state
    module procedure read_atm_ens     ! read atmospheric ensemble
    module procedure read_atm_grid    ! read atmospheric grid
  end interface

  interface pprint
     module procedure print_2d_rfield, print_2d_lfield, print_3d_rfield
  end interface

  interface write_grib
  !---------------------------------------------------
  ! generic interface to write model or ensemble state
  ! or grid invariant fields to a GRIB file
  !---------------------------------------------------
    module procedure write_atm_grib   ! write atmospheric state
    module procedure write_atm_ens    ! write atmospheric state ensemble
    module procedure write_atm_grib_p ! write atmospheric states in parallel
    module procedure write_grid_grib  ! write atmospheric grid
    module procedure write_var_grib   ! write atmospheric variable
  end interface write_grib

  interface construct
  !------------------------------
  ! construct inventory container
  !------------------------------
    module procedure construct_ctr_slot
  end interface construct

  interface destruct
  !-----------------------------
  ! destruct inventory container
  !-----------------------------
    module procedure destruct_ctr_slot
  end interface destruct

  interface read_invt
  !--------------------------------
  ! read inventories into container
  !--------------------------------
    module procedure read_inv_ctr_slot
    module procedure read_inv_ctr
    module procedure read_inv_files
  end interface read_invt

!==============================================================================
  !----------
  ! constants
  !----------
  integer ,parameter :: BCAST   = 1
  integer ,parameter :: SCATTER = 2
  !-------------------
  ! error return codes
  !-------------------
  integer ,parameter :: IE_OK        =  0  ! no error
  integer ,parameter :: IE_MISSING   = -1  ! requested field is missing
  integer ,parameter :: IE_AMBIGUOUS = -2  ! requested field is ambiguous
  integer ,parameter :: IE_FILE      = -3  ! file not present or empty
  integer ,parameter :: IE_MEMBER    = -4  ! requested member is missing
  integer ,parameter :: IE_LEVELS    = -5  ! some levels of field are missing
  integer ,parameter :: IE_INV_LEVEL = -6  ! invalid level
!==============================================================================
contains
!==============================================================================
  subroutine read_atm_state (state, file, invt, fields, runtype, time,  &
                             reftime, month, member, expid, leveltypes, &
                             optionals, unsp_type, ierr                 )
  !----------------------------------------------------------------------------
  ! read atmospheric state variables from GRIB file
  !
  ! atm:       Atmospheric state variable to set. Must be an array in case
  !            of an ensemble. The variable must be initialised already
  !            (call construct).
  ! file:      Name of GRIB file to read. In case of an ensemble the name
  !            is extended by '.eee' (3 digits ensemble number). In a
  !            parallel environment Ensemble members are read in parallel
  !            from different processors.
  ! invt:      Inventory of the GRIB file. May be provided so that it is
  !            not set up twice unnecessarily.
  ! fields:    List of fields to read, separated by blanks.
  ! runtype:   Run-type to select ('forecast','analysis',...)
  ! time:      Verification-time to select.
  ! reftime:   Reference time to select.
  ! member:    Ensemble member to select.
  ! optionals: List of optional fields to read (subset of 'fields').
  !            The routine does not abort if any of these fields
  !            is missing in the file.
  ! unsp_type: optional List of fields with unspecified runtype
  !            (subset of 'fields'). These fields are read even if their
  !            run-type does not correspond with the parameter 'runtype'.
  ! ierr:      error return code
  !----------------------------------------------------------------------------
  type (t_atm)     ,intent(inout)        :: state         ! atm. state variable
  character(len=*) ,intent(in)           :: file          ! grib file name
  type(t_inventory),pointer    ,optional :: invt(:)       ! inventory
  character(len=*) ,intent(in) ,optional :: fields        ! fields to read
  character(len=*) ,intent(in) ,optional :: runtype       ! forecast, analysis
  type(t_time)     ,intent(in) ,optional :: time          ! verification time
  type(t_time)     ,intent(in) ,optional :: reftime       ! reference time
  integer          ,intent(in) ,optional :: month         ! month of verif. time
  integer          ,intent(in) ,optional :: member        ! ensemble member
  integer          ,intent(in) ,optional :: expid         ! experiment id
  integer          ,intent(in) ,optional :: leveltypes(:) ! level types to read
  character(len=*) ,intent(in) ,optional :: optionals     ! optional fields
  character(len=*) ,intent(in) ,optional :: unsp_type     ! unspecified runtype
  integer          ,intent(out),optional :: ierr          ! error return code

    !---------------------------------------------------
    ! additional local parameters passed to read_atm_gen
    !---------------------------------------------------
    type (t_atm) :: atm(1)  ! local copy (vector of size 1)
    integer      :: ie      ! ensemble member index to read on this PE
    integer      :: pio(1)  !
    !---------------------------------
    ! set p_readgrib if not set so far
    !---------------------------------
    if (p_readgrib < 0) then
      p_readgrib = 0
      if (dace% lpio) &
        write(6,'(/a,i4/)') 'read_atm_state: p_readgrib =',p_readgrib
    endif
    !-------------------------------------------------------
    ! set additional parameters to be passed to read_atm_gen
    !-------------------------------------------------------
    atm(1) = state
    pio(1) = p_readgrib
    ie     = 1;   if (.not.dace% lpio) ie = -1
    !------------------------------------------------------
    ! call generic routine to read single atmospheric state
    !------------------------------------------------------
    call read_atm_gen (atm, ie, pio, file, invt, fields, runtype, time,      &
                       reftime, month, member, expid, leveltypes, optionals, &
                       unsp_type, ierr)
    state = atm (1)

  end subroutine read_atm_state
!------------------------------------------------------------------------------
  subroutine read_atm_ens (atm, file, invt, fields, runtype, time, reftime,   &
                           members, expid, leveltypes, optionals, unsp_type,  &
                           ierr, ctr                                          )
  !---------------------------------------------------
  ! read atmospheric ensemble variables from GRIB file
  !---------------------------------------------------
  type (t_atm)     ,intent(inout)        :: atm(:)        ! atm. state variable
  character(len=*) ,intent(in)           :: file          ! grib file name
  type(t_inventory),pointer    ,optional :: invt(:)       ! inventory
  character(len=*) ,intent(in) ,optional :: fields        ! fields to read
  character(len=*) ,intent(in) ,optional :: runtype       ! forecast, analysis
  type (t_time)    ,intent(in) ,optional :: time          ! verification time
  type (t_time)    ,intent(in) ,optional :: reftime       ! reference time
  integer          ,intent(in) ,optional :: members(:)    ! ensemble members to read
  integer          ,intent(in) ,optional :: expid         ! experiment id
  integer          ,intent(in) ,optional :: leveltypes(:) ! level types to read
  character(len=*) ,intent(in) ,optional :: optionals     ! optional fields
  character(len=*) ,intent(in) ,optional :: unsp_type     ! unspecified runtype
  integer          ,intent(out),optional :: ierr          ! error return code
  type(t_ctr_inv)  ,intent(in) ,optional :: ctr(:)        ! inventory container

    type(t_inventory) ,pointer :: inv(:)
    integer                    :: i         ! start index of chunk to read at a time
    integer                    :: n         ! size        of chunk to read at a time
    integer                    :: ie        ! ensemble member index in chunk
    integer                    :: ke        ! ensemble member index
    integer                    :: nchunk    ! number of chunks to read
    integer                    :: mcsize    ! maximum chunk size
    integer                    :: nens      ! ensemble size
    integer, allocatable       :: pio(:)    ! processor elements to use for reading
    integer                    :: lerr      ! error return variable

    if (present (invt)) then
      inv => invt
    else
      nullify (inv)
    endif
    nens = size (atm)
    allocate (pio(nens))
!   call set_iope_ens (pio, nchunk, mcsize, stride=+1)  ! Tight I/O pattern
    call set_iope_ens (pio, nchunk, mcsize, stride=-1)  ! Optimized I/O stride
    if (present(ctr)) then
      if (.not.all(ctr(:)% global)) call check_iope_ens( ctr(:)% pe, pio, 'read_atm_ens')
    end if
    !-----------------------------------------------------------
    ! call read_atm_gen to read 1 member per chosen PE at a time
    !-----------------------------------------------------------
    lerr = IE_OK
    do i = 1, nens, mcsize              ! loop over bunches to read
      n  = min (nens-i+1, mcsize)       ! number of members to read
      ie = -1
      ke = -1
      if (any (pio(i:i+n-1) == dace% pe)) then
        ie = minloc (abs (pio(i:i+n-1) - dace% pe),1) ! index within chunk
        ke = i - 1 + ie                               ! read member on this PE
        if (present(members)) ke = members(ke)
        if (present(ctr)) inv => ctr(ke)% invt
        !print *, "### read_atm_ens: p_pe,ie,ke =", dace% pe, ie, ke
      endif
      if (dace% lpio) then
        write(6,*)
        write(6,*) ' reading ensemble members',i,' to',i+n-1
        write(6,*)
      endif
      if (present(ierr)) ierr = IE_OK
      call read_atm_gen (atm(i:i+n-1), ie, pio(i:i+n-1),             &
                         path_file ('',file, iens=ke), invt=inv,     &
                         fields=fields, runtype=runtype, time=time,  &
!                        reftime=reftime, member=ke, expid=expid,    &! check encoding in GRIB-file
                         reftime=reftime,            expid=expid,    &! do not check
                         leveltypes=leveltypes,                      &
                         optionals=optionals, unsp_type=unsp_type,   &
                         ierr=ierr                                   )
      if (present(ierr)) lerr = min (lerr, ierr)
    end do
    if (present(ierr)) ierr = lerr
    if (present(invt) .or. present(ctr)) then
      nullify (inv)
    else
      if (associated(inv)) deallocate (inv)
    end if

  end subroutine read_atm_ens
!------------------------------------------------------------------------------
  subroutine set_iope_ens (pio, nchunk, mcsize, stride, offset, roundrobin)
    integer, intent(out)           :: pio(:)     ! I/O pe for members
    integer, intent(out)           :: nchunk     ! Number of chunks
    integer, intent(out) ,optional :: mcsize     ! Max. chunk size
    integer, intent(in)  ,optional :: stride     ! PE stride
    integer, intent(in)  ,optional :: offset     ! Initial PE offset
    logical, intent(in)  ,optional :: roundrobin ! Round-robin?
    !-----------------------------------------------
    ! Set I/O processor pattern for ensemble members
    !-----------------------------------------------
    integer :: nens     ! Ensemble size
    integer :: off      ! pe offset
    integer :: inc      ! pe increment
    integer :: ie       ! member index
    integer :: pe       ! processor index
    logical :: rr       ! round-robin with interleave?
    integer :: chsize   ! chunk size estimate

    nens = size (pio)
    inc  = -1; if (present (stride)) inc = stride
    if (inc < 0) inc = max (1, dace% npe / nens)

    off = dace% pio; if (present (offset)) off = offset
    off = max (0, min (off, dace% npe-1))

    rr = .false.
    if (present (roundrobin)) rr = roundrobin
    if (inc == 0)             rr = .false.

    chsize = 1; if (inc > 0) chsize = min (dace% npe / inc, nens)

    pe     = huge (0)
    nchunk = 0
    do ie = 1, nens
       if (pe >= dace% npe .or. inc == 0) then
          pe     = off
          nchunk = nchunk + 1
          if (rr) then
             off = off + 1
             if (off >= inc) off = 0
          end if
       end if
       pio(ie) = pe
       pe      = pe + inc
    end do
!   if (present (mcsize)) mcsize = (nens + nchunk-1) / nchunk ! Avg. chunksize
    if (present (mcsize)) mcsize = chsize                     ! Max. chunksize
!   print *
!   print *, "set_iope_ens: nchunk=", nchunk
!   print *, "set_iope_ens: mcsize=", chsize
!   print *, "set_iope_ens:    pio=", pio
  end subroutine set_iope_ens
!------------------------------------------------------------------------------
  subroutine check_iope_ens (pe, pio, caller)
  !-------------------------------
  ! Compare I/O processor patterns
  !-------------------------------
  integer                      ,intent(in)     :: pe (:)   ! I/O pe pattern
  integer                      ,intent(in)     :: pio(:)   ! reading pe's
  character (len=*) ,optional  ,intent(in)     :: caller   ! name of calling function

    character (len=512)         :: error        ! error message
    logical                     :: lerr         ! error?
    integer                     :: nens         ! ensemble size
    integer                     :: i            ! loop index

    nens = size (pe)
    lerr = .false.
    if (nens /= size(pio)) then
      error = 'size (pe) /= size(pio)'
      lerr = .true.
    else
      error = 'member not found on correct pe'
      do i = 1, nens
        if ( pe(i) /= pio(i) ) then
          error = trim(error)//'; mem='//char3(i)//', pe='//char3(pe(i))//', pio='//char3(pio(i))
          lerr = .true.
        end if
      end do
    end if
    if (lerr) then
      if (present(caller)) error='called from '//caller//': '//trim(error)
      call finish ('check_iope_ens',trim(error))
    end if

  end subroutine check_iope_ens
!------------------------------------------------------------------------------
  subroutine read_atm_gen (state, ie, pio, file, invt, fields, runtype, time, &
                           reftime, month, member, expid, leveltypes,         &
                           optionals, unsp_type, ierr                         )
  type (t_atm)     ,intent(inout)        :: state(:)      ! atm. state variable
  integer          ,intent(in)           :: ie            ! member index to read
  integer          ,intent(in)           :: pio(:)        ! PEs to use for input
  character(len=*) ,intent(in)           :: file          ! grib file name
  type(t_inventory),pointer    ,optional :: invt(:)       ! inventory
  character(len=*) ,intent(in) ,optional :: fields        ! fields to read
  character(len=*) ,intent(in) ,optional :: runtype       ! forecast, analysis
  type(t_time)     ,intent(in) ,optional :: time          ! verification time
  type(t_time)     ,intent(in) ,optional :: reftime       ! reference time
  integer          ,intent(in) ,optional :: month         ! month of verif. time
  integer          ,intent(in) ,optional :: member        ! ensemble member
  integer          ,intent(in) ,optional :: expid         ! experiment id
  integer          ,intent(in) ,optional :: leveltypes(:) ! level types to read
  character(len=*) ,intent(in) ,optional :: optionals     ! optional fields
  character(len=*) ,intent(in) ,optional :: unsp_type     ! unspecified runtype
  integer          ,intent(out),optional :: ierr          ! error return code
  !-----------------------
  ! read atmospheric state
  !-----------------------

    !----------------
    ! local variables
    !----------------
    type (t_inventory) ,pointer :: inv (:)
    type (t_grid)      ,pointer :: grd       ! grid information
    type (t_grib1)              :: grib      ! GRIB record
    integer                     :: ier       ! error return code
    integer                     :: ier_sav   ! error return code
    integer                     :: pi        ! PE to read from
    character(len=8)            :: runtypg   ! runtype to use in general
    character(len=8)            :: runtyp    ! runtype for current variable
    character(len=16)           :: name
    character(len=16)           :: gribname
    integer                     :: lt, l, j
    character(len=16)           :: flds (64)
    character(len=16)           :: urun (64)
    character(len=16)           :: optn (64)
    logical                     :: first
    integer                     :: i
    integer                     :: id        ! record index read
    integer                     :: ierror (size (state))
    !--------------------
    ! set local variables
    !--------------------
    first = .true.
    ier   = 0
    if (present (ierr  )) ierr = IE_OK
    if (present (fields)) then
      call split (flds, fields, lt)
      if (lt<0) call finish('read_atm_state','increase array "flds"')
    else
      call split (flds, 'ps psr t u v q', lt)
    endif
    j = 0
    optn = ' '; if (present(optionals)) call split (optn, optionals, j)
                if (j<0) call finish('read_atm_state','increase "optn"')
    urun = ' '; if (present(unsp_type)) call split (urun, unsp_type, j)
                if (j<0) call finish('read_atm_state','increase "urun"')
    grd => state(1)% grid
    pi = -1; if (ie > 0) pi = pio(ie)
    !--------------
    ! get inventory
    !--------------
    if (present (invt)) then
      inv => invt
    else
      nullify (inv)
    endif
    if (.not.associated(inv) .and. dace% pe==pi) then
!     call get_inventory (inv, file, pio=pi, comm=ci)
      call get_inventory (inv, file, pio=dace% pe, comm=self% comm)
    endif
    !----------------------------------------------------
    ! check if inventory is derived and member is present
    !----------------------------------------------------
    if (dace% pe==pi .and. .not. associated (inv))  ier = IE_FILE
    if (error (ier, 'no inventory: '//trim (file))) return
    if (dace% pe==pi) then
      if (size (inv) == 0)                          ier = IE_FILE
    end if
    if (error (ier, 'no inventory: '//trim (file))) return
    if (present (member)) then
      if (dace% pe==pi) then
        if (.not. any (inv(:)% en% no == member))   ier = IE_MEMBER
      end if
      if (error (ier, 'member not in inventory: '//trim (string(member))//&
                      '  file: '//trim (file)   )) return
    endif
    !---------------
    ! open GRIB file
    !---------------
    if (dace% pe == pi) call open_gribfile (grib, file, 'r', kret=ier)
    if (ier /= 0) ier = IE_FILE
    if (error (ier, 'cannot open: '//trim (file))) return
    !---------------
    ! select runtype
    !---------------
    if (dace% pe == pio(1)) runtypg = inv(1)%pa% runtype
    call p_bcast (runtypg, pio(1))
    if (present(runtype)) then
      select case (runtype)
      case ('first')        ! take from 1st record in file
      case default
        runtypg = runtype   ! explicitly specified
      end select
    endif
    !----------------------------
    ! loop over variables to read
    !----------------------------
    do l = 1, lt
      name = flds (l)
      if (name=='') cycle
      runtyp = runtypg
      if (runtyp    == 'any') runtyp = ''
      if (any (urun == name)) runtyp = ''
      !--------------
      ! read variable
      !--------------
      if (dace% lpio) write(6,*) 'read_atm_state: reading ',trim (name)
      gribname = name
      ier_sav = 0
      select case (name)
      case default
        call allocate (state(:), name, ierr=ierror)
        if (any (ierror/=0)) then
          if (all (optn /= name)) then
             call finish ('read_atm_state',&
                          'allocation failed for mandatory: '//trim(name))
          else
             call message('read_atm_state',&
                          'allocation failed for optional: '//trim(name))
             cycle
          end if
        end if
        call read_sm (name, name, id, ier, member=member, expid=expid)
        !----------------------------------
        ! second try for specific variables
        !----------------------------------
        if (ier/=0) then
          select case (name)
          case ('clcl')
            call read_sm (name, name, id, ier, member=member, expid=expid, level=800)
          case ('clcm')
            call read_sm (name, name, id, ier, member=member, expid=expid, level=400)
          case ('clch')
            call read_sm (name, name, id, ier, member=member, expid=expid, level=000)
          case ('vis')
            call read_sm (name, name, id, ier, member=member, expid=expid, level=000)
          case ('pp')
            if (dace% lpio) write(6,*) 'read_atm_state:  read P instead of PP'
            call read_sm ('p', name, id, ier, member=member, expid=expid)
            if(ier==0) gribname = 'p'
            if(ier==0) then
              do i=1, size(state)
                if (.not.associated(state(i)% grid% p0)) &
                  call finish ('read_atm_state','pp = p - p0: p0 not associated')
                state(i)% pp = state(i)% pp - state(i)% grid% p0
                where (state(i)%m%i%name == 'pp')
                  state(i)%m%i%code  = 255
                  state(i)%m%i%table = 255
                end where
              end do
            endif
          case ('ph')
            if (dace% lpio) write(6,*) 'read_atm_state:  read P instead of PH'
            call read_sm ('p', name, id, ier, member=member, expid=expid)
            if(ier==0) gribname = 'p'
          case ('pf')
            if (dace% lpio) write(6,*) 'read_atm_state:  read P instead of PF'
            call read_sm ('p', name, id, ier, member=member, expid=expid)
            if(ier==0) gribname = 'p'
          case ('o3')
            if (dace% lpio) write(6,*) 'read_atm_state:  read GO3 instead of O3'
            call read_sm ('go3', name, id, ier, member=member, expid=expid)
            if(ier==0) gribname = 'go3'
          end select
        endif
      case ('ps','psr')
        !-----------------------------------------------------------
        ! special handling for surface pressure / reference pressure
        !-----------------------------------------------------------
        call allocate (state(:), name)
        call read_sm ('ps', name, id, ier, member=member, expid=expid)
        if(ier==0) gribname = 'ps'
        if (ier /= 0) then
          !--------------------------------------------------------
          ! Next attempt: ECMWF's sp is model level 1, leveltyp=109
          !--------------------------------------------------------
          if (ier == IE_AMBIGUOUS) then
            print *, "read_atm_state:  Warning: ps is ambiguous!"
          end if
          ier_sav = min (ier_sav, ier)
          if (dace% lpio) print *, "read_atm_state:  now trying read_m(ps)"
          call read_sm ('ps', name, id, ier, member=member, expid=expid, &
                        force_m=.true.                                    )
          if(ier==0) gribname = 'ps'
        endif
        !----------------------------------------
        ! alternatively read log surface pressure
        !----------------------------------------
        if (ier/=0) then
          if (ier == IE_AMBIGUOUS) then
            print *, "read_atm_state:  Warning: ps is ambiguous!"
          end if
          ier_sav = min (ier_sav, ier)
          if (dace% lpio) &
            print *, "read_atm_state:  read_s(ps) failed, trying lnsp instead"
          call read_sm ('lnsp', name, id, ier, member=member, expid=expid)
          if (ier==0) gribname = 'lnsp'
          if (ier /= 0) then
            if (ier == IE_AMBIGUOUS) then
              print *, "read_atm_state:  Warning: lnsp is ambiguous!"
            end if
            ier_sav = min (ier_sav, ier)
            !------------------------------------------------------------
            ! last attempt: ECMWF's ln(sp) is model level 1, leveltyp=109
            !------------------------------------------------------------
            if (dace% lpio) print *,                                          &
              "read_atm_state:  read_s(lnsp) failed, now trying read_m(ln(sp))"
            call read_sm ('lnsp', name, id,ier, member=member, expid=expid,  &
                          force_m=.true.                                      )
            if (ier==0) gribname = 'lnsp'
          end if
          !-------------------------------
          ! if ln(p) was read convert to p
          !-------------------------------
          if (ier == 0) then
            do i=1,size(state)
              if (name == 'ps') then
                state(i)% ps  = exp (state(i)% ps)
              else
                state(i)% psr = exp (state(i)% psr)
              endif
            end do
          endif
        endif
      end select
      !---------------------------------------------
      ! report missing, continue for optional fields
      !---------------------------------------------
      if (ier/=0) then
        if (ier == IE_AMBIGUOUS) then
           print *, "read_atm_state:  Warning: field is ambiguous!"
        end if
        ier_sav = min (ier_sav, ier)
        if (all (optn /= name)) then
          call error_printout
          if (ier_sav == IE_AMBIGUOUS) then
            if (error (ier_sav,                             &
             "field '"//trim(name)//"' is ambiguous")) return
           else if (ier_sav == IE_LEVELS) then
            if (error (ier_sav,                             &
             "missing levels for '"//trim(name)//"'")) return
          else
            if (error (IE_MISSING,                          &
              "cannot find field '"//trim(name)//"'")) return
          end if
        else
          call deallocate (state(:), name)
        endif
      endif
      !--------------------------------------
      ! set meta data from first valid field:
      !   time, reftime
      !   ensemble member, size
      !   experiment id, runclass
      !--------------------------------------
      if (ier==0 .and. first) then
        first = .false.
        if (dace% pe == pi) then
          !-------------------------------------
          ! use current inventory index if valid
          !-------------------------------------
          i = id
          if (id==0) then
            do i = 1, size (inv)
               if (inv(i)% pa% iname == gribname) exit
            end do
          endif
          if (i > size (inv)) i = 1 ! Fallback if gribname /= statename
          !--------------
          ! set meta data
          !--------------
          state(ie)% time          = inv(i)% ti% ver_time
          state(ie)% ref_time      = inv(i)% ti% ref_time
          if(runtyp == runtypg) then
            state(ie)% runtype     = inv(i)% pa% runtype
            runtypg                = inv(i)% pa% runtype
          endif
          state(ie)% runclass      = inv(i)% pa% runclass
          state(ie)% expid         = inv(i)% pa% expid
          state(ie)% member        = inv(i)% en% no
          state(ie)% members       = inv(i)% en% size
          state(ie)% ensemble_id   = inv(i)% en% id
        endif
        if (present(member))  state(:)% member   = member
        if (present(time))    state(:)% time     = time
        if (present(reftime)) state(:)% ref_time = reftime
        if (p_parallel) then
          do i=1,size(state)
           call p_bcast (state(i)% time,        pio(i))
           call p_bcast (state(i)% ref_time,    pio(i))
           call p_bcast (state(i)% member,      pio(i))
           call p_bcast (state(i)% members,     pio(i))
           call p_bcast (state(i)% runclass,    pio(i))
           call p_bcast (state(i)% expid,       pio(i))
           call p_bcast (state(i)% runtype,     pio(i))
           call p_bcast (state(i)% ensemble_id, pio(i))
           call p_bcast (          runtypg,     pio(i))
          end do
        endif
      endif
    end do
    !--------
    ! cleanup
    !--------
    call set_pointers (state)
    if (dace% pe == pi) call close_gribfile (grib)

    if (.not. present(invt) .and. associated(inv)) deallocate (inv)
  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine error_printout
    !---------------------------------------------------------
    ! print out search criteria if field was not found in file
    !---------------------------------------------------------
      if (.not.dace% lpio) return
      select case (ier_sav)
      case (IE_AMBIGUOUS)
       write (6,*)  'read_atm_state:  ambiguous field    ',trim(name)
      case (IE_LEVELS)
       write (6,*)  'read_atm_state:  missing levels for ',trim(name)
      case default
       write (6,*)  'read_atm_state:  cannot find field  ',trim(name)
      end select
      write (6,*)                          '  file       = ',trim(file)
      write (6,*)                          '  ierr       = ',ier
      if (present(fields    )) write (6,*) '  fields     = ',trim(fields)
      if (present(runtype   )) write (6,*) '  runtype    = ',runtype
      write (6,*)                          '  runtyp     =>',runtyp
      if (present(time      )) write (6,*) '  valid time = ',cyyyymmddhhmm(time)
      if (present(reftime   )) write (6,*) '  reftime    = ',cyyyymmddhhmm(reftime)
      if (present(member    )) write (6,*) '  member     =', member
      if (present(expid     )) write (6,*) '  expid      = ',expid
      if (present(leveltypes)) write (6,*) '  leveltypes = ',leveltypes
      if (present(optionals )) write (6,*) '  optionals  = ',trim(optionals)
      if (present(unsp_type )) write (6,*) '  unsp_type  = ',trim(unsp_type)

    end subroutine error_printout
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function error (ie, text)
    integer          ,intent(in) :: ie
    character(len=*) ,intent(in) :: text
    logical                      :: error

      if (present (ierr)) then
        ierr  = p_min (ie)
        error = ierr < IE_OK
      else
        if (ie < 0) then
          if (any (dace% pe == pio)) then
            write (6, *) 'read_atm_gen:',trim(text)//' in file '//trim(file)
            call finish ('read_atm_gen' ,trim(text)//' in file '//trim(file))
          else
            ! Non-I/O pe:
            call finish ('read_atm_gen' ,trim(text))
          end if
        endif
        error = .false.
      endif

    end function error
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine read_sm (gribname, statename, id, ierr, &
                        member, expid, force_m, level)
    !------------------------------------
    ! read single or multi-level variable
    !------------------------------------
    character(len=*) ,intent(in)            :: gribname
    character(len=*) ,intent(in)            :: statename
    integer          ,intent(out)           :: id
    integer          ,intent(out)           :: ierr
    integer          ,intent(in)  ,optional :: member
    integer          ,intent(in)  ,optional :: expid    ! experiment id
    logical          ,intent(in)  ,optional :: force_m  ! Enforce read_m
    integer          ,intent(in)  ,optional :: level    ! specific level

      integer :: i    ! variable loop index
      integer :: is   ! state loop index
      logical :: f_m  ! local copy of force_m

      f_m = .false.; if (present (force_m)) f_m = force_m
      ierr = 0
      id   = 0
      !-----------------------------------
      ! find variable in atmospheric state
      !-----------------------------------
      do i=1,size (state(1)%m)
        if   (state(1)% m(i)%i% name  == statename) then
          !-------------------------------------------
          ! read either single or multi-level variable
          !-------------------------------------------
          if (state(1)% m(i)%i% ub(3) ==            &
              state(1)% m(i)%i% lb(3) .and..not. f_m) then
            call read_s (gribname, statename, state% m(i), id, ierr, &
                         member=member, expid=expid, level=level     )
          else
            call read_m (gribname, statename, state% m(i), id, ierr, &
                         member=member, expid=expid                  )
          end if
          !-------------------------------------------
          ! transform hPa -> Pa for pressure deviation
          !-------------------------------------------
          if (ierr==0 .and. gribname == 'pp') then
            do is = 1, size(state)
              state(is)% m(i)% ptr = state(is)% m(i)% ptr * 100._wp
            end do
          endif
          !-------------------------------------
          ! transform % -> 1 for 2m rel.humidity
          !-------------------------------------
          if (ierr==0 .and. gribname == 'rh2m') then
            do is = 1, size(state)
              state(is)% m(i)% ptr = state(is)% m(i)% ptr * 0.01_wp
            end do
          endif
          return
        end if
      end do
      !----------------------------------------------------
      ! error handling, field not present in state variable
      !----------------------------------------------------
      call finish ('read_sm',&
        "cannot find field '"//trim(statename)//"' in model state.")
    end subroutine read_sm
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine read_s (gribname, statename, m, id, ierr, member, expid, level)
    character(len=*) ,intent(in)            :: gribname
    character(len=*) ,intent(in)            :: statename
    type(t_m)        ,intent(inout)         :: m(:)
    integer          ,intent(out)           :: id
    integer          ,intent(out)           :: ierr
    integer          ,intent(in)  ,optional :: member
    integer          ,intent(in)  ,optional :: expid      ! experiment id
    integer          ,intent(in)  ,optional :: level

      integer           :: code, table
      integer           :: dis, cat, num
      integer           :: i
      real(wp) ,pointer :: xg (:,:,:)
      real(wp) ,pointer :: x4 (:,:,:,:)
      logical           :: decode
      integer           :: pmode

      !---------------------------------
      ! find record in GRIB file, decode
      !---------------------------------
      decode = .false.; if(ie>0) decode = dace% pe == pio(ie)
      pmode  = SCATTER
      code   = 255; table = 255
      dis    = 255; cat   = 255; num = 255
      id     = 0
      ierr   = 0
      nullify (xg, x4)

      if (decode) then
        call decode_single (id, inv, ierr=ierr, code=code, table=table,   &
                            dis=dis, cat=cat, num=num,                    &
                            name=gribname, runtype=runtyp, time=time,     &
                            reftime=reftime, member=member, expid=expid,  &
                            grd=grd, grib=grib, scanmode=WMO8_J_POSITIVE, &
                            level=level, x4=x4)
        if (ierr == 0) xg => x4(:,:,1,:)
      endif
      !---------------
      ! gather / bcast
      !---------------
      ierr = p_min (ierr)
      if (ierr == 0) then
        select case (pmode)
        case (SCATTER)
!         do i = 1, size(m)
!           call scatter_level (xg, m(i)% ptr(:,:,1,:), grd%dc, pio(i))
!         end do
          if (size(m) == 1) then
!             call scatter_level (xg, m(1)% ptr(:,:,1,:), grd%dc, pio(1))
              call scatter_multi (x4, m(1)% ptr,          grd%dc, pio(1))
          else
            call alltoallv_multi (x4, m, grd%dc, pio(:))
          endif
        case (BCAST)
          do i = 1, size(m)
            if (i==ie)    m(i)% ptr(:,:,1,:) = xg
            call p_bcast (m(i)% ptr(:,:,:,:), pio(i))
          end do
        case (0)
          m(1)% ptr(:,:,1,:) = xg
        end select
      end if

      !---------------
      ! error handling
      !---------------
      if      (associated(x4))  then
         deallocate (x4)
      else if (associated(xg))  then
         deallocate (xg)
      endif
      if (ierr /= 0) return

      !--------------------------------
      ! keep GRIB code and table number
      !--------------------------------
      m% i% code  = code
      m% i% table = table
      m% i% dis   = dis
      m% i% cat   = cat
      m% i% num   = num

    end subroutine read_s
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine read_m (gribname, statename, m, id, ierr, member, expid)
    character(len=*) ,intent(in)            :: gribname
    character(len=*) ,intent(in)            :: statename
    type(t_m)        ,intent(inout)         :: m(:)
    integer          ,intent(out)           :: id
    integer          ,intent(out)           :: ierr
    integer          ,intent(in)  ,optional :: member
    integer          ,intent(in)  ,optional :: expid      ! experiment id

      integer           :: code, table
      integer           :: dis, cat, num
      integer           :: i
      real(wp) ,pointer :: xg (:,:,:,:)
      logical           :: decode
      integer           :: pmode
      integer           :: k1, k2        ! lowest, highest level number
      integer           :: kk, kt, dk, n ! Auxiliaries for level blocking
      integer(i8)       :: lev_size      ! Bytes / level

      !---------------------------------
      ! find record in GRIB file, decode
      !---------------------------------
      decode = .false.; if(ie>0) decode = dace% pe == pio(ie)
      pmode  = SCATTER
      code   = 255; table = 255
      dis    = 255; cat   = 255; num = 255
      ierr   = 0

      k1 = lbound (m(1)% ptr,3)
      k2 = ubound (m(1)% ptr,3)
      dk = k2 - k1 + 1
      if (pmode == SCATTER .and. io_max_scatter >= 0) then
         lev_size = (grd% ubg(1)-grd% lbg(1)+1) * (grd% ubg(2)-grd% lbg(2)+1)
         lev_size = lev_size                    * (grd% ubg(4)-grd% lbg(4)+1)
         lev_size = lev_size * (storage_size (0._wp) / 8)
         dk       = min (max (int (io_max_scatter / lev_size), 1), dk)
      end if
      !----------------------------
      ! Level blocking (outer loop)
      !----------------------------
      do kk = k1, k2, dk
         kt = min (kk+dk-1, k2)
         n  = kt - kk + 1

      nullify (xg)
      if (decode)                                                         &
        call decode_multi (xg, inv, grd, grib, n, decode,                 &
                           id,                                            &
                           ierr=ierr,                                     &
                           scanmode=WMO8_J_POSITIVE,                      &
                           code=code, table=table,                        &
                           dis=dis, cat=cat, num=num,                     &
                           name=gribname,                                 &
                           runtype=runtyp,                                &
                           time=time, reftime=reftime, month=month,       &
                           kf=kk,                                         &
                           member=member,                                 &
                           expid=expid,                                   &
                           leveltypes=leveltypes                          )
      !---------------
      ! gather / bcast
      !---------------
      ierr = p_min (ierr)
      if (ierr == 0) then
        if (dace% lpio .and. kk == k1 .and. size(m) > 1) &
             write(6,*) 'scatter fields'
FTRACE_BEGIN("read_m:scatter_level")
        select case (pmode)
        case (SCATTER)
!         do i = 1, size(m)
!           call scatter_multi (xg, m(i)% ptr, grd%dc, pio(i))
!         end do
          if (size(m) == 1) then
            call scatter_multi (xg, m(1)% ptr(:,:,kk:kt,:), grd%dc, pio(1))
          else
            call alltoallv_multi (xg, m, grd%dc, pio(:), lb3=kk, ub3=kt)
          endif
        case (BCAST)
          do i = 1, size(m)
            if (i==ie)    m(i)% ptr = xg
            call p_bcast (m(i)% ptr, pio(i))
          end do
        case (0)
            m(1)% ptr = xg
        case default
        end select
FTRACE_END  ("read_m:scatter_level")
      endif

      !---------------
      ! error handling
      !---------------
      if (associated (xg)) deallocate (xg)
      if (ierr /= 0) return

      end do    ! kk loop (level blocking)

      !--------------------------------
      ! keep GRIB code and table number
      !--------------------------------
      m% i% code  = code
      m% i% table = table
      m% i% dis   = dis
      m% i% cat   = cat
      m% i% num   = num

    end subroutine read_m
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine read_atm_gen
!==============================================================================
  subroutine read_atm_grid (grid, file, invt, invar, nproc1, nproc2, comm, &
                            geosp, lsm, hhl, optionals,                    &
                            scanmode, member, gridfile,                    &
                            rf, rni, nzr, g_coarse                         )
  type(t_grid)     ,intent(out)          :: grid     ! grid variable
  character(len=*) ,intent(in)           :: file     ! grib file name
  type(t_inventory),pointer    ,optional :: invt(:)  ! inventory
  character(len=*) ,intent(in) ,optional :: invar    ! invariant fields file
  integer          ,intent(in) ,optional :: nproc1   ! number of PEs in x
  integer          ,intent(in) ,optional :: nproc2   ! number of PEs in y
  integer          ,intent(in) ,optional :: comm     ! MPI communicator
  logical          ,intent(in) ,optional :: geosp    ! Add orography
  logical          ,intent(in) ,optional :: lsm      ! Add land-sea-mask
  logical          ,intent(in) ,optional :: hhl      ! Add half levels (ICON)
  character(len=*) ,intent(in) ,optional :: optionals! optional fields
  integer          ,intent(in) ,optional :: scanmode ! Forced scanmode
  integer          ,intent(in) ,optional :: member   ! Ensemble member
  character(len=*) ,intent(in) ,optional :: gridfile ! grid metadata file
  !----------------------------------------------------
  ! optional arguments for coarse LETKF grid generation
  !----------------------------------------------------
  integer          ,intent(in) ,optional :: rf       ! coarsening factor
  integer          ,intent(in) ,optional :: rni      ! GMI resolution
  integer          ,intent(in) ,optional :: nzr      ! number of levels
  type(t_grid)     ,pointer    ,optional :: g_coarse ! coarse grid

  !----------------------------------------------------------------------
  ! Derive the model grid from a GRIB file.
  ! The variable 'grid' must not be initialized (don't call constuct before).
  !
  ! The GDS of a suitable record in the file 'file' is evaluated.
  ! The variable 'grid' is set accordingly.
  ! Additional fields (orography, land sea mask) are read from 'file' and
  ! stored within 'grid'. If 'invar' is given, the invariant fields are
  ! read from a seperate file with this name. If 'geosp' or 'lsm' are passed
  ! with value .false. the routine will not abort if the respective fields are
  ! missing.
  !
  ! If the inventory 'invt' of the GRIB file is given, it is not
  ! derived from the file. If 'invt' is passed but not associated, the
  ! inventory is provided on output.
  !
  ! If the model fields shall be distributed in a parallel environment
  ! 'nproc1', 'nproc2' must be provided so that nproc1*nproc2 equals the
  ! total number of processors in the group. Otherwise the complete model
  ! fields will be allocated on each processor. If 'comm' is not present
  ! MPI_COMM_SELF will be assumed for the processor group.
  !----------------------------------------------------------------------
    !----------------
    ! local variables
    !----------------
    integer                     :: i, j, k
    integer                     :: irec      ! record in inventory
    integer                     :: gridtype  ! grid type
    type (t_inventory) ,pointer :: inv (:)   ! inventory
    type (t_inventory) ,pointer :: iinv (:)  ! inventory for invariant fields
    type (t_grib1)              :: grib      ! GRIB record
    real(wp)           ,pointer :: a(:),b(:) ! vertical coordinate parameters
    real(wp)           ,pointer :: d(:)      ! depth below surface (cm)
    integer                     :: nv        ! no. of vert. coord. parameters
    integer                     :: ke        ! number of vertical levels + 1
    integer                     :: nx        ! number of grid points in x dir.
    integer                     :: ny        ! number of grid points in y dir.
    integer                     :: ni        ! number of sections /diamond side
    integer                     :: nd        ! number of diamonds
    integer                     :: ns        ! number of soil layers
    integer                     :: ierr      ! error return code
    integer                     :: lat1, lon1, levtyp
    real(wp)                    :: dx, dy    ! grid spacing
    integer                     :: latr      ! latitude  of S.pole of rotation
    integer                     :: lonr      ! longitude of S.pole of rotation
    logical                     :: lgeosp    ! copy of optional parameter geosp
    logical                     :: llsm      ! copy of optional parameter lsm
    logical                     :: linvar    ! copy of optional parameter invar
    logical                     :: lhhl      ! copy of optional parameter hhl
    character(len=128)          :: lopt      ! copy of optional parameter optionals
    integer                     :: latx      ! temporary
    character                   :: arakawa   ! 'A' or 'C'
    integer                     :: scan      ! Scanmode
    real(wp)           ,pointer :: x3(:,:,:) ! 3d auxiliary field
!   real(wp), allocatable       :: lon(:,:,:)! Gridpoint longitudes
!   real(wp), allocatable       :: lat(:,:,:)! Gridpoint latitudes
    type(t_grid)                :: gtmp      ! temporary grid (for ICON hor.)
    integer                     :: nproc
    character(len=1)            :: uuid_v(16)! UUID of vertical grid
    logical                     :: chk_icon  ! check order of neighbours
    !------------------------------------------
    ! variables to interface with get_vertcoord
    !------------------------------------------
    type(t_vcoord)    :: vc              ! COSMO vertical coordinate parameters
    type(t_refatm)    :: ra              ! COSMO reference atmosphere parameters
    integer           :: ivctype   = -1  ! COSMO vertical coordinate  type
    integer           :: ke1       =  0  ! number of half levels

    !---------------------------------------------
    ! process optional parameters
    ! decide which PE is used for reading the file
    !---------------------------------------------
    llsm   = .true. ; if (present (lsm))       llsm   = lsm
    linvar = .false.; if (present (invar))     linvar = invar /= ''
    lhhl   = .true. ; if (present (hhl))       lhhl   = hhl
    lopt   = ''     ; if (present (optionals)) lopt   = optionals
    lgeosp = index (lopt,'geosp') == 0
                      if (present (geosp))     lgeosp = geosp
    if (p_readgrib < 0) p_readgrib = 0
    if (dace% lpio) then
      write(6,'(/a,i0)') 'read_atm_grid: p_readgrib = ',p_readgrib
      write(6,'(a,a /)') 'read_atm_grid: file       = ',trim(file)
    endif
    !--------------
    ! get inventory
    !--------------
    if (present (invt)) then
      inv => invt
    else
      nullify (inv)
    endif
    if (.not. associated (inv))  call get_inventory (inv, file, pio=p_readgrib)
    if (.not.associated (inv))                                    &
         call finish('read_atm_grid','no inventory: '//trim (file))
    if (size(inv) == 0)                                           &
         call finish('read_atm_grid','no inventory: '//trim (file))
    !------------------------------------
    ! Read ICON horizontal grid metadata;
    ! for now only warn if not present.
    !------------------------------------
    gtmp% size = 0
    nullify (gtmp% icongrid)
    if (any (inv(:)% gr% gridtype == DWD6_ICON)) then
       ni = maxval (inv(:)% gr% ni, mask = (inv(:)% gr% gridtype == DWD6_ICON))
       nx = maxval (inv(:)% gr% nxny,                                &
                    mask = (inv(:)% gr% gridtype == DWD6_ICON) .and. &
                           (inv(:)% gr% grid_ref ==     1    )       ) ! Cells
       i  = minloc (inv(:)% gr% ni, mask = (inv(:)% gr% ni == ni), dim=1)
       gtmp% gridtype = DWD6_ICON
       gtmp% ni       = ni
       gtmp% nxny     = nx
       if (present (gridfile)) then
          if (gridfile /= "") then
             if (dace% lpio) &
                  write(*,*) "ICON grid requesting uuid = ", &
                  byte2hex (inv(i)% gr% uuid)
             call read_icon_metadata (gtmp, gridfile, pio=p_readgrib)
          end if
       end if
       if (.not. associated (gtmp% icongrid)) &
            call finish ("read_atm_grid", "no metadata for ICON grid")
       call verify_grid_num (inv, gtmp% icongrid% grid_num)
       !-------------------------------------------------------------
       ! check neighbour relationships in icon grid for correct order
       !-------------------------------------------------------------
       select case (chk_icon_grid)
       case (1)
         chk_icon = .true.
       case (0)
         chk_icon = .false.
       case default
         !-----------------------------------------------------
         ! default: try to estimate broken grids for nests only
         !-----------------------------------------------------
         if (gtmp% icongrid% global) then
           chk_icon = .false.
         else
           chk_icon = .true.
         endif
       end select
       if (chk_icon) then
         if (dace% lpio) then
           write (6,*)
           write (6,*) 'check neighbour relationships in ICON grid (num, global)'&
                       , gtmp% icongrid% grid_num, gtmp% icongrid% global
           write (6,*)
         endif
         call check_neighbours (gtmp% icongrid)
       end if
    end if
    !---------------------------------------
    ! open GRIB file, process on one PE only
    !---------------------------------------
    nullify (a)
    nullify (b)
    nullify (d)
    uuid_v = achar (0)
    levtyp = -1
    ke1    =  0
    if (dace% pe == p_readgrib) then
      call open_gribfile (grib, file, 'r')
      !-----------------------------------------------
      ! scan the file for GDS from hybrid level record
      !-----------------------------------------------
      irec    = 0
      !-------------------------------------------------
      ! 1st try, do not use wind fields (staggered grid)
      !-------------------------------------------------
      do i=1,size(inv)
        select case (inv(i)%lv% leveltype)
        case (WMO3_HYBRID, WMO3_HYBRIDB, WMO3_HHYBRID, WMO3_GENV)
          select case (inv(i)%pa% shortname)
          case('U','V','W','u','v','w')
            cycle
          case default
            irec   = i
            levtyp = inv(i)%lv% leveltype
            exit
          end select
        end select
      end do
      !------------------------------
      ! 2nd try, now take u,v as well
      !------------------------------
      if (irec == 0) then
        do i=1,size(inv)
          select case (inv(i)%lv% leveltype)
          case (WMO3_HYBRID, WMO3_HYBRIDB, WMO3_HHYBRID, WMO3_GENV)
          irec   = i
          levtyp = inv(i)%lv% leveltype
          exit
        end select
        end do
      endif
      !--------------------------------------------------
      ! if hybrid level is found: read record, decode GDS
      !--------------------------------------------------
      if (irec /= 0) then
        call pbseek      (grib, inv(irec)%ds% ptr)
!       call pbgrib      (grib)
!       call gribex      (grib, 'J')
        call read_gribex (grib, 'J')
        nv = grib% isec2 (12)
        !------------------------------------
        ! check for GME/HRM hybrid pressure
        !        or COSMO vertical coordinate
        !------------------------------------

print *,'derive vertical coordinates grom GDS:'
print *,'center      =',grib% isec1% center
print *,'sub_center  =',grib% isec1% sub_center
print *,'process     =',grib% isec1% process
print *,'local_ident =',grib% isec1% local_ident
print *,'repres      =',grib% isec2  (1)
print *,'nv          =',nv
do i=1,nv
   print '(i8,f18.8)',i,grib% rsec2% vcp(i)
end do

!       if ((grib% isec2     (1) == WMO6_ROTLL  .and. &
!            grib% isec2    (12) == 45          .and. &
!            grib% rsec2% vcp(1) >= 100000._dp)  .or. &
!           (grib% isec2     (1) == WMO6_ROTLL  .and. &
!            grib% isec2    (12) == 55          .and. &
!            grib% rsec2% vcp(1) >= 99999._dp))  then

        if (     grib% isec2     (1) == WMO6_ROTLL          .and. &
            (abs(grib% rsec2% vcp(1) -  100000._dp) < 1._dp .or.  &
             abs(grib% rsec2% vcp(3) -  100000._dp) < 1._dp     ) ) then
          !----------------------------------------
          ! assume COSMO vertical coordinate system
          !----------------------------------------
          call get_vertcoord_new (nv, grib% rsec2% vcp, vc, ra)
          ivctype = vc% ivctype
          ke1     = vc% nlevels

print *,'COSMO vertical coordinate parameters:'
print *,'ke1         =',ke1
print *,'ivctype     =',ivctype
print *,'irefatm     =',ra% irefatm
print *,'p0sl        =',ra% p0sl
print *,'t0sl        =',ra% t0sl
print *,'dt0lp       =',ra% dt0lp
print *,'vcflat      =',vc% vcflat
print *,'delta_t     =',ra% delta_t
print *,'h_scal      =',ra% h_scal

        else if (levtyp == WMO3_GENV) then

          if (inv(irec)%lv% levels(2) == inv(irec)%lv% levels(1)+1) then
             ke1 = inv(irec)% nlev + 1  ! Derive from number of full-levels or layers
          else
             ke1 = inv(irec)% nlev      ! Derive from number of half-levels
          end if

          if (inv(irec)%lv% nlev > -1 .and. inv(irec)%lv% nlev /= ke1) then
             print *, "read_atm_grid: Warning: nlev /= ke+1"
             ke1 = inv(irec)%lv% nlev
          end if
          print *,'Generalized vertical coordinate, nlev = nhalf =', ke1

          ivctype = inv(irec)%lv% grid_num
          uuid_v  = inv(irec)%lv% uuid
          select case (ivctype)
          case (1:3) ! COSMO
             print *, "COSMO: numberOfVGridUsed =",  ivctype
             print *, "       uuidOfVGrid       = ", byte2hex (uuid_v)
          case (IVCTYPE_ICON)
             print *, "ICON:  numberOfVGridUsed =",  ivctype
             print *, "       uuidOfVGrid       = ", byte2hex (uuid_v)
          case default
             write(0,*) "ivctype: unsupported numberOfVGridUsed =", ivctype
             call finish ("read_atm_grid","unsupported ivctype")
          end select
          nv = 0
          allocate (a(0), b(0))
          vc% vc_uuid = uuid_v
          vc% ivctype = ivctype

!TODO: check uuidOfHGrid, uuidOfVGrid

        else

          !-------------------------------------------
          ! assume GME/HRM hybrid pressure coordinates
          !-------------------------------------------
          allocate (a (nv/2), b (nv/2))
          a = grib% rsec2% vcp (     1:nv/2)
          b = grib% rsec2% vcp (nv/2+1:nv  )

          select case (levtyp)
          case (WMO3_HYBRID, WMO3_HYBRIDB)
            if (a(1) < a(2)) then
print *,'GME/HRM vertical coordinate'
              ivctype = 0
            else
print *,'ICON vertical coordinate, disguised as hybrid?'
              levtyp  = WMO3_HHYBRID ! or WMO3_GENV
              ivctype = IVCTYPE_ICON
            end if
          case (WMO3_HHYBRID)        ! Height-based hybrid levels
print *,'ICON vertical coordinate?'
            ivctype = IVCTYPE_ICON
          case default
            call finish('read_atm_grid','unknown/unsupported hybrid levels')
          end select
        endif
      else
      !----------------------------------------------
      ! if hybrid level is not found:
      ! get grid type, grid size from isobaric levels
      !----------------------------------------------
        do i=1,size(inv)
          select case (inv(i)%lv% leveltype)
          case (WMO3_ISOBARIC)
            !-----------------------------------------------------------
            ! Skip isobaric layers (for which second surface is present)
            !-----------------------------------------------------------
            if (inv(i)% ds% edition   >  1 .and. &
                inv(i)% lv% levels(2) >= 0       ) cycle
            irec = i
            levtyp = inv(i)%lv% leveltype
            exit
          case default
            levtyp = inv(i)%lv% leveltype
            if (i == size(inv)) print *, "read_atm_grid:&
              & No hybrid/isobaric levels present; using leveltype", levtyp
          end select
        end do
        if (irec /= 0) then
          nv = 2* count (inv(irec)%pa == inv(:)%pa &
                   .and. inv(irec)%ti == inv(:)%ti &
                   .and. inv(irec)%en == inv(:)%en)
          allocate (a (nv/2), b (nv/2))
          b = 0._wp
          k = 0
          do j=irec,size(inv)
            if     (inv(irec)%pa == inv(j)%pa &
              .and. inv(irec)%ti == inv(j)%ti &
              .and. inv(irec)%en == inv(j)%en ) then
              k=k+1
              a(k)=inv(j)%lv% levelvalue * 100._wp
            endif
          end do
          if (a(1)>a(nv/2)) a = a((/(j,j=nv/2,1,-1)/))
          call pbseek      (grib, inv(irec)%ds% ptr)
!         call pbgrib      (grib)
!         call gribex      (grib, 'J')
          call read_gribex (grib, 'J')
        endif
      endif
      !----------------------------
      ! level type is not supported:
      ! report surface type
      !----------------------------
      if (irec == 0) then
!       call pbgrib      (grib)
!       call gribex      (grib, 'J')
        call read_gribex (grib, 'J')
        nv     = 0
        allocate (a (nv/2), b (nv/2))
      endif
      !--------------------------
      ! extract subsurface levels
      !--------------------------
      ns = 0
      allocate (d(ns))
      do i=1,size(inv)
        select case (inv(i)%lv% leveltype)
        case (WMO3_BELOWSUR)
          k = count (inv(i)%pa == inv(:)%pa            &
               .and. inv(i)%ti == inv(:)%ti            &
               .and. inv(i)%en == inv(:)%en            &
               .and. 0         /= inv(:)%lv% levelvalue)
          !-------------------------------------------------------
          ! field has more levels than any one before: recalculate
          !-------------------------------------------------------
          if (k > ns) then
            ns = k
            if (ns /= size(d)) then
              deallocate (d)
              allocate (d(ns))
            endif
            d = 0
            k = 0
            do j=i,size(inv)
              if     (inv(i)%pa == inv(j)%pa &
                .and. inv(i)%ti == inv(j)%ti &
                .and. inv(i)%en == inv(j)%en ) then
                !-----------------------------------------------
                ! Skip level==0 here, since it likely equals T_G
                !-----------------------------------------------
                if (inv(j)%lv% levelvalue == 0) cycle
                k=k+1
                d(k)=inv(j)%lv% levelvalue
              endif
            end do
            ns = k
          endif
        end select
      end do
      if (ns>0) d(1:ns) = d (index (d(1:ns)))
      !--------------------------------------
      ! extract parameters from GRIB sections
      !--------------------------------------
      nx = 0; ny = 0; ni = 0; lat1 = 0; lon1 = 0
      dx = 0._wp; dy = 0._wp; latr = 0; lonr = 0
      select case (grib% isec2 (1))
      case (WMO6_GAUSSIAN)
        nx   = grib% gauss%  ni
        ny   = grib% gauss%  nj
        if (nx /= 2*ny) then
           write(0,*) "read_atm_grid: gauss: nx, ny =", nx, ny
        end if
        lon1 = grib% gauss% lon_first
        if (lon1 /= 0) call finish ("read_atm_grid","gaussian, but lon1/=0")
!       scan = grib% gauss% scan_mode               ! Retain scanmode from GRIB
        scan = WMO8_J_POSITIVE                      ! DWD internal use: +j
        if (present (scanmode)) scan = scanmode
        select case (scan)
        case (WMO8_J_POSITIVE)
        case (0)
          call message ("read_atm_grid","scanmode 0 for gaussian grid!")
          call finish ('read_atm_grid','not fully implemented, aborting.')
        case default
          write (0,*) "scanmode =", scan
          call finish ('read_atm_grid','invalid scanmode')
        end select
      case (WMO6_LATLON, WMO6_ROTLL)
        nx   = grib% latlon% ni
        ny   = grib% latlon% nj
        if (ibits(grib%latlon% increments, 7, 1) == 1) then ! increments given
          dx = grib% latlon% di
          dy = grib% latlon% dj
        else
          dx = (real(grib%latlon%lon_last,wp) - grib%latlon%lon_first)/(nx-1)
          dy = (real(grib%latlon%lat_last,wp) - grib%latlon%lat_first)/(ny-1)
        endif
        if (present (scanmode)) then
           if (scanmode /= WMO8_J_POSITIVE) then    ! Scan +j is DWD default
              call finish ('read_atm_grid','scanmode only +j for latlon grids')
           end if
        end if
        select case (grib% latlon% scan_mode)
        case (WMO8_J_POSITIVE)
        case (0)
          dy = abs (dy)
          latx                  = grib%latlon%lat_last
          grib%latlon%lat_last  = grib%latlon%lat_first
          grib%latlon%lat_first = latx
        case default
          call finish ('read_atm_grid','invalid scan_mode')
        end select
        lat1 = grib% latlon% lat_first
        lon1 = grib% latlon% lon_first
        latr = grib% latlon% lat_rot
        lonr = grib% latlon% lon_rot
      case (DWD6_ICOSAHEDRON)
        ni   = grib% tri%    ni
        nd   = grib% tri%    nd
      case (DWD6_ICON)
        ni   = grib% tri%    ni
      case default
        call finish ('read_atm_grid','unknown representation')
      end select
    endif
    !--------------------------------------
    ! broadcast to other processor elements
    !--------------------------------------
    if (p_parallel) then
      call p_bcast     (irec,    p_readgrib)
      call p_bcast     (nx,      p_readgrib)
      call p_bcast     (ny,      p_readgrib)
      call p_bcast     (ni,      p_readgrib)
      call p_bcast     (nv,      p_readgrib)
      call p_bcast     (ke1,     p_readgrib)
      call p_bcast     (ns,      p_readgrib)
      call p_bcast     (lat1,    p_readgrib)
      call p_bcast     (lon1,    p_readgrib)
      call p_bcast     (latr,    p_readgrib)
      call p_bcast     (lonr,    p_readgrib)
      call p_bcast     (dx,      p_readgrib)
      call p_bcast     (dy,      p_readgrib)
      call p_bcast     (nd,      p_readgrib)
      call p_bcast     (levtyp,  p_readgrib)
      call p_bcast_ptr (a,       p_readgrib)
      call p_bcast_ptr (b,       p_readgrib)
      call p_bcast_ptr (d,       p_readgrib)
      call p_bcast     (scan,    p_readgrib)
      call p_bcast     (ivctype, p_readgrib)
      call p_bcast     (ra,      p_readgrib)
      call p_bcast     (vc,      p_readgrib)
    endif
    !---------------------
    ! construct basic grid
    !---------------------
    gridtype = inv(1)%gr% gridtype
    uuid_v   = vc% vc_uuid
    select case (gridtype)
    case (DWD6_ICOSAHEDRON)
      select case (levtyp)
      case (WMO3_ISOBARIC)
        call construct (grid, gridtype = gridtype, &
                                   ny  = ny,       &
                                   nx  = nx,       &
                                   ni  = ni,       &
                                levtyp = levtyp,   &
                                   akf = a,        &
                                   bkf = b,        &
                                    ns = ns,       &
                                   dbs = d,        &
                                nproc1 = nproc1,   &
                                nproc2 = nproc2,   &
                                  comm = comm)
      case default
        call construct (grid, gridtype = gridtype, &
                                   ny  = ny,       &
                                   nx  = nx,       &
                                   ni  = ni,       &
                                levtyp = levtyp,   &
                                   ak  = a,        &
                                   bk  = b,        &
                                    ns = ns,       &
                                   dbs = d,        &
                                nproc1 = nproc1,   &
                                nproc2 = nproc2,   &
                                  comm = comm)
      end select
    case (WMO6_LATLON, WMO6_ROTLL)
      select case (levtyp)
      case (WMO3_ISOBARIC)
      call construct (grid, gridtype = gridtype,     &
                                  ny = ny,           &
                                  nx = nx,           &
                                 akf = a,            &
                                 bkf = b,            &
                              levtyp = levtyp,       &
                                 lo1 = lon1/1000._wp,&
                                 la1 = lat1/1000._wp,&
                                  di = dx  /1000._wp,&
                                  dj = dy  /1000._wp,&
                                 lor = lonr/1000._wp,&
                                 lar = latr/1000._wp,&
                                  ns = ns,           &
                                 dbs = d,            &
                              nproc1 = nproc1,       &
                              nproc2 = nproc2,       &
                                comm = comm)
      case default
      !----------------------
      ! arakawa C for HRM, LM
      !----------------------
!     if (dx * nx < 359000 .and. dy * ny < 179000) then     ! ??????
      if (ra% irefatm >= 0 .or. ivctype > 0) then
        arakawa = 'C'
      else
        arakawa = 'A'
      endif

      if (ivctype == 0 .or. levtyp == WMO3_HHYBRID) then
        !----
        ! HRM
        !----
        call construct (grid, gridtype = gridtype,     &
                                    ny = ny,           &
                                    nx = nx,           &
                                   ak  = a,            &
                                   bk  = b,            &
                                levtyp = levtyp,       &
                                uuid_v = uuid_v,       &
                                   lo1 = lon1/1000._wp,&
                                   la1 = lat1/1000._wp,&
                                    di = dx  /1000._wp,&
                                    dj = dy  /1000._wp,&
                                   lor = lonr/1000._wp,&
                                   lar = latr/1000._wp,&
                                    ns = ns,           &
                                   dbs = d,            &
                                nproc1 = nproc1,       &
                                nproc2 = nproc2,       &
                                  comm = comm,         &
                               arakawa = arakawa       )

      else if (levtyp == WMO3_GENV .and. ivctype == IVCTYPE_ICON) then
        !---------------------------------------
        ! Generalized vertical coordinate (ICON)
        !---------------------------------------
        call construct (grid, gridtype = gridtype,     &
                                    ny = ny,           &
                                    nx = nx,           &
                                    ke = ke1-1,        &
                                levtyp = levtyp,       &
                               ivctype = ivctype,      &
                                uuid_v = uuid_v,       &
                                   lo1 = lon1/1000._wp,&
                                   la1 = lat1/1000._wp,&
                                    di = dx  /1000._wp,&
                                    dj = dy  /1000._wp,&
                                   lor = lonr/1000._wp,&
                                   lar = latr/1000._wp,&
                                    ns = ns,           &
                                   dbs = d,            &
                                nproc1 = nproc1,       &
                                nproc2 = nproc2,       &
                                  comm = comm,         &
                               arakawa = arakawa       )
      else
        !------
        ! COSMO
        !------
        call construct (grid, gridtype = gridtype,      &
                                    ny = ny,            &
                                    nx = nx,            &
                                    ke = ke1-1,         &
                                levtyp = levtyp,        &
                                refatm = ra,            &
                                vcoord = vc,            &
                               ivctype = ivctype,       &
                                uuid_v = uuid_v,        &
                                   lo1 = lon1/1000._wp, &
                                   la1 = lat1/1000._wp, &
                                    di = dx  /1000._wp, &
                                    dj = dy  /1000._wp, &
                                   lor = lonr/1000._wp, &
                                   lar = latr/1000._wp, &
                                    ns = ns,            &
                                   dbs = d,             &
                                nproc1 = nproc1,        &
                                nproc2 = nproc2,        &
                                  comm = comm,          &
                               arakawa = arakawa        )
      endif
      end select
    case (WMO6_GAUSSIAN)
      select case (levtyp)
      case (WMO3_ISOBARIC)
      call construct (grid, gridtype = gridtype,     &
                                  ngl= ny,           &
                                  nx = nx,           &
                                 akf = a,            &
                                 bkf = b,            &
                              levtyp = levtyp,       &
                                  ns = ns,           &
                                 dbs = d,            &
                              nproc1 = nproc1,       &
                              nproc2 = nproc2,       &
                                comm = comm,         &
                            scanmode = scan)
      case default
      call construct (grid, gridtype = gridtype,     &
                                  ngl= ny,           &! was ny=ny?
                                  nx = nx,           &
                                 ak  = a,            &
                                 bk  = b,            &
                              levtyp = levtyp,       &
                                  ns = ns,           &
                                 dbs = d,            &
                              nproc1 = nproc1,       &
                              nproc2 = nproc2,       &
                                comm = comm,         &
                            scanmode = scan)
      end select
    case (DWD6_ICON)
      nx = gtmp% nxny
!     ny = nx
!     allocate (lon (nx,1,1), lat (ny,1,1))
!     if (associated (gtmp% icongrid)) then
!        lon(:,:,1) = gtmp% icongrid% patch% cells% center(:,:)% lon * d2r
!        lat(:,:,1) = gtmp% icongrid% patch% cells% center(:,:)% lat * d2r
!     else
!        lon = -HUGE(0._wp)
!        lat = -HUGE(0._wp)
!     end if
      !------------------------------------------------------------
      ! Note: the current ICON grid decomposition requires nproc2=1
      !------------------------------------------------------------
      if (present (nproc1) .and. present (nproc2)) then
        nproc = max (nproc1*nproc2, 1)
      else
        nproc = dace% npe
      end if
      select case (levtyp)
      case (WMO3_ISOBARIC)
!       call finish('read_atm_grid','ICON+isobaric not yet tested')
        call construct (grid, gridtype = gridtype,      &
                                    ni = ni,            &
                                    nx = nx,            &
                              icongrid = gtmp% icongrid,&
                                levtyp = levtyp,        &
                                   akf = a,             &
                                   bkf = b,             &
                                    ns = ns,            &
                                   dbs = d,             &
                                nproc1 = nproc,         &
                                nproc2 = 1,             &
                                  comm = comm)
      case (WMO3_HHYBRID)
        call construct (grid, gridtype = gridtype,      &
                                    ni = ni,            &
                                    nx = nx,            &
                              icongrid = gtmp% icongrid,&
                                levtyp = levtyp,        &
                                   ak  = a,             &
                                   bk  = b,             &
                                    ns = ns,            &
                                   dbs = d,             &
                                nproc1 = nproc,         &
                                nproc2 = 1,             &
                                  comm = comm)
      case (WMO3_GENV)
!print *, "Before construct ICON/GENV"
        call construct (grid, gridtype = gridtype,      &
                                    ni = ni,            &
                                    nx = nx,            &
                                    ke = ke1-1,         &
                              icongrid = gtmp% icongrid,&
                                levtyp = levtyp,        &
                               ivctype = ivctype,       &
                                uuid_v = uuid_v,        &
                                    ns = ns,            &
                                   dbs = d,             &
                                nproc1 = nproc,         &
                                nproc2 = 1,             &
                                  comm = comm)
        ! Do we need to pass additional grid information?
      case default
        if (dace% lpio) then
          write(6,'(/a/)') "   setting up grid without vertical coordinate !"
        endif

        call construct (grid, gridtype = gridtype,      &
                                    ni = ni,            &
                                    nx = nx,            &
                                    ke = 0,             &
                              icongrid = gtmp% icongrid,&
                                levtyp = 0,             &
                                    ns = ns,            &
                                   dbs = d,             &
                                nproc1 = nproc,         &
                                nproc2 = 1,             &
                                  comm = comm)
      end select
!     deallocate (lon, lat)
      !--------------------------------------
      ! Add metadata for ICON horizontal grid
      !--------------------------------------
      grid% icongrid => gtmp% icongrid
    case default
      call finish('read_atm_grid','unknown gridtype')
    end select

    !----------------------------------------------------
    ! set permutation indices if gridpoints are reordered
    ! to fit LETKF coarse grid
    !----------------------------------------------------
    if (present (rni)) then
      if (.not. present (rf)) call finish ('read_atm_grid','rf is missing')
      if (rni > 1 .or. rf > 1) then
        if (dace% lpio) write(6,'(/a,2i6/)') &
             'read_atm_grid: set up reduced grid; rni, rf =', rni, rf
        if (.not. present (g_coarse)) &
             call finish ('read_atm_grid','g_coarse is missing')
        call reduced_grid (rf, rni, nzr, grid, g_coarse)
      endif
    endif
    if (present (g_coarse)) then
      if (associated (g_coarse)) then
        if (gridtype == DWD6_ICON) &
          call distribute_icon_gme (grid, g_coarse)
      endif
    endif
    !-----------------------------
    ! Derive number of full levels
    !-----------------------------
    if (ke1 > 0) then
      ke = ke1 - 1
    else
      ke = max (nv/2-1, 0)
    end if

    !----------------------------------------
    ! open optional file for invariant fields
    !----------------------------------------
    if (linvar) then
      if (dace% lpio) then
        write(6,'(/a,a /)') 'read_atm_grid: invar      = ',trim(invar)
      endif
      nullify (iinv)
      call get_inventory (iinv, invar, pio=p_readgrib)
      if (.not.associated (iinv))                                  &
         call finish('read_atm_grid','no inventory: '//trim (invar))
      if (size(iinv) == 0)                                         &
         call finish('read_atm_grid','no inventory: '//trim (invar))
      if (dace% lpio) then
        write(6,'()')
        call print_inventory (iinv, first=.true.)
        call print_inventory (iinv, first=.true., liname=.true.)
        write(6,'()')
      endif
      if (dace% pe == p_readgrib) then
        call close_gribfile (grib)
        call open_gribfile  (grib, invar, 'r')
      endif
    else
      iinv => inv
    endif
    if (associated (grid% icongrid)) &
         call verify_grid_num (iinv, grid% icongrid% grid_num)
    !-------------------------------
    ! read additional surface fields
    !-------------------------------
    if (irec/=0) then
!print *,"read_atm_grid: irec,inv(irec)%lv%leveltype =", irec,inv(irec)%lv%leveltype
!     select case (inv(irec)%lv% leveltype)
      select case (levtyp)
      case (WMO3_HHYBRID, WMO3_GENV)
        !-------------------------------------
        ! For ICON hybrid coordinates, try HHL
        !-------------------------------------
        call allocate (grid, 'hsurf')
        if (lhhl) then
          call allocate (grid, 'hhl')
          call read_multi (grid% hhl, iinv, grid, grib, name='hhl',     &
                           pio=p_readgrib, pmode=SCATTER, ierr=ierr,    &
                           scanmode=WMO8_J_POSITIVE, member=member      )
          if (ierr == IE_MISSING) then
            call read_multi (grid% hhl, iinv, grid, grib, name='hsurf', &
                           pio=p_readgrib, pmode=SCATTER, ierr=ierr,    &
                           scanmode=WMO8_J_POSITIVE, member=member      )
          end if
          if (ierr /= 0) then
            if (index (lopt,'hhl') /= 0) then
              lhhl = .false.
              call deallocate (grid, 'hhl')
            else
              call error ('read_atm_grid','hhl',ierr)
            endif
          else
            if (p_parallel) then
              x3 => grid% hsurf(:,:,1,:)
              call gather_level (x3, grid% hhl(:,:,ke+1,:), grid% dc, p_readgrib)
              call p_bcast      (grid% hsurf, p_readgrib)
            else
              grid% hsurf(:,:,1,:) = grid% hhl(:,:,ke+1,:)
            end if
          endif
        endif
        if (.not. lhhl) then
          call read_single (grid% hsurf, iinv, grid, grib, name='hsurf', &
                            pio=p_readgrib, pmode=BCAST, ierr=ierr,      &
                            scanmode=WMO8_J_POSITIVE, member=member      )
          if (ierr /= 0) then
            if (index (lopt,'hsurf') /= 0) then
              call deallocate (grid, 'hsurf')
            else
              call error ('read_atm_grid','hsurf',ierr)
            endif
          endif
        endif
        if (associated (grid% hsurf)) then
          call allocate (grid, 'geosp')
          grid% geosp = grid% hsurf * gacc
        endif
        if (lhhl .and. model == 'COSMO') then
           call cosmo_ref_atm (grid, lnew_hhl=.false.)
        end if
      case (WMO3_HYBRID, WMO3_HYBRIDB, WMO3_ISOBARIC)
        if (ivctype /= 0 .and. levtyp /= WMO3_ISOBARIC) then
          !-----------------------------------------------------------
          ! read HSURF field for COSMO, calculate reference atmosphere
          !-----------------------------------------------------------
          call allocate (grid, 'hsurf')
          call read_single (grid% hsurf, iinv, grid, grib, name='hsurf', &
                            pio=p_readgrib, pmode=BCAST, ierr=ierr,      &
                            scanmode=WMO8_J_POSITIVE, member=member      )
          if (ierr == IE_MISSING) then
            call read_single (grid% hsurf, iinv, grid, grib, name='geosp', &
                            pio=p_readgrib, pmode=BCAST, ierr=ierr,        &
                            scanmode=WMO8_J_POSITIVE, member=member        )
            grid% hsurf = grid% hsurf / gacc
          endif
          if (ierr /= 0) then
            if (index (lopt,'hsurf') /= 0) then
              call deallocate (grid, 'hsurf')
            else
              call error ('read_atm_grid','hsurf',ierr)
            endif
          endif
          if (associated (grid% hsurf)) then
            call allocate (grid, 'geosp')
            grid% geosp = grid% hsurf * gacc
          endif
          if (lhhl) call cosmo_ref_atm (grid, lnew_hhl=lhhl)
        else
          !------------------------------------------------
          ! read geopotential height at the surface [m2/s2]
          !------------------------------------------------
          call allocate (grid, 'geosp')
          call read_single (grid% geosp, iinv, grid, grib, name='geosp', &
                            pio=p_readgrib, pmode=BCAST, ierr=ierr,      &
                            scanmode=WMO8_J_POSITIVE, member=member      )
          if (ierr == IE_MISSING) &
            call read_multi(grid% geosp, iinv, grid, grib, name='geosp',     &
                            pio=p_readgrib, pmode=BCAST, kf=ke+1, ierr=ierr, &
                            scanmode=WMO8_J_POSITIVE, member=member          )
          if (ierr == IE_MISSING) then
            call read_multi(grid% geosp, iinv, grid, grib, name='geosp',   &
                            pio=p_readgrib, pmode=BCAST, kf=ke, ierr=ierr, &
                            scanmode=WMO8_J_POSITIVE, member=member        )
          endif
          if (ierr == IE_MISSING .and. lgeosp) then
            if (dace% lpio) print *, "Fallback: trying to read geosp as FI(k=1)"
            call read_multi(grid% geosp, iinv, grid, grib, name='geoh',    &
                            pio=p_readgrib, pmode=BCAST, kf=1, ierr=ierr,  &! IFS
                            scanmode=WMO8_J_POSITIVE, member=member        )
          endif
          if (ierr /= 0) then
            if (lgeosp) call error ('read_atm_grid','geosp (orography)',ierr)
            grid% geosp = 0._wp
            if (levtyp == WMO3_ISOBARIC) call deallocate (grid, 'geosp')
          endif
!         if (ierr == 0 .and. lgeosp) then
!           call allocate (grid, 'hsurf')
!           grid% hsurf = grid% geosp / gacc    ! Emulate hsurf
!         end if
        endif
      end select
    else
      !--------------------------------------------
      ! unknown level type, still try to read hsurf
      !--------------------------------------------
      call allocate (grid, 'hsurf')
      call read_single (grid% hsurf, iinv, grid, grib, name='hsurf', &
                        pio=p_readgrib, pmode=BCAST, ierr=ierr,      &
                        scanmode=WMO8_J_POSITIVE, member=member      )
      if (ierr /= 0) then
        if (index (lopt,'hsurf') /= 0) then
          call deallocate (grid, 'hsurf')
        else
          call error ('read_atm_grid','hsurf',ierr)
        endif
      endif
      if (associated (grid% hsurf)) then
        call allocate (grid, 'geosp')
        grid% geosp = grid% hsurf * gacc
      endif
    endif
    !-----------------
    ! derive model top
    !-----------------
    if (lhhl) then
      call set_ptopf (grid)
    else
      grid% ptopf = 0._wp
    endif
    !-----------------------------------------
    ! estimate unperturbed model-level heights
    !-----------------------------------------
    if (lhhl) then
       call set_zlev_ref (grid)
    end if
    !--------------------------------
    ! set some pressure level indices
    !--------------------------------
    call set_plev_indices (grid)
    !------------------------------------------------------------------------
    ! read surface fields (independent from model vertical coordinate system)
    ! land sea mask:
    !------------------------------------------------------------------------
    if (llsm) then
      call allocate (grid, 'lsm')
      call read_single (grid% lsm, iinv, grid, grib, name='lsm', &
                        pio=p_readgrib, pmode=BCAST, ierr=ierr,  &
                        scanmode=WMO8_J_POSITIVE, member=member  )
      if (ierr /= 0) then
        if (index (lopt,'lsm') /= 0) then
          call deallocate (grid, 'lsm')
        else
          call error ('read_atm_grid','lsm (land-sea mask)',ierr)
        endif
      endif
    endif
    !----------------------
    ! soiltyp (if present):
    !----------------------
    call allocate (grid, 'soiltyp')
    call read_single (grid% soiltyp, iinv, grid, grib, name='soiltyp', &
                      pio=p_readgrib, pmode=BCAST, ierr=ierr,          &
                      scanmode=WMO8_J_POSITIVE, member=member          )
    if (ierr /= 0) call deallocate (grid, 'soiltyp')
    !------------------------------
    ! flake meta data (if present):
    !------------------------------
    call allocate (grid, 'fr_lake')
    call read_single (grid% fr_lake, iinv, grid, grib, name='fr_lake', &
                      pio=p_readgrib, pmode=BCAST, ierr=ierr,          &
                      scanmode=WMO8_J_POSITIVE, member=member          )
    if (ierr /= 0) call deallocate (grid, 'fr_lake')

    call allocate (grid, 'depth_lk')
    call read_single (grid% depth_lk, iinv, grid, grib, name='depth_lk', &
                      pio=p_readgrib, pmode=BCAST, ierr=ierr,            &
                      scanmode=WMO8_J_POSITIVE, member=member            )
    if (ierr /= 0) call deallocate (grid, 'depth_lk')
    !--------------------
    ! sso_* (if present):
    !--------------------
    call allocate (grid, 'sso_stdh')
    call read_single (grid% sso_stdh,  iinv, grid, grib, name='sso_stdh',  &
                      pio=p_readgrib, pmode=SCATTER, ierr=ierr,            &!
                      scanmode=WMO8_J_POSITIVE, member=member              )
    if (ierr /= 0) call deallocate (grid, 'sso_stdh')
    call allocate (grid, 'sso_gamma')
    call read_single (grid% sso_gamma, iinv, grid, grib, name='sso_gamma', &
                      pio=p_readgrib, pmode=SCATTER, ierr=ierr,            &!
                      scanmode=WMO8_J_POSITIVE, member=member              )
    if (ierr /= 0) call deallocate (grid, 'sso_gamma')
    call allocate (grid, 'sso_theta')
    call read_single (grid% sso_theta, iinv, grid, grib, name='sso_theta', &
                      pio=p_readgrib, pmode=SCATTER, ierr=ierr,            &!
                      scanmode=WMO8_J_POSITIVE, member=member              )
    if (ierr /= 0) call deallocate (grid, 'sso_theta')
    call allocate (grid, 'sso_sigma')
    call read_single (grid% sso_sigma, iinv, grid, grib, name='sso_sigma', &
                      pio=p_readgrib, pmode=SCATTER, ierr=ierr,            &!
                      scanmode=WMO8_J_POSITIVE, member=member              )
    if (ierr /= 0) call deallocate (grid, 'sso_sigma')
    call allocate (grid, 'emis_rad')
    call read_single (grid% emis_rad,  iinv, grid, grib, name='emis_rad',  &
                      pio=p_readgrib, pmode=SCATTER, ierr=ierr,            &!
                      scanmode=WMO8_J_POSITIVE, member=member              )
    if (ierr /= 0) call deallocate (grid, 'emis_rad')
    !-------------------------------------------------
    ! COSMO specific external parameters (if present):
    !-------------------------------------------------
    call allocate (grid, 'prs_min')
    call read_single (grid% prs_min,   iinv, grid, grib, name='prs_min', &
                      pio=p_readgrib, pmode=SCATTER, ierr=ierr,          &!
                      scanmode=WMO8_J_POSITIVE, member=member            )
    if (ierr /= 0) call deallocate (grid, 'prs_min')
    call allocate (grid, 'skc')
    call read_single (grid% skc,       iinv, grid, grib, name='skc',     &
                      pio=p_readgrib, pmode=SCATTER, ierr=ierr,          &!
                      scanmode=WMO8_J_POSITIVE, member=member            )
    if (ierr /= 0) call deallocate (grid, 'skc')
    !--------
    ! cleanup
    !--------
    if (linvar) then
      if (associated (iinv)) deallocate (iinv)
    endif
    if (.not. present (invt)) then
      deallocate (inv)
    else
      invt => inv
    endif
    if (p_readgrib == dace% pe) call close_gribfile (grib)
    if (associated (a))      deallocate (a)
    if (associated (b))      deallocate (b)
    if (associated (d))      deallocate (d)
    call destruct  (vc)
    call destruct  (gtmp)
  end subroutine read_atm_grid
!==============================================================================
  subroutine verify_grid_num (inv, grid_num)
    type(t_inventory),pointer    :: inv(:)   ! inventory
    integer          ,intent(in) :: grid_num ! numberOfGridUsed

    integer :: mig, mag  ! min/max of numberOfGridUsed

    if (grid_num <= 0 .or. all (inv(:)% gr% gridtype /= DWD6_ICON)) return
    mig = minval (inv(:)% gr% grid_num, inv(:)% gr% gridtype == DWD6_ICON)
    mag = maxval (inv(:)% gr% grid_num, inv(:)% gr% gridtype == DWD6_ICON)
    if (mig /= mag .or. mig /= grid_num) then
       if (dace% lpio) then
          write(6,*) "Incompatible unstructured horizontal grids:"
          write(6,*) "numberOfGridUsed (GRIB):", mig, "...", mag
          write(6,*) "<> number_of_grid_used =", grid_num
       end if
       call finish ("verify_grid_num","multiple or incompatible ICON grids")
    end if
  end subroutine verify_grid_num
!==============================================================================

  subroutine read_single (x, inv, grd, grib, pos, code, table, dis, cat, num, &
                          name, runtype, time, reftime, pio, pmode,           &
                          isp, ierr, scanmode, member, level, expid           )
  !----------------------------------------------------------------------------
  ! Read a single level field from a GRIB file.
  ! Distribute over processor elements.
  !----------------------------------------------------------------------------
  real(wp)          ,intent(out)            :: x(:,:,:,:) ! field to read
  type(t_inventory) ,intent(in)             :: inv(:)     ! gribfile inventory
  type(t_grid)      ,intent(in)             :: grd        ! grid information
  type(t_grib1)     ,intent(inout)          :: grib       ! GRIB record buffer
  integer           ,intent(in)   ,optional :: pos        ! position in inv
  integer           ,intent(inout),optional :: code       ! grib code to get
  integer           ,intent(inout),optional :: table      ! grib table to use
  integer           ,intent(inout),optional :: dis        ! grib2 discipline
  integer           ,intent(inout),optional :: cat        ! grib2 param.category
  integer           ,intent(inout),optional :: num        ! grib2 param.number
  character(len=*)  ,intent(in)   ,optional :: name       ! name of field to get
  character(len=*)  ,intent(in)   ,optional :: runtype    ! forecast or analysis
  type (t_time)     ,intent(in)   ,optional :: time       ! verification time
  type (t_time)     ,intent(in)   ,optional :: reftime    ! reference    time
  integer           ,intent(in)   ,optional :: pio        ! pe to use for reading
  integer           ,intent(in)   ,optional :: pmode      ! parallel mode
  integer           ,intent(in)   ,optional :: isp        ! grid spacing
  integer           ,intent(out)  ,optional :: ierr       ! error return code
  integer           ,intent(in)   ,optional :: scanmode   ! scan mode
  integer           ,intent(in)   ,optional :: member     ! ensemble member
  integer           ,intent(in)   ,optional :: level      ! levelvalue
  integer           ,intent(in)   ,optional :: expid      ! experiment id

    integer           :: id            ! grib record number to decode
    integer           :: lmode         ! local copy of pmode
    integer           :: lpio          ! local copy of pio
    logical           :: lio           ! true for I/O processor
!   integer           :: lsp           ! grid spacing
    real(wp) ,pointer :: xg (:,:,:)    ! temporary global field

    lmode = 0;          if(present(pmode)) lmode = pmode
    id    = 0;          if(present(pos))   id    = pos
!   lsp   = 1;          if(present(isp))   lsp   = isp    ! currently not used
    lpio  = p_readgrib; if(present(pio))   lpio  = pio
    lio   = (lpio == dace% pe)

    nullify (xg)
    !----------------------------------------------------------------
    ! find index of single level field in GRIB file inventory, decode
    !----------------------------------------------------------------
    call decode_single (id, inv, ierr, code, table, dis, cat, num,      &
                        name, runtype, time, reftime, member, level,    &
                        grd, grib, scanmode, xg, decode=lio, expid=expid)

    !---------------
    ! error handling
    !---------------
    if(present(ierr)) then
      if (ierr /= 0) return
    endif

    !---------------------------
    ! scatter / broadcast fields
    !---------------------------
    select case (lmode)
    case (SCATTER)
      call scatter_level (xg, x(:,:,1,:), grd%dc, lpio)
    case (BCAST)
      if (lio) x(:,:,1,:) = xg
      call p_bcast (x(:,:,:,:), lpio)
    case (0)
      x(:,:,1,:) = xg
    end select
    if (associated(xg))  deallocate (xg)

  end subroutine read_single
!------------------------------------------------------------------------------
  subroutine decode_single (id, inv, ierr, code, table, dis, cat, num,   &
                            name, runtype, time, reftime, member, level, &
                            grd, grib, scanmode, x, x4, decode, expid)
  !-----------------------------------------------------------
  ! find a single level field in a GRIB file inventory,
  ! optionally decode the record
  ! Optional error return code:
  ! ierr = -1 : field not found
  ! ierr = -2 : field is ambiguous
  !-----------------------------------------------------------
  integer           ,intent(inout)        :: id         ! grib record index
  type(t_inventory) ,intent(in)           :: inv(:)     ! gribfile inventory
  integer           ,intent(out),optional :: ierr       ! error return code
  !--------------------
  ! selection criteria:
  !--------------------
  integer         ,intent(inout),optional :: code       ! grib code to get
  integer         ,intent(inout),optional :: table      ! grib table to use
  integer         ,intent(inout),optional :: dis        ! grib2 discipline
  integer         ,intent(inout),optional :: cat        ! grib2 param.category
  integer         ,intent(inout),optional :: num        ! grib2 param.number
  character(len=*),intent(in)   ,optional :: name       ! name of field to get
  character(len=*),intent(in)   ,optional :: runtype    ! eg forecast, analysis
  type (t_time)   ,intent(in)   ,optional :: time       ! verification time
  type (t_time)   ,intent(in)   ,optional :: reftime    ! reference    time
  integer         ,intent(in)   ,optional :: member     ! ensemble member
  integer         ,intent(in)   ,optional :: level      ! levelvalue
  integer         ,intent(in)   ,optional :: expid      ! experiment id
  !-------
  ! decode
  !-------
  type(t_grid)    ,intent(in)   ,optional :: grd        ! grid information
  type(t_grib1)   ,intent(inout),optional :: grib       ! GRIB record buffer
  integer         ,intent(in)   ,optional :: scanmode   ! scan mode
  real(wp)        ,pointer      ,optional :: x(:,:,:)   ! decoded field
  real(wp)        ,pointer      ,optional :: x4(:,:,:,:)! alternative for x
  logical         ,intent(in)   ,optional :: decode     ! decode on this PE

    integer           :: i             ! grib record number loop index
    logical           :: ldec          ! local copy of decode
    logical           :: ok1
    logical           :: ok2

    !-----------------------------------------------------------------
    ! if position in inventory is not prescribed search suitable entry
    !-----------------------------------------------------------------
    if (id <= 0) then
      !-----------------------------------------------------------------
      ! for given name/process search gribtable for table/code/leveltype
      !-----------------------------------------------------------------
      ok1 = .false.
      ok2 = .false.
      if (present(code) .and. present(table)) then
        if (code /= 255 .and. table/= 255) then
          ok1 = .true.
        endif
      endif
      if (present(dis) .and. present(cat) .and. present(num)) then
        if (dis /= 255 .and. cat/= 255 .and. num/=255) then
          ok2 = .true.
        endif
      endif
      if (present (x)) nullify (x)
      !----------------------------
      ! loop over inventory entries
      !----------------------------
      id = 0
      do i=1,size(inv)
        !----------------------------------------------
        ! check correspondance with optional parameters
        !----------------------------------------------
        if (present(runtype)) then
          if (runtype /= '' .and. inv(i)%pa% runtype /= runtype) cycle
        endif
        if (present(expid)) then
          if (expid > 0 .and. inv(i)%pa% expid /= expid)         cycle
        endif
        if (present(time)) then
          if (inv(i)%ti% ver_time /= time) cycle
        endif
        if (present(reftime)) then
          if (inv(i)%ti% ref_time /= reftime) cycle
        endif
        if (present(member)) then
          if (inv(i)%en% no /= member) cycle
        endif
        if (ok1) then
          if (inv(i)%pa% code /= code)  cycle
          if (inv(i)%pa% table/= table) cycle
        endif
        if (ok2) then
          if (inv(i)%pa% discipline /= dis) cycle
          if (inv(i)%pa% category   /= cat) cycle
          if (inv(i)%pa% number     /= num) cycle
        endif
        if (present(level)) then
          if (inv(i)%lv% levelvalue /= level) cycle
          select case (inv(i)%lv% leveltype)
          case (100,109,110)
          case default
            cycle
          end select
        else
          select case (inv(i)%lv% leveltype)
          case (1:9,102,103,105,111,162,166,208)
          case default
            cycle
          end select
        endif
        if (.not.(ok1.or.ok2).and.present(name)) then
          if (inv(i)%pa% iname /= name) cycle
          if (present(code))  code  = inv(i)%pa% code
          if (present(table)) table = inv(i)%pa% table
          if (present(dis))   dis   = inv(i)%pa% discipline
          if (present(cat))   cat   = inv(i)%pa% category
          if (present(num))   num   = inv(i)%pa% number
        endif
        !-----------------------------------
        ! error return if field is ambiguous
        !-----------------------------------
        if (id /= 0) then
          id = IE_AMBIGUOUS ! -2
          exit
        endif
        id   = i
      end do
      if (id==0) id = IE_MISSING ! -1

      !---------------
      ! error handling
      !---------------
      if(present(ierr)) ierr   = 0
      if (id <= 0) then
        if (present (ierr)) then
          ierr = id
          return
        else
          call error ('read_single', name, id)
        endif
      endif
    endif

    !-------
    ! decode
    !-------
    if (id <= 0) return
    ldec = .true.
    if (present (decode)) ldec = decode
    if (ldec) then
      if (.not.  present (grd)  .or.    &
          .not.  present (grib) .or.    &
          .not. (present (x)    .or.    &
                 present (x4)        )) &
        call finish('decode_single','grd,grib or x not present')

      call pbseek (grib, inv(id)%ds% ptr)
      if (present (x4)) then
        allocate (x4(grd%lbg(1) : grd%ubg(1) , &
                     grd%lbg(2) : grd%ubg(2) , &
                              1 : 1          , &
                     grd%lbg(4) : grd%ubg(4) ) )
        call read_gribex (grib, 'R', x4(:,:,1,:), scanmode=scanmode)
        call permut_read (grd, x4=x4)
        if (inv(id)%pa% factor /= 1._wp)                  &
                               x4 = x4 * inv(id)%pa% factor ! scaling factor
        if (present (x)) x => x4 (:,:,1,:)
      else
        allocate (x (grd%lbg(1) : grd%ubg(1) , &
                     grd%lbg(2) : grd%ubg(2) , &
                     grd%lbg(4) : grd%ubg(4) ) )
        call read_gribex (grib, 'R', x, scanmode=scanmode)
        call permut_read (grd, x=x)
        if (inv(id)%pa% factor /= 1._wp)                &
                               x = x * inv(id)%pa% factor   ! scaling factor
      endif
    endif

  end subroutine decode_single
!==============================================================================
  subroutine read_multi (x, inv, grd, grib, code, table, dis, cat, num,  &
                         name, runtype, time, reftime, month,            &
                         pio, pmode, kf, ierr, scanmode, member, expid   )
  real(wp)          ,intent(out)            :: x(:,:,:,:) ! field to read
  type(t_inventory) ,intent(in)             :: inv(:)     ! gribfile inventory
  type(t_grid)      ,intent(in)             :: grd        ! grid information
  type(t_grib1)     ,intent(inout)          :: grib       ! GRIB record buffer
  character(len=*)  ,intent(in)   ,optional :: name       ! name of field toget
  integer           ,intent(inout),optional :: code       ! grib code to get
  integer           ,intent(inout),optional :: table      ! grib table to use
  integer           ,intent(inout),optional :: dis        ! grib2 discipline
  integer           ,intent(inout),optional :: cat        ! grib2 param.category
  integer           ,intent(inout),optional :: num        ! grib2 param.number
  character(len=*)  ,intent(in)   ,optional :: runtype    ! forecast, analysis
  type (t_time)     ,intent(in)   ,optional :: time       ! verification time
  type (t_time)     ,intent(in)   ,optional :: reftime    ! reference    time
  integer           ,intent(in)   ,optional :: month      ! month of verif. time
  integer           ,intent(in)   ,optional :: pio        ! pe to read from
  integer           ,intent(in)   ,optional :: pmode      ! parallel mode
  integer           ,intent(in)   ,optional :: kf         ! lowest level number
  integer           ,intent(out)  ,optional :: ierr       ! error return code
  integer           ,intent(in)   ,optional :: scanmode   ! scan mode
  integer           ,intent(in)   ,optional :: member     ! ensemble member
  integer           ,intent(in)   ,optional :: expid      ! experiment id

    integer              :: lmode         ! local copy of pmode
    integer              :: lpio          ! local copy of pio
    logical              :: decode        ! true for I/O processor
    integer              :: id            ! inventory index of record
    real(wp) ,pointer    :: xg (:,:,:,:)  ! temporary global field
    integer              :: k1, k2        ! lowest, highest level number
    integer              :: kk, kt, dk, n ! Auxiliaries for level blocking
    integer(i8)          :: lev_size      ! Bytes / level

    lmode  = 0         ;if(present(pmode)) lmode = pmode
    lpio   = dace% pio ;if(present(pio))   lpio  = pio
    decode = (lpio == dace% pe)

    k1     = 1         ;if(present(kf))    k1    = kf
    k2     = size (x,3) - k1 + 1
    dk     = k2 - k1 + 1

    if (lmode == SCATTER .and. io_max_scatter >= 0) then
       lev_size = (grd% ubg(1)-grd% lbg(1)+1) * (grd% ubg(2)-grd% lbg(2)+1)
       lev_size = lev_size                    * (grd% ubg(4)-grd% lbg(4)+1)
       lev_size = lev_size * (storage_size (0._wp) / 8)
       dk       = min (max (int (io_max_scatter / lev_size), 1), dk)
    end if
    !----------------------------
    ! Level blocking (outer loop)
    !----------------------------
    do kk = k1, k2, dk
       kt = min (kk+dk-1, k2)
       n  = kt - kk + 1
       nullify (xg)
    !-------------------------------------------
    ! seek record in GRIB file inventory, decode
    !-------------------------------------------
    call decode_multi (xg, inv, grd, grib, n, decode, id, ierr,         &
                       scanmode, code, table, dis, cat, num, name,      &
                       runtype, time, reftime, month, kk, member, expid )
    !---------------
    ! error handling
    !---------------
    if (present (ierr)) then
      if (ierr /= 0) then
        if (associated (xg)) deallocate (xg)
        return
      end if
    endif

    !--------------------
    ! scatter / broadcast
    !--------------------
FTRACE_BEGIN("read_multi:scatter_level")
    select case (lmode)
    case (SCATTER)
      call scatter_multi (xg, x(:,:,kk:kt,:), grd%dc, lpio)
    case (BCAST)
      if (decode) x = xg
      call p_bcast (x, lpio)
    case (0)
      x = xg
    case default
    end select
FTRACE_END  ("read_multi:scatter_level")
    if (associated (xg)) deallocate (xg)

    end do    ! kk loop (level blocking)

  end subroutine read_multi
!------------------------------------------------------------------------------
  subroutine decode_multi (x, inv, grd, grib, nk, decode, id, ierr, scanmode, &
                           code, table, dis, cat, num, name, runtype,         &
                           time, reftime, month, kf, member, expid, leveltypes)
  !-------------------------------------------------------------
  ! find records of a multilevel field in a GRIB file inventory.
  ! optionally decode the record.
  !-------------------------------------------------------------
  real(wp)          ,pointer                :: x(:,:,:,:) ! field to read
  type(t_inventory) ,intent(in)             :: inv(:)     ! gribfile inventory
  type(t_grid)      ,intent(in)             :: grd        ! grid information
  type(t_grib1)     ,intent(inout)          :: grib       ! GRIB record buffer
  integer           ,intent(in)             :: nk         ! number of levels
  logical           ,intent(in)             :: decode     ! decode on this PE
  integer           ,intent(out)            :: id         ! record index
  integer           ,intent(out)  ,optional :: ierr       ! error return code
  integer           ,intent(in)   ,optional :: scanmode   ! scan mode
  !---------------------
  ! selection parameters
  !---------------------
  integer         ,intent(inout),optional :: code         ! grib code to get
  integer         ,intent(inout),optional :: table        ! grib table to use
  integer         ,intent(inout),optional :: dis          ! grib2 discipline
  integer         ,intent(inout),optional :: cat          ! grib2 param.category
  integer         ,intent(inout),optional :: num          ! grib2 param.number
  character(len=*),intent(in)   ,optional :: name         ! name of field toget
  character(len=*),intent(in)   ,optional :: runtype      ! forecast, analysis
  type (t_time)   ,intent(in)   ,optional :: time         ! verification time
  type (t_time)   ,intent(in)   ,optional :: reftime      ! reference    time
  integer         ,intent(in)   ,optional :: month        ! month of verif. time
  integer         ,intent(in)   ,optional :: kf           ! first level number
  integer         ,intent(in)   ,optional :: member       ! ensemble member
  integer         ,intent(in)   ,optional :: expid        ! experiment id
  integer         ,intent(in)   ,optional :: leveltypes(:)! level types

    !----------------
    ! local variables
    !----------------
    integer              :: i             ! grib record number
    integer              :: k             ! level
    integer              :: j             ! index
    logical, allocatable :: found (:)     ! flag field (levels found)
    integer              :: k1, k2        ! lowest, highest level number
    integer              :: k0            ! vertical index
    logical              :: ok1
    logical              :: ok2
    integer              :: ierror        ! local copy of ierr
    integer              :: nfound        ! # levels found

    id     = 0
    ierror = 0
    k1     = 1    ;if(present(kf))    k1    = kf
    k2     = k1 + nk - 1
    !---------------------------------------
    ! Allocate global array on I/O processor
    !---------------------------------------
    if (decode) then
      allocate (x (grd%lbg(1) : grd%ubg(1) , &
                   grd%lbg(2) : grd%ubg(2) , &
                            1 : k2-k1+1    , &
                   grd%lbg(4) : grd%ubg(4) ) )
    else
      nullify (x)
    end if
    !-----------------------------------------------------------------
    ! for given name/process search gribtable for table/code/leveltype
    !-----------------------------------------------------------------
    ok1 = .false.
    ok2 = .false.
    if (present(code) .and. present(table)) then
      if (code /= 255 .and. table/= 255) then
        ok1 = .true.
      endif
    endif
    if (present(dis) .and. present(cat) .and. present(num)) then
      if (dis /= 255 .and. cat/= 255 .and. num/=255) then
        ok2 = .true.
      endif
    endif
    !----------------------------
    ! loop over inventory entries
    !----------------------------
    allocate (found (k1:k2))
    found = .false.
    do i=1,size(inv)
      !----------------------------------------------
      ! check correspondence with optional parameters
      !----------------------------------------------
      if (present(runtype)) then
        if (runtype /= '' .and. inv(i)%pa% runtype /= runtype) cycle
      endif
      if (present(expid)) then
        if (expid > 0 .and. inv(i)%pa% expid /= expid)         cycle
      endif
      if (present(time)) then
        if (inv(i)%ti% ver_time /= time)    cycle
      endif
      if (present(reftime)) then
        if (inv(i)%ti% ref_time /= reftime) cycle
      endif
      if (present(month)) then
        if (imm(inv(i)%ti% ver_time) /= month) cycle
      end if
      if (present(member)) then
        if (inv(i)%en% no /= member)        cycle
      endif
      if (ok1) then
        if (inv(i)%pa% code  /= code)  cycle
        if (inv(i)%pa% table /= table) cycle
      endif
      if (ok2) then
        if (inv(i)%pa% discipline /= dis) cycle
        if (inv(i)%pa% category   /= cat) cycle
        if (inv(i)%pa% number     /= num) cycle
      endif
      select case (inv(i)%lv% leveltype)
      case (WMO3_ISOBARIC, WMO3_BELOWSUR, WMO3_HYBRID, WMO3_HYBRIDB, &
            WMO3_HHYBRID,  WMO3_GENV)
       if (present (leveltypes)) then
         if (.not. any(leveltypes == inv(i)%lv% leveltype)) cycle
       endif
      case default
        cycle
      end select
      if (.not. (ok1.or.ok2) .and. present(name)) then
        if (inv(i)%pa% iname /= name) cycle
        if (present(code))  code  = inv(i)%pa% code
        if (present(table)) table = inv(i)%pa% table
        if (present(dis))   dis   = inv(i)%pa% discipline
        if (present(cat))   cat   = inv(i)%pa% category
        if (present(num))   num   = inv(i)%pa% number
      endif
      !----------------------------------------------
      ! skip ambiguous t_s if t_so is requested
      ! obsolete for GRIB2 as level type is different
      !----------------------------------------------
      if(present(name)) then
        if(name == 't_so' .and. inv(i)% pa% iname == 't_s') then
          if (dace% lpio) call message('decode_multi',   &
            'skipped ambigous T_S when t_so was requested')
          cycle
        endif
      endif
      !-----------------------
      ! determine level number
      !-----------------------
      k = -1
      select case (inv(i)%lv% leveltype)
      case (WMO3_HYBRID, WMO3_HYBRIDB, WMO3_HHYBRID, WMO3_GENV)
        k = inv(i)%lv% levels(1)
        if (k == 0) k = inv(i)%lv% levels(2)    ! used for hhl from ICON
      case (WMO3_ISOBARIC)
        do j=1,grd% nz
          if (inv(i)%lv% levelvalue * 100._wp == grd% akf(j)) then
            k = j
            exit
          endif
        end do
      case (WMO3_BELOWSUR)
        do j=1,grd% ns
          if (inv(i)%lv% levelvalue == grd% dbs(j)) then
            k = j
            exit
          endif
        end do
        if (k < 0 .and. inv(i)%lv% levelvalue == 0) k = 0  ! for 't_so(0)'
      end select
!     if (k/=0 .and. (k<k1 .or. k>k2)) then    ! allow k=0 for 't_so'
      if (k<k1 .or. k>k2) then
        if (k >= 0) cycle
!       cycle
        ierror = IE_INV_LEVEL; goto 999
      endif
      if (found(k)) then
        ierror = IE_AMBIGUOUS; goto 999
      endif
      !-------------
      ! decode field
      !-------------
      k0 = k-k1+1
      if (decode) then
        call pbseek      (grib, inv(i)%ds% ptr)
        call read_gribex (grib, 'R', x(:,:,k0,:), scanmode=scanmode)
        call permut_read (grd, x(:,:,k0,:))
        if (inv(i)%pa% factor /= 1._wp)                  &
             x(:,:,k0,:) = x(:,:,k0,:) * inv(i)%pa% factor  ! scaling factor
      end if
      found(k) = .true.
      id       = i
    end do
    !---------------
    ! error handling
    !---------------
    nfound = count (found(k1:k2))
    if (nfound == 0) then
       ierror = IE_MISSING
    else if (nfound /= k2-k1+1) then
       ierror = IE_LEVELS
!print *, "decode_multi: found=", count (found(k1:k2)), ' expected=', k2-k1+1
    end if
999 continue
    if (present (ierr)) then
       ierr = ierror
    else
       if (ierror/=0) call error ('read_multi', name, ierror)
    end if

  end subroutine decode_multi
!------------------------------------------------------------------------------
  subroutine error (routine, name, ierror)
  character(len=*) ,intent(in)           :: routine
  character(len=*) ,intent(in) ,optional :: name
  integer          ,intent(in)           :: ierror
    select case (ierror)
    case (IE_MISSING)
      if (present(name)) then
        call finish (routine,'field not found: '//name)
      else
        call finish (routine,'field not found')
      endif
    case (IE_AMBIGUOUS)
      if (present(name)) then
        call finish (routine,'multiple fields found: '//name)
      else
        call finish (routine,'multiple fields found')
      endif
    case (IE_FILE)
      if (present(name)) then
        call finish (routine,' not found in grib file: '//name)
      else
        call finish (routine,' field not found in grib file')
      endif
    case (IE_INV_LEVEL)
      call finish (routine,'invalid level')
    case (IE_LEVELS)
      if (present(name)) then
        call finish (routine,'some levels missing: '//name)
      else
        call finish (routine,'some levels missing')
      endif
    case default
      if (present(name)) then
         write(0,*) trim (name), ", Error code:", ierror
      else
         write(0,*) "Error code:", ierror
      end if
      call finish (routine,'invalid errorcode')
    end select
  end subroutine error
!==============================================================================
  subroutine print_2d_rfield (x, name, minv, maxv)
  real(wp)        ,intent(in)           :: x(:,:)
  character(len=*),intent(in)           :: name
  real(wp)        ,intent(in), optional :: minv, maxv
    integer  :: lb(2), ub(2), i,j
    real(wp) :: maxx, minx, delt
    lb = lbound(x)
    ub = ubound(x)
    minx = minval(x)
    maxx = maxval(x)
    if (present(minv)) minx = minv
    if (present(maxv)) maxx = maxv
    delt = (maxx-minx)/10
    if (delt==0._wp) delt = 1._wp
    print *
    print *, name, ': ',minx,' .. ', maxx
    print *
    do j=ub(2), lb(2), -1
      write(*,'(i3,1x,200i1)') j,(int((x(i,j)-minx)/delt),i=lb(1), ub(1))
    end do
    print *
    if (ub(1)>99) write(*,'(4x,200i1)') (mod(i/100,10),i=lb(1), ub(1))
    if (ub(1)> 9) write(*,'(4x,200i1)') (mod(i/ 10,10),i=lb(1), ub(1))
    if (ub(1)> 0) write(*,'(4x,200i1)') (mod(i    ,10),i=lb(1), ub(1))
    print *
  end subroutine print_2d_rfield
!------------------------------------------------------------------------------
  subroutine print_2d_lfield (x, name)
  logical         ,intent(in)           :: x(:,:)
  character(len=*),intent(in)           :: name
    integer  :: lb(2), ub(2), i,j
    lb = lbound(x)
    ub = ubound(x)
    print *
    print *, name, ': '
    print *
    do j=ub(2), lb(2), -1
      write(*,'(i3,1x,200l1)') j,(x(i,j),i=lb(1), ub(1))
    end do
    print *
    if (ub(1)>99) write(*,'(4x,200i1)') (mod(i/100,10),i=lb(1), ub(1))
    if (ub(1)> 9) write(*,'(4x,200i1)') (mod(i/ 10,10),i=lb(1), ub(1))
    if (ub(1)> 0) write(*,'(4x,200i1)') (mod(i    ,10),i=lb(1), ub(1))
    print *
  end subroutine print_2d_lfield
!------------------------------------------------------------------------------
  subroutine print_3d_rfield (x, name, minv, maxv)
  real(wp)        ,intent(in)           :: x(:,:,:)
  character(len=*),intent(in)           :: name
  real(wp)        ,intent(in), optional :: minv, maxv
    integer k
    do k=lbound(x,3), ubound(x,3)
      call print_2d_rfield(x(:,:,k),                                  &
                           trim(name)//'(:,:,'//trim(string(k))//')', &
                           minv, maxv)
    end do
  end subroutine print_3d_rfield
!------------------------------------------------------------------------------
  pure function string(i) result (c)
  !-----------------------------------------------------------------
  ! generic call : string(i)
  ! returns the left adjusted character representation of an integer
  !-----------------------------------------------------------------
  integer, intent(in) :: i
  character (len=10)  :: c
    integer :: k, l
    logical :: negative
    negative = i < 0
    l = abs(i)
    k = len(c)
    c = ' '
    c(k:k) = '0'
    k = k + 1
    do
      if (l>0) then
        k = k - 1
        if (k<1) then
          c = '***'
          return
        endif
        c(k:k) = char(mod(l,10)+ichar('0'))
        l = l / 10
      else
        exit
      endif
    end do
    if (negative) then
      k = k - 1
      if (k<1) then
        c = '***'
        return
      endif
      c(k:k) = '-'
    endif
    c = adjustl (c)
  end function string
!==============================================================================
  !--------------------------------------------------------
  ! Write the fields of an atmospheric state to a GRIB-file
  !--------------------------------------------------------
  subroutine write_atm_grib (atm, file, mode, grid, ref, edition)
  type (t_atm)     ,intent(inout)        :: atm     ! atmosph. state to write
  character(len=*) ,intent(in)           :: file    ! file to write to
  character        ,intent(in) ,optional :: mode    ! 'w' or 'a' (append)
  character(len=*) ,intent(in) ,optional :: grid    ! grid-parameters to write
  logical          ,intent(in) ,optional :: ref     ! write 'reference' fields
  integer          ,intent(in) ,optional :: edition ! GRIB edition to write

    !----------------------------------------
    ! local variables passed to write_atm_gen
    !----------------------------------------
    type (t_atm) :: la  (1) ! local copy  (array of size 1 to fake ensemble)
    integer      :: pio (1) ! processors to use for I/O
    integer      :: ie      ! element to write on this processor element

    la (1) = atm
    pio(1) = dace% pio
    ie     = 1;   if (.not.dace% lpio) ie = -1

    call write_atm_gen (la, file, pio, ie, mode, grid, ref, edition)

  end subroutine write_atm_grib
!------------------------------------------------------------------------------
  !------------------------------------------------------
  ! Write an ensemble of atmospheric states to GRIB-files
  !------------------------------------------------------
  subroutine write_atm_ens (atm, file, mode, grid, ref, edition)
  type (t_atm)     ,intent(inout)        :: atm(:)  ! ensemble to write
  character(len=*) ,intent(in)           :: file    ! file to write to
  character        ,intent(in) ,optional :: mode    ! 'w' or 'a' (append)
  character(len=*) ,intent(in) ,optional :: grid    ! grid-parameters to write
  logical          ,intent(in) ,optional :: ref     ! write reference fields
  integer          ,intent(in) ,optional :: edition ! GRIB edition to write

    integer              :: i         ! start index of chunk to write at a time
    integer              :: n         ! size        of chunk to write at a time
    integer              :: ie        ! ensemble member index in chunk
    integer              :: ke        ! ensemble member index
    integer              :: nchunk    ! number of chunks to write
    integer              :: mcsize    ! maximum chunk size
    integer              :: nens      ! ensemble size
    integer, allocatable :: pio(:)    ! processor elements to use for writing

    nens = size (atm)
    allocate (pio(nens))
!   call set_iope_ens (pio, nchunk, mcsize, stride=0)   ! Sequential write
!   call set_iope_ens (pio, nchunk, mcsize, stride=+1)  ! Tight I/O pattern
    call set_iope_ens (pio, nchunk, mcsize, stride=-1)  ! Optimized I/O stride
    !-------------------------------------------------------------
    ! call write_atm_gen to write 1 member per chosen PE at a time
    !-------------------------------------------------------------
    do i = 1, nens, mcsize              ! loop over bunches to write
      n  = min (nens-i+1, mcsize)       ! number of members to write
      ie = -1
      ke = -1
      if (any (pio(i:i+n-1) == dace% pe)) then
        ie = minloc (abs (pio(i:i+n-1) - dace% pe),1) ! index within chunk
        ke = i - 1 + ie                               ! write member on this PE
        !print *, "### write_atm_ens: p_pe,ie,ke =", dace% pe, ie, ke
      endif
      if (dace% lpio) then
        write(6,*)
        write(6,*) ' writing ensemble members',i,' to',i+n-1
        write(6,*)
      endif
      call write_atm_gen (atm(i:i+n-1), trim(file)//'.'//char3(ke), &
                          pio(i:i+n-1), ie, mode, grid, ref, edition)
    end do

  end subroutine write_atm_ens
!------------------------------------------------------------------------------
  !---------------------------------------------------
  ! Write atmospheric states to GRIB-files in parallel
  !---------------------------------------------------
  subroutine write_atm_grib_p (atm, file, mode, grid, ref, edition, &
                               pid, ptype, mask, verbose            )
  type (t_atm)     ,intent(inout)        :: atm(:)   ! ensemble to write
  character(len=*) ,intent(in)           :: file(:)  ! file to write to
  character        ,intent(in) ,optional :: mode     ! 'w' or 'a' (append)
  character(len=*) ,intent(in) ,optional :: grid     ! grid-parameters to write
  logical          ,intent(in) ,optional :: ref      ! write reference fields
  integer          ,intent(in) ,optional :: edition  ! GRIB edition to write
  integer          ,intent(in) ,optional :: pid(:)   ! process id of file
  integer          ,intent(in) ,optional :: ptype(:) ! process type of file
  logical          ,intent(in) ,optional :: mask(:)  ! subset to write
  logical          ,intent(in) ,optional :: verbose  ! verbose diagnostics

    integer              :: i         ! start index of chunk to write at a time
    integer              :: ii        ! start index in presence of mask
    integer              :: n         ! size        of chunk to write at a time
    integer              :: ie        ! ensemble member index in chunk
    integer              :: ke        ! ensemble member index
    integer              :: nchunk    ! number of chunks to write
    integer              :: mcsize    ! maximum chunk size
    integer              :: nens      ! ensemble size
    integer, allocatable :: pio(:)    ! processor elements to use for writing
    integer              :: k         ! index
    character(len(file)) :: fname     ! temporary file name
    logical              :: verb      ! local copy of verbose

    !----------------------
    ! check input arguments
    !----------------------
    nens = size (atm)
    if ( nens < 1 ) return
    if ( nens /= size (file) ) &
      call finish ('write_atm_grib_p','dim atm /= dim file')
    if (present(pid)) then
      if ( nens /= size (pid) ) &
        call finish ('write_atm_grib_p','dim atm /= dim pid')
    end if
    if (present(ptype)) then
      if ( nens /= size (ptype) ) &
        call finish ('write_atm_grib_p','dim atm /= dim ptype')
    end if
    verb = .false.; if (present (verbose)) verb = verbose

    if (present (mask)) then
       if (size (mask) /= nens) &
          call finish ('write_atm_grib_p','size (mask) /= size (atm)')
       if (count (mask) == 0) return
    end if
    !------------------------------------------------
    ! workaround for different allocation status
    ! set dummy flag, allocate fields with fill value
    !------------------------------------------------
    do k = 1, size(atm(1)% m)
      if ( any(atm% m(k)% i% alloc) ) then
        do i = 1, nens
          if ( .not. atm(i)% m(k)% i% alloc ) then
            call allocate (atm(i), atm(i)% m(k)% i% name)
            atm(i)% m(k)% ptr      = 0._wp
            atm(i)% m(k)% i% dummy = .true.
          end if
        end do
      end if
    end do

    allocate (pio(nens))
!   call set_iope_ens (pio, nchunk, mcsize, stride=0)   ! Sequential write
!   call set_iope_ens (pio, nchunk, mcsize, stride=+1)  ! Tight I/O pattern
    call set_iope_ens (pio, nchunk, mcsize, stride=-1)  ! Optimized I/O stride
    if (present (mask)) then
      where (.not. mask) pio = -1
    end if
    !-------------------------------------------------------------
    ! call write_atm_gen to write 1 file per chosen PE at a time
    !-------------------------------------------------------------
    do i = 1, nens, mcsize              ! loop over bunches to write
      n  = min (nens-i+1, mcsize)       ! number of states  to write
      ie = -1
      ke = -1
      ii = -1
      do ii = i, i+n-1
        if (pio(ii) >= 0) exit          ! search first state with mask=.true.
      end do
      if (ii >= i+n) cycle
      if (any (pio(i:i+n-1) == dace% pe)) then
        ie = minloc (abs (pio(i:i+n-1) - dace% pe),1) ! index within chunk
        ke = i - 1 + ie                               ! file to write on this PE
        !print *, "### write_atm_ens: p_pe,ie,ke =", dace% pe, ie, ke
      endif
      if (dace% lpio) then
        write(6,'(/,A,I3,A,/)') &
             ' writing ', count (pio(i:i+n-1) >=0 ), ' GRIB files in parallel'
        if (verb) then
          do k = i, i+n-1
            if (pio(k) < 0) cycle
            write(6,'(A,I3,A)') '  pe ', pio(k), ': '//trim(file(k))
          end do
        end if
      endif
      if (ke > 0) then
        if (present(pid))   call set_defaults (process      = pid(ke)  )
        if (present(ptype)) call set_defaults (process_type = ptype(ke))
        fname = file(ke)
      else
        fname = ""
      end if
      call write_atm_gen (atm(ii:i+n-1), trim (fname),               &
                          pio(ii:i+n-1), ie, mode, grid, ref, edition)
    end do

    !------------------------
    ! deallocate dummy fields
    !------------------------
    do i = 1, nens
      do k = 1, size(atm(1)% m)
        if ( atm(i)% m(k)% i% dummy ) then
            call deallocate(atm(i), atm(i)% m(k)% i% name)
            atm(i)% m(k)% i% dummy = .false.
        end if
      end do
    end do

  end subroutine write_atm_grib_p
!------------------------------------------------------------------------------
  !----------------------------------------------------
  ! Write atmospheric fields to a GRIB-file
  ! Generic routine to write a single atmospheric state
  !                       or an ensemble
  !----------------------------------------------------
  subroutine write_atm_gen (atm, file, pio, ie, mode, grid, ref, edition)
  type (t_atm)     ,intent(inout)        :: atm (:) ! atmosph. state to write
  character(len=*) ,intent(in)           :: file    ! file to write to
  integer          ,intent(in)           :: pio (:) ! processors to use for I/O
  integer          ,intent(in)           :: ie      ! ensemble index or -1
  character        ,intent(in) ,optional :: mode    ! 'w' or 'a' (append)
  character(len=*) ,intent(in) ,optional :: grid    ! grid-parameters to write
  logical          ,intent(in) ,optional :: ref     ! write reference fields
  integer          ,intent(in) ,optional :: edition ! GRIB edition to write

    type (t_grib1)         :: grib   ! GRIB buffers
    character              :: mod    ! local copy of 'mode'
    integer                :: i      ! loop index
    logical                :: lref   ! copy of 'ref'

    lref = .true.; if (present (ref)) lref = ref
    !----------------------------------------
    ! check for mode, must be write or append
    !----------------------------------------
    mod = 'w'; if (present(mode)) mod = mode
    select case (mod)
    case default
      call finish ('write_atm_gen','invalid mode: '//mod)
    case ('w','a')
    end select
    !---------------------------------------------------
    ! open file, allocate buffers, set verification time
    !---------------------------------------------------
    if (ie > 0) then
      call open_gribfile  (grib, file, mod, edition)
      call set_ref_time   (grib, atm(ie)% ref_time)
      call set_ver_time   (grib, atm(ie)% time)
      if      (local_dwd   (grib)) then
!!!     if (atm(ie)% member >= 0) then
        call set_dwd_ens (grib, atm(ie)% members, atm(ie)% member,&
                                atm(ie)% ensemble_id              )
!!!     end if
      else if (local_ecmwf (grib)) then
        !--------------------------------
        ! Always fill in ECMWF local part
        !--------------------------------
        call set_ecmwf_ens (grib, atm(ie)% members, atm(ie)% member)
      end if
    endif
    !-----------------------------------------------
    ! Set GRIB1/2 codes for all fields to be written
    !-----------------------------------------------
    if (ie > 0) then
       call set_gribcodes (atm(ie), force=.false., verbose=.false.)
!      call set_gribcodes (atm(ie), force=.false., verbose=.true.)
    end if
    !-------------------------------
    ! loop over the variables, write
    !-------------------------------
    do i=1,size (atm(1)% m)
!     !--------------------------------------------------------------
!     ! time range indicator = 'analysis' for certain constant fields
!     !--------------------------------------------------------------
!     if (ie > 0) then
!       if (atm(ie)% ref_time == atm(ie)% time) then
!         select case (atm(ie)% m(i)%i %name)
!         case ('lai','plcov','rootdp','vio3','hmo3')
!           call set_ver_time (grib, atm(ie)% time, range=WMO5_FORECAST)
!         case default
!           call set_ver_time   (grib, atm(ie)% time)
!         end select
!       endif
!     endif
      !------------------
      ! write GRIB record
      !------------------
      if (lref .or. .not. atm(1)% m(i)% i% ref)                     &
        call write_var_gen (grib, atm% m(i), pio, ie, atm(1)% grid, &
                            edition=edition                         )
    end do
    !-------------------------------
    ! optionally write the grid info
    !-------------------------------
    if (present (grid)) then
      !----------------------------------------------------------------
      ! Write invariant data with analysis date (hard-coded) for COSMO,
      ! but use 'classical date' 010101 for other models.
      !----------------------------------------------------------------
      if (model == 'COSMO') then
         call write_grib (atm(1)% grid, file, edition=edition, &
                          time=atm(1)% time, fields=grid,      &
                          grib=grib, pio=pio, ie=ie            )
      else
         call write_grib (atm(1)% grid, file, edition=edition, &
                          fields=grid,                         &
                          grib=grib, pio=pio, ie=ie            )
      end if
    endif
    !--------------------------------
    ! close the file, release buffers
    !--------------------------------
    if (ie > 0) then
      call close_gribfile (grib)
    endif
  end subroutine write_atm_gen
!------------------------------------------------------------------------------
  subroutine write_grid_grib (grid, file, mode, time, grib, pio, ie, edition, &
                              fields                                          )
  type (t_grid)    ,intent(in)             :: grid    ! atmospheric grid
  character(len=*) ,intent(in)             :: file    ! file to write to
  character        ,intent(in)   ,optional :: mode    ! 'w' or 'a' (append)
  type (t_time)    ,intent(in)   ,optional :: time    ! time to write to GRIB
  type (t_grib1)   ,intent(inout),optional :: grib    ! GRIB file descriptor
  integer          ,intent(in)   ,optional :: pio(:)  ! processors for I/O
  integer          ,intent(in)   ,optional :: ie      ! element to write
  integer          ,intent(in)   ,optional :: edition ! GRIB edition to write
  character(len=*) ,intent(in)   ,optional :: fields  ! parameters to write
  !-------------------------------------------------------------
  ! Write the grid to a GRIB-file.
  ! If 'grib' is present, it is used as the GRIB file descriptor
  ! and the GRIB-file is not opened and closed by this routine.
  !-------------------------------------------------------------

    target                  :: grib
    type (t_grib1) ,pointer :: gribp   ! local pointer to GRIB buffers
    character               :: mod     ! local copy of 'mode'
    integer                 :: i       ! loop index
    type(t_time)            :: tconst  ! time for constant fields
    logical                 :: lw      ! write on this PE
    type(t_m)  ,allocatable :: m(:)
    logical                 :: wlev    ! write levels
    character(len=256)      :: lf      ! local copy of 'fields'
    character(len=16)       :: flds (size(grid% m))
    integer                 :: lt

    !----------------------------------------
    ! check for mode, must be write or append
    !----------------------------------------
    mod = 'w'; if (present(mode)) mod = mode
    select case (mod)
    case default
      call finish ('write_grid_grib','invalid mode: '//mod)
    case ('w','a')
    end select

    !-------------------------------------------
    ! check for fields to be written in "fields"
    ! default: all fields despite hhl
    !-------------------------------------------
    if (present (fields)) then
      if (fields == '') return
    endif
    call   concat      (lf, grid% m% i% name)
    call   eval_string (lf,'- hhl rlon rlat geo_s geoid xnglob p0 dp0 rho0')
    if (present (fields)) then
      call eval_string (lf, fields)
    endif
    call split (flds, lf, lt)
    if (lt<0) call finish('write_grid_grib','increase flds')

    lw = dace% lpio
    if (present (ie)) then
      lw = ie>0
      allocate (m(size(pio)))
    endif

    if (lw) then
      !--------------------------------------------------
      ! either use GRIB file derived type variable passed
      ! or     set up local variable, open file
      !--------------------------------------------------
      if (present (grib)) then
        gribp => grib
      else
        allocate (gribp)
        !----------------------------
        ! open file, allocate buffers
        !----------------------------
        call open_gribfile (gribp, file, mod, edition)
      endif

      !----------------------
      ! set verification time
      !----------------------
      if (present (time)) then
        tconst = time
      else
        call init_time (tconst,1,1,1)
      endif
      call set_ref_time   (gribp, tconst)
      call set_ver_time   (gribp, tconst)
    else
      allocate (gribp)
    endif
    !-------------------------------
    ! loop over the variables, write
    !-------------------------------
    do i=1,size (grid% m)
      select case (grid% m(i)%i% name)
      case default
        !----------------------------------------------------------------------
        ! 'lsm','geosp','soiltyp','fr_lake','depth_lk','hsurf','hhl','sso_*'
        !----------------------------------------------------------------------
        wlev = (grid% m(i)%i% name == 'hhl')
        if (any (flds(1:lt) == grid% m(i)%i% name)) then
          if (present(ie)) then
            m = grid% m(i)
            call write_var_gen  (gribp,       m(:), pio, ie, grid, wlev=wlev, &
                                 edition=edition                              )
          else
            call write_var_grib (gribp, grid% m(i),          grid, wlev=wlev, &
                                 edition=edition                              )
          endif
        endif
      case ('rlon','rlat','geo_s','geoid','xnglob','p0','dp0','rho0')
        !---------------------------------------
        ! Internal use (not intended for export)
        !---------------------------------------
      end select
    end do
    !--------------------------------
    ! close the file, release buffers
    !--------------------------------
    if (lw) then
      if (.not. present(grib)) then
        call close_gribfile (gribp)
        deallocate (gribp)
      endif
    else
      deallocate (gribp)
    endif

  end subroutine write_grid_grib
!==============================================================================
  !---------------------------------
  ! write one field to the GRIB file
  !---------------------------------
  subroutine write_var_grib (grib, m, g, wlev, edition)
  type (t_grib1) ,intent(inout)        :: grib    ! GRIB record / file
  type (t_m)     ,intent(in)           :: m       ! field to write
  type (t_grid)  ,intent(in)           :: g       ! corresponding grid
  logical        ,intent(in) ,optional :: wlev    ! write levels to section 2
  integer        ,intent(in) ,optional :: edition ! GRIB edition to write

    !---------------------------------------------
    ! additional arguments passed to write_var_gen
    !---------------------------------------------
    type (t_m) :: lm  (1) ! local copy  (array of size 1 to fake ensemble)
    integer    :: pio (1) ! processors to use for I/O
    integer    :: ie      ! element to write on this processor element

    lm (1) = m
    pio(1) = dace% pio
    ie     = 1
    if (.not.dace% lpio) ie = -1

    call write_var_gen (grib, lm, pio, ie, g, wlev, edition)

  end subroutine write_var_grib
!------------------------------------------------------------------------------
  !-----------------------------------------------------------------
  ! write fields to the GRIB file
  ! generic routine to write a single field
  !                          or an ensemble of fields
  ! All ensemble members must have same grid, allocation status, etc
  !-----------------------------------------------------------------
  subroutine write_var_gen (grib, m, pio, ie, g, wlev, edition)
  type (t_grib1) ,intent(inout)        :: grib    ! GRIB record / file
  type (t_m)     ,intent(in)    TARGET :: m   (:) ! [ensemble of] fields
  integer        ,intent(in)           :: pio (:) ! processors to use for I/O
  integer        ,intent(in)           :: ie      ! member to write on this PE
  type (t_grid)  ,intent(in)           :: g       ! corresponding grid
  logical        ,intent(in) ,optional :: wlev    ! write levels to section 2
  integer        ,intent(in) ,optional :: edition ! GRIB edition to write

    type (ar_des)     :: des           ! GRIB descriptor
    integer           :: k             ! level index variable
    integer           :: code          ! grib code  number
    integer           :: table         ! grib table number
    integer           :: dis, cat, num ! GRIB2 discipline, category, number
    integer           :: ctyp          ! GRIB2 constituentType
    real(wp) _POINTER :: xg(:,:,:,:)   ! buffer to gather decomposed field
    integer           :: nz            ! number of vertical levels
    integer           :: levtyp
    integer           :: ffs, sfs      ! grib2 surface types
    logical           :: ldec
    real(wp)          :: lo1, la1
    real(wp)          :: lor, lar
    real(wp)          :: vcp  (128)
    real(wp)          :: lev           ! Level value
    integer           :: nvcp
    integer           :: zen           ! Zusatzelementnummer
    integer           :: ed            ! Edition
    integer           :: k1, k2        ! Level bounds
    integer           :: kk, kt, dk    ! Auxiliaries for level blocking
    integer(i8)       :: lev_size      ! Bytes / level
    integer, parameter:: realsize = size (transfer (0._wp, ['*'])) ! Bytes/wp

    !---------------------------------------------------
    ! Half-levels of soil model, hard-coded, for GRIB1
    ! The surface level is bounded by two half levels 0.
    !---------------------------------------------------
    real(wp) ,parameter :: soil_half_level_cm(0:9) = &
         (/ 0, 0, 1, 3, 9, 27, 81, 243, 729, 2187 /)
    !-----------------------------------------------------------
    ! Full-levels of soil model (rounded to nearest value in cm)
    !-----------------------------------------------------------
!   real(wp) ,parameter :: soil_full_level_cm(0:8) = &
!        nint ((soil_half_level_cm(0:8) + soil_half_level_cm(1:9)) / 2)
    !+++ Workaround for sxf90 stone age:
    real(wp)            :: soil_full_level_cm(0:8)
    !------------------------------------------------------------
    ! Half- and full-levels, converted to SI units [m], for GRIB2
    !------------------------------------------------------------
    real(wp) ,parameter :: soil_half_level(0:9) = soil_half_level_cm * 0.01_wp
    real(wp) ,parameter :: soil_full_level(0:8) = &
         (soil_half_level(0:8) + soil_half_level(1:9)) / 2

    !+++ Workaround for sxf90 stone age:
    soil_full_level_cm(0:8) = &
         nint ((soil_half_level_cm(0:8) + soil_half_level_cm(1:9)) / 2)

    nz = g% nz
    ed = grib% isec0% edition
    if (present (edition)) ed = edition

    if (all (m% i% alloc)) then
      !----------------------------------------
      ! process only gridpoint fields up to now
      !----------------------------------------
      select case (m(1)% i% rep)
      case ('gg')
        if (ie > 0) then
          !-----------------------------------------------------------------
          ! estimate grib code and table number, if not assigned to variable
          !-----------------------------------------------------------------
          code = m(ie)%i% code; table = m(ie)%i% table
          dis  = m(ie)%i% dis ; cat   = m(ie)%i% cat  ; num = m(ie)%i% num
          ctyp = m(ie)%i% ctyp
          if (code == 255) then
            des   = search (iname = m(ie)% i% name)
            code  = des% ee
            table = des% tabtyp
            if (code == 0) then
              write(0,*) m(ie)% i% name, grib% isec1% center,  &
                                         grib% isec1% sub_center
              call finish ('write_var_grib',                     &
                           'grib code not found: '//m(ie)% i% name)
            endif
          endif
          !---------------------------
          ! set lat-lon representation
          !---------------------------
          select case (g% gridtype)
          case (WMO6_LATLON, WMO6_ROTLL)
            !-----------------------------------------------
            ! shift in lat/lon for wind components on C-grid
            !-----------------------------------------------
            lo1 = g%lo1
            la1 = g%la1
            if      ((m(ie)%i%name == 'u').and.(g%arakawa == 'C')) then
              lo1 = g%lo1 + 0.5_wp*g%di
            else if ((m(ie)%i%name == 'v').and.(g%arakawa == 'C')) then
              la1 = g%la1 + 0.5_wp*g%dj
            endif
            !--------------------------------------------
            ! pass coordinates of S-pole for rotated grid
            !--------------------------------------------
            if (g% rot) then
              lar = - g% dlatr
              lor =   g% dlonr - 180._wp
              if (lor < -180._wp) lor = lor + 360._wp
              call set_grid_latlon (grib          ,&
                                    ni    = g% nx ,&
                                    nj    = g% ny ,&
                                    dlatf =    la1,&
                                    dlonf =    lo1,&
                                    di    = g% di ,&
                                    dj    = g% dj ,&
                                    dlatr =    lar,&
                                    dlonr =    lor )
            else
              !-----------------
              ! non-rotated grid
              !-----------------
#ifdef __ICON__
!DIR$ NOINLINE
#endif
              call set_grid_latlon (grib          ,&
                                    ni    = g% nx ,&
                                    nj    = g% ny ,&
                                    dlatf =    la1,&
                                    dlonf =    lo1,&
                                    di    = g% di ,&
                                    dj    = g% dj  )
            endif
          case (DWD6_ICOSAHEDRON)
            call set_grid_tri   (grib,         &
                                 ni    = g% ni,&
                                 nd    = g% nd )
          case (WMO6_GAUSSIAN)
            call set_grid_gauss (grib          ,&
                                 ni    = g% nx ,&
                                 nj    = g% ny ,&
                                 dlatf = g% la1,&
                                 dlonf = g% lo1,&
                                 di    = g% di, &
                                 scan  = WMO8_J_POSITIVE)    ! C.f. set_data
          case (DWD6_ICON)
            call set_grid_icon  (grib,                             &
                                 npts     = g% nxny,               &
                                 grid_num = g% icongrid% grid_num, &
                                 uuid     = g% icongrid% uuid      )
          case default
            call finish ('write_var_grib','gridtype not supported')
          end select
          !------------------------------------------------
          ! set leveltype in dependence of number of levels
          ! supply vct for hybrid coordinates
          !  1. determine level type of the variable
          !------------------------------------------------
          if (m(ie)% i% lb(3)==m(ie)% i% ub(3)) then
            select case (m(ie)% i% name)
            case default
              select case (g% levtyp)
              case default
                levtyp = WMO3_SURFACE
              case (WMO3_ISOBARIC)
                levtyp = g% levtyp      ! Special case: single pressure level?
              end select
            case ('t2m','rh2m','td2m','u_10m','v_10m', &
                  't2m_land', 'rh2m_land', 'td2m_land' )
              levtyp = WMO3_ABOVESUR
            case ('t_s')
              if (ed == 1) then
                levtyp = WMO3_BELOWSUR  ! GRIB1 convention (sic!)
              else
                levtyp = WMO3_SURFACE   ! GRIB2 convention
              end if
            end select
          else
            select case (g% levtyp)
            case (WMO3_HYBRID, WMO3_HYBRIDB, WMO3_HHYBRID, WMO3_GENV)
              if (m(ie)% i% ub(3) == g% nz) then
                levtyp = WMO3_HYBRIDB
              else if (m(ie)% i% ub(3) == g% nz + 1) then
                levtyp = WMO3_HYBRID
              else if (m(ie)% i% ub(3) == g% ns) then
                levtyp = WMO3_BELOWSUR
              else
                call finish('write_var_grib','invalid number of levels')
              endif
            case (WMO3_ISOBARIC)
              levtyp = g% levtyp
            case default
              if (m(ie)% i% ub(3) == g% nz .or. &
                  m(ie)% i% ub(3) == g% nz + 1  ) then
                print *, 'invalid level-type =', g% levtyp
                call finish ('write_var_grib',                             &
                             'invalid grid level-type for '//m(ie)% i% name)
              else if (m(ie)% i% ub(3) == g% ns) then
                levtyp = WMO3_BELOWSUR
              else
                call finish ('write_var_grib',                                &
                  'invalid number of levels / level-type for '//m(ie)% i% name)
              endif
            end select
          endif
          !-----------------------------------------------------------
          ! 2. call set_grid_vert in dependence of grid     level type
          !                                    and variable level type
          !-----------------------------------------------------------
          select case (g% levtyp)
          case (WMO3_HYBRID, WMO3_HYBRIDB, WMO3_HHYBRID, WMO3_GENV)
            select case (levtyp)
            case (WMO3_HYBRID, WMO3_HYBRIDB, WMO3_HHYBRID, WMO3_GENV)
              if (g% levtyp == WMO3_GENV) then
                !--------------------------------
                ! Generalized vertical coordinate
                !--------------------------------
                call set_grid_vert (grib,    WMO3_GENV,                &
                                             nlev    = nz+1,           &
                                             gridnum = g% vc% ivctype, &
                                             uuid    = g% vc% vc_uuid, &
                                             wlev    = .false.         )
              else if (g% vc% ivctype == IVCTYPE_GME) then
                !-----------------------------------
                ! GME/HRM hybrid pressure coordinate
                !-----------------------------------
                call set_grid_vert (grib,    levtyp,     &
                                          g% ak(1:nz+1), &
                                          g% bk(1:nz+1), &
                                             wlev=.true. )
              else if (g% vc% ivctype == IVCTYPE_ICON) then
                !------------------------------
                ! ICON height hybrid coordinate
                !------------------------------
                call set_grid_vert (grib,    levtyp,     &
                                          g% ak(1:nz+1), &
                                          g% bk(1:nz+1), &
                                             wlev=.true. )
              else
                !--------------------------
                ! COSMO hybrid z-coordinate
                !--------------------------
                call cosmo_vcp (vcp, nvcp, g)
                call set_grid_vert (grib,    levtyp,     &
                                             vcp(1:nvcp),&
                                             wlev=.true. )
              endif
            case default
              if (g% levtyp == WMO3_GENV) then
                !--------------------------------
                ! Generalized vertical coordinate
                !--------------------------------
                call set_grid_vert (grib,    levtyp,     &
                                             wlev=.false.)
              else if (g% vc% ivctype == IVCTYPE_GME) then
                !-----------------------------------
                ! GME/HRM hybrid pressure coordinate
                !-----------------------------------
                call set_grid_vert (grib,    levtyp,     &
                                          g% ak(1:nz+1), &
                                          g% bk(1:nz+1), &
                                             wlev=wlev )
              else if (g% vc% ivctype == IVCTYPE_ICON) then
                !------------------------------
                ! ICON height hybrid coordinate
                !------------------------------
                call set_grid_vert (grib,    levtyp,     &
                                          g% ak(1:nz+1), &
                                          g% bk(1:nz+1), &
                                             wlev=wlev )
              else
                !--------------------------
                ! COSMO hybrid z-coordinate
                !--------------------------
                call cosmo_vcp (vcp, nvcp, g)
                call set_grid_vert (grib,    levtyp,     &
                                             vcp(1:nvcp),&
                                             wlev=wlev )
              endif
            end select
          case (WMO3_ISOBARIC)
             call set_grid_vert (grib,    levtyp,    &
                                       g% akf(1:nz), &
                                          wlev=wlev  )
          case default
!            call finish ('write_var_grib','invalid grid level-type')
             call set_grid_vert (grib,    levtyp,     &
                                          wlev=.false.)
          end select
          !------------------------
          ! set grib code and table
          !------------------------
          zen = 0
          !++++++++++++++++++++++++++++++++++++++++++++
          ! set Zusatzelementnummer for GME W_SO to 255
          !++++++++++++++++++++++++++++++++++++++++++++
          if (g% gridtype == DWD6_ICOSAHEDRON .and. m(ie)% i% name == 'w_so') &
               zen = 255
          call set_code (grib, code, table, dis, cat, num, element_no=zen, &
                         ctyp=ctyp,  bits=m(ie)% i% bits)
        endif
        !------------------------------
        ! gather field on I/O processor
        !------------------------------
!       ldec = g% ldec
        ldec = m(1)%i% lb(1) /= g%lbg(1) .or. m(1)%i% lb(2) /= g%lbg(2) .or. &
               m(1)%i% ub(1) /= g%ubg(1) .or. m(1)%i% ub(2) /= g%ubg(2)
        k1   = m(1)%i% lb(3)
        k2   = m(1)%i% ub(3)
        !--------------------------------
        ! Set max.no. of levels per chunk
        !--------------------------------
        if ((ldec .or. g% d_gme(1) >= 0) .and. io_max_gather > 0) then
           lev_size = (g%ubg(1)-g%lbg(1)+1) * (g%ubg(2)-g%lbg(2)+1) * realsize
           dk = max (int (io_max_gather / lev_size), 1)
        else
           dk = k2 - k1 + 1
        end if
        !----------------------------
        ! Level blocking (outer loop)
        !----------------------------
        do kk = k1, k2, dk
          kt = min (kk+dk-1, k2)
          !-----------------
          ! Allocate buffers
          !-----------------
          nullify (xg)
          if (ldec .or. g% d_gme(1) >= 0) then
            if (ie > 0) then
              allocate (xg (g%lbg(1):g%ubg(1),  & ! Array slice on receiver
                            g%lbg(2):g%ubg(2),  &
                                  kk:kt      ,  & ! m(1)%i%lb(3):m(1)%i%ub(3)
                            g%lbg(4):g%ubg(4) ) )
            else
              allocate (xg (0, 0, kk:kt, 0))      ! Dummy slice on sender
            end if
          endif
          !-----------------------
          ! provide data and write
          !-----------------------
          if (ldec) then
!           !+++++++++++++++++++++++++++++++++++++++++++++++
!           ! intermediate use of gather_multi (for testing)
!           !+++++++++++++++++++++++++++++++++++++++++++++++
!           do i = 1, size (m)
!             call gather_multi (xg, m(i)% ptr (:,:,:,:), g% dc, pio(i))
!           end do
            !-----------------------------------------------
            ! gather fields using MPI alltoall communication
            !-----------------------------------------------
            call alltoallv_multi_w (xg, m, g% dc,  pio(:))
            if (ie > 0) call permut_write (g, x4=xg)
          else
            if (ie > 0) then
              if (g% d_gme(1) >= 0) then
                xg =  m(ie)% ptr (:,:,kk:kt,:)
                call permut_write (g, x4=xg)
              else
                xg => m(ie)% ptr
              endif
            endif
          endif
          !
          if (ie > 0) then
           !--------------------------
           ! do not write dummy fields
           !--------------------------
           if (m(ie)% i% dummy) then
             if (ldec .or. g% d_gme(1) >= 0) deallocate (xg)
             cycle
           end if
           !-----------------
           ! loop over levels
           !-----------------
!          do k = m(1)% i% lb(3), m(1)% i% ub(3)
           do k = kk, kt
           !--------------------------------------------
           ! set level in dependence of number of levels
           !--------------------------------------------
            !--------------------
            ! single level fields
            !--------------------
            if (m(ie)% i% lb(3) == m(ie)% i% ub(3)) then
              select case (m(ie)% i% name)
              case default
!print*,"### single level???", m(ie)% i% name, levtyp
                select case (levtyp)
                case default
                  lev = 0._wp
                case (WMO3_ISOBARIC)
                  lev = g% akf(k)       ! Special case: single pressure level?
                end select
                if (ed == 1) then
                  call set_level (grib,    lev, levtyp)
                else
                  ffs = m(ie)% i% ffs
                  sfs = m(ie)% i% sfs
                  !--------------------------------------------------------------
                  ! Catch potential multi-layer fields that are single-level here
                  !--------------------------------------------------------------
                  if (ffs == -1) ffs = levtyp
!print*,"### call set_level???", ffs, sfs, lev
                  call set_level (grib,    lev, levtyp=ffs, levtyp2=sfs)
                end if
              case ('hsurf')
                call set_level   (grib,  0._wp, WMO3_SURFACE, &
                                  levtyp2   =   WMO3_SEALEVEL )
              case ('pmsl')
                call set_level   (grib,  0._wp, WMO3_SEALEVEL)
              case ('t2m','rh2m','td2m')
                call set_level   (grib,  2._wp, WMO3_ABOVESUR)
              case ('t2m_land','rh2m_land','td2m_land')
                call set_level   (grib,  2._wp, WMO3_ABOVESUR, &
                                  level2    =   0._wp,         &
                                  levtyp2   =   WMO3_TILE_LAND )
              case ('u_10m','v_10m')
                call set_level   (grib, 10._wp, WMO3_ABOVESUR)
              case ('t_so')
                !--------------------------------------
                ! special case for "single level" t_so(0)
                ! +++ temporarily hard coded levels +++
                ! must distinguish grib1 from grib2.
                !--------------------------------------
                if (ed == 1) then
                  call set_level (grib, soil_full_level_cm(k)*0.01_wp,      &
                                                               WMO3_BELOWSUR)
                else
                  call set_level (grib, soil_full_level(k),    WMO3_BELOWSUR)
                end if
              end select
            else
              select case (m(ie)% i% name)
              case ('t_so','w_so','w_so_ice','smi')
                !--------------------------------------
                ! below surface fields
                ! +++ temporarily hard coded levels +++
                ! must distinguish grib1 from grib2.
                !--------------------------------------
                if (ed == 1) then
                  call set_level   (grib, soil_full_level_cm(k)*0.01_wp, &
                                                                WMO3_BELOWSUR)
                else
                  if (m(ie)% i% name == 't_so') then
                    call set_level (grib, soil_full_level(k),   WMO3_BELOWSUR)
                  else
                    call set_level (grib, soil_half_level(k),   WMO3_BELOWSUR, &
                                          soil_half_level(k+1), WMO3_BELOWSUR)
                  end if
                end if
              case default
                !-------------------
                ! atmospheric fields
                !-------------------
                select case (g% levtyp)
                case (WMO3_ISOBARIC)
                  call set_level (grib, g% akf(k), WMO3_ISOBARIC)
                case default
                  if (m(ie)% i% name == 'hhl') then
                    call set_level (grib, real (k,wp), levtyp, &
                                    levtyp2 = WMO3_SEALEVEL)
                  else
                    call set_level (grib, real (k,wp), levtyp)
                  end if
                end select
              end select
            endif
            !------
            ! write
            !------
            call set_data (grib, xg(:,:,k,:))     ! assumes orientation S->N
            if (m(ie)% i% name == 'pp') &         ! for pressure deviation ..
              grib% rsec4 = grib% rsec4 * 0.01_wp ! .. convert Pa -> hPa
            if (m(ie)% i% name == 'rh2m') &       ! for 2m rel.hum.
              grib% rsec4 = grib% rsec4 * 100._wp ! .. convert 1 -> %
            call write_gribrecord (grib, edition)
           end do ! levels
          endif   ! ie > 0
          if (ldec .or. g% d_gme(1) >= 0) deallocate (xg)
        end do    ! kk loop (level blocking)
      case default
        call finish ('write_var_grib','representation not supported')
      end select
    endif
  end subroutine write_var_gen

!==============================================================================
  subroutine set_gribcodes (state, force, verbose)
    type(t_atm),  target, intent(inout) :: state    ! Atmospheric state
    logical,    optional, intent(in)    :: force    ! Forcing consistent codes
    logical,    optional, intent(in)    :: verbose  ! Verbose conversion

    logical               :: lforce  ! local copy of force
    logical               :: vrb     ! local copy of verbose
    integer               :: i       ! loop index
    type(t_grid), pointer :: g       ! grid pointer
    type(t_mi),   pointer :: mi      ! metadata pointer

    lforce = .false.; if (present (force))   lforce = force
    vrb    = .false.; if (present (verbose)) vrb    = verbose

    g => state% grid
    do i = 1, size (state% m)
       mi => state% m(i)% i
       if (.not. mi% alloc) cycle
       call set_codes (mi)
    end do
    do i = 1, size (g% m)
       mi => g% m(i)% i
       if (.not. mi% alloc) cycle
       select case (mi% name)
       case default
          !----------------------------------------------------------------------
          ! 'lsm','geosp','soiltyp','fr_lake','depth_lk','hsurf','hhl','sso_*'
          !----------------------------------------------------------------------
          call set_codes (mi)
       case('rlon','rlat')
          !---------------------------------------------------------------
          ! Internal use (not intended for export), GRIB encoding possible
          !---------------------------------------------------------------
       case('geo_s','geoid','xnglob','p0','dp0','rho0')
          !--------------------------------------------------
          ! Internal use only, not intended for GRIB encoding
          !--------------------------------------------------
       end select
    end do
  contains
    subroutine set_codes (mi)
      type(t_mi),   pointer :: mi      ! metadata pointer

      logical               :: sl      ! single level field
      integer               :: ed      ! source grib edition
      integer               :: levtyp  ! level type
      integer               :: level   ! level value
      integer               :: T_FFS   ! typeOfFirstFixedSurface
      integer               :: T_SFS   ! typeOfSecondFixedSurface
      type(ar_des)          :: des     ! GRIB descriptor (old)
      type(t_par_grib1)     :: code1   ! GRIB1 parameters
      type(t_par_grib2)     :: code2   ! GRIB2 parameters
      character(len=10)     :: iname   ! internally used name
      logical               :: ok1     ! GRIB 1 codes are set
      logical               :: ok2     ! GRIB 2 codes are set
      !---------------------------------------------
      ! Specify level type / typeOfFirstFixedSurface
      ! for disambiguation of table lookup.
      !---------------------------------------------
      levtyp = -1
      T_FFS  = -1
      level  = -1
      T_SFS  = empty2% T_SFS
      sl     = (mi% lb(3) == mi% ub(3))
      if (sl) then
         select case (mi% name)
         case default
            select case (g% levtyp)
            case default
               levtyp = WMO3_SURFACE
               T_FFS  = 1
               level  = 0
            case (WMO3_ISOBARIC)
               ! Defer for single pressure levels
            end select
         case ('t2m','rh2m','td2m')
            levtyp = WMO3_ABOVESUR
            T_FFS  = 103
            level  = 2
         case ('t2m_land','rh2m_land','td2m_land')
            levtyp = WMO3_ABOVESUR
            T_FFS  = 103
            level  = 2
            T_SFS  = WMO3_TILE_LAND  ! GRIB2 only: land-tile average
         case ('u_10m','v_10m')
            levtyp = WMO3_ABOVESUR
            T_FFS  = 103
            level  = 10
         case ('t_s')
            levtyp = WMO3_BELOWSUR
            T_FFS  = 1          ! GRIB1: below surface, GRIB2: surface!
            level  = 0
         case ('pmsl')
            levtyp = WMO3_SEALEVEL
            T_FFS  = 101
            level  = 0
         end select
      else
        select case (mi% name)
        case ('t_so', 'w_so', 'w_so_ice', 'smi')  ! exclude soil fields
        case default
          select case (g% levtyp)
          case (WMO3_HYBRID, WMO3_HYBRIDB, WMO3_HHYBRID, WMO3_GENV)
            if (mi% ub(3) == g% nz) then
               select case (mi% name)
               case ('geof')               ! Distinguish geof(FIF)/geoh(FI)
                  levtyp = WMO3_HYBRIDB
                  T_FFS  = 105
               end select
            else if (mi% ub(3) == g% nz + 1) then
               levtyp = WMO3_HYBRID
            else if (mi% ub(3) == g% ns) then
!              levtyp = WMO3_BELOWSUR
            else
               call finish('set_gribcodes','invalid number of levels')
            endif
          case (WMO3_ISOBARIC)
!           levtyp = WMO3_ISOBARIC
!           T_FFS  = WMO3_ISOBARIC
          case default
!           print *, 'invalid level-type =', g% levtyp
            call finish ('set_gribcodes','unsupported grid level-type for ' &
                         //trim (mi% name)//' : '//char3 (g% levtyp)        )
          end select
        end select
      end if
      !--------------------------------------------------
      ! Guess "edition" where grib codes were already set
      ! (This may need to be revised for non-DWD gribs).
      !--------------------------------------------------
      ed = -1
      ok1 = all ( (/mi% code ,mi% table          /) /= 255)
      ok2 = all ( (/mi% dis  ,mi% cat   ,mi% num /) /= 255)
      if (ok1) ed = 1
      if (ok2) ed = 2
      !----------------------------------------
      ! Delete old ECMWF code/table assignments
      ! we don't handle mixed grib1/grib2 cruft
      !----------------------------------------
      if (mi% table == 128) then
         mi% code  = 255
         mi% table = 255
         ed        = -1
         ok1       = .false.
      end if
      !------------------------
      ! Last resort: guess code
      !------------------------
      if (ed < 0) then
         iname = mi% name
         if (iname == "psr") iname = "ps"
         des = search (iname = iname)
         if (des% ee == 0) &
              call finish ('set_gribcodes','grib code not found: '//mi% name)
         mi% code  = des% ee
         mi% table = des% tabtyp
         ed = 1
         if (vrb) &
              print *, "set_gribcodes: fallback table,code=", mi%table,mi%code
      end if
      if (ed == 1) then
         !---------------------------
         ! Check for full match first
         !---------------------------
         if (levtyp == WMO3_ABOVESUR .or. levtyp == WMO3_BELOWSUR) then
            code1 = search_grib1 (code=mi% code, table=mi% table, levtyp=levtyp,&
                                                                  level =level  )
         else
            code1 = search_grib1 (code=mi% code, table=mi% table, levtyp=levtyp)
         end if
         !------------------------
         ! Last resort: any levtyp
         !------------------------
         if (code1% paramID < 0) then
            code1 = search_grib1 (code=mi% code, table=mi% table)
         end if
         if (vrb .or. code1% paramID < 0) then
            print '(1x,a,a,99(1x,a,i0))', "grib1 <- ", mi% name, &
                 "table=", mi% table, "code=", mi% code, &
                 "paramID=", code1% paramID, "levtyp=", code1% levtyp
         end if
         if (code1% paramID < 0) call message ("set_gribcodes",&
              "parameter not in GRIB1 table:"//trim (mi% name))
         if (lforce .or. .not. ok2) then
            code2 = search_grib2 (paramID = code1% paramID)
            if (code2% paramID < 0) then
               write(0,*) "paramID=", code1% paramID, &
                    "shortname=",trim (mi% name)
               if (lforce) then
                 call finish  ("set_gribcodes",                       &
                      "parameter not in GRIB2 table:"//trim (mi% name))
               else
                 call message ("set_gribcodes",                       &
                      "parameter not in GRIB2 table:"//trim (mi% name))
               endif
            end if
            if (vrb) then
               print '(1x,a,a,99(1x,a,i0))', "grib2 -> ",  code2% shortname, &
                    "dis=", code2% discipline, "cat=", code2% category, &
                    "num=", code2% number, "T_FFS=", code2% T_FFS
            end if
            if (code2% paramID > 0) then
               mi% dis = code2% discipline
               mi% cat = code2% category
               mi% num = code2% number
               mi% ffs = code2% T_FFS
               mi% sfs = code2% T_SFS
            end if
         end if
         mi% parID = code1% paramID
      else ! edition /= 1
         !------------------------------------------------
         ! Validate entry in DWD table of "local concepts"
         !------------------------------------------------
         code2 = search_grib2 (dis    =mi% dis, &
                               cat    =mi% cat, &
                               num    =mi% num, &
                               levtyp1=T_FFS,   &
                               levtyp2=T_SFS)
         !------------------------
         ! Last resort: any levtyp
         !------------------------
         if (code2% paramID < 0) then
            code2 = search_grib2 (dis=mi% dis, &
                                  cat=mi% cat, &
                                  num=mi% num  )
         end if
         if (vrb .or. code2% paramID < 0) then
            print '(1x,a,a,99(1x,a,i0))', "grib2 <- ", mi% name, &
                 "dis=", mi% dis, "cat=", mi% cat, "num=", mi% num, &
                 "paramID=", code2% paramID, "T_FFS=", code2% T_FFS
         end if
         if (code2% paramID < 0) call message ("set_gribcodes",&
              "parameter not in GRIB2 table:"//trim (mi% name))
         if (lforce .or. .not. ok1) then
            code1 = search_grib1 (paramID = code2% paramID)
            if (code1% paramID < 0) then
               call message ("set_gribcodes",&
                    "parameter not in GRIB1 table:"//trim (mi% name))
            else
               mi% code  = code1% code
               mi% table = code1% table
            end if
            if (vrb) then
               print '(1x,a,a,99(1x,a,i0))', "grib1 -> ",  code1% shortname, &
                    "table=", code1% table, "code=", code1% code, &
                    "levtyp=", code1% levtyp
            end if
         end if
         mi% parID = code2% paramID
         mi% ffs   = code2% T_FFS
         mi% sfs   = code2% T_SFS
      end if
      !----------------------------------------------------------------------
      ! set constituentType for Category == Atmospheric chemical constituents
      !----------------------------------------------------------------------
      if (mi% cat == 20) then
        code2 = search_grib2 (name=toupper(mi% name), &
                              dis =        mi% dis,   &
                              cat =        mi% cat,   &
                              num =        mi% num    )
        mi% ctyp = code2% cTyp
      endif
    end subroutine set_codes
  end subroutine set_gribcodes
!==============================================================================
  elemental subroutine set_bits (state, fields, bits)
    type(t_atm), target, intent(inout) :: state   ! atmospheric state to write
    character(len=*),    intent(in)    :: fields  ! names of fields
    integer,             intent(in)    :: bits    ! encoding precision
    !-----------------------------------------
    ! Set grib bits/value for specified fields
    !-----------------------------------------
    character(len(fields)+2) :: names
    integer                  :: i, k
    type(t_mi),      pointer :: mi      ! metadata pointer

    if (fields == "") return
    names = " " // trim (fields) // " "

    do i = 1, size (state% m)
       mi => state% m(i)% i
       if (.not. mi% alloc) cycle
       k = index (names, " "//trim (mi% name)//" ")
       if (k > 0) mi% bits = bits
    end do

    do i = 1, size (state% grid% m)
       mi => state% grid% m(i)% i
       if (.not. mi% alloc) cycle
       if (      mi% dummy) cycle
       k = index (names, " "//trim (mi% name)//" ")
       if (k > 0) mi% bits = bits
    end do
  end subroutine set_bits
!==============================================================================
  subroutine permut_read (grid, x, x4)
  !-----------------------------------------------
  ! permutate global fields after reading
  ! used for re-ordering of ICON ensemble in LETKF
  !-----------------------------------------------
  type(t_grid) ,intent(in)              :: grid
  real(wp)     ,intent(inout) ,optional :: x  (:,:,:)
  real(wp)     ,intent(inout) ,optional :: x4 (:,:,:,:)

    if (grid% d_gme(1) >= 0 .and. associated (grid% marr)) then
      if (present (x )) x ( grid% marr(2,:,1,1), 1,    1) = x  (:,1,  1)
      if (present (x4)) x4( grid% marr(2,:,1,1), 1, 1, 1) = x4 (:,1,1,1)
    endif

  end subroutine permut_read
!==============================================================================
  subroutine permut_write (grid, x, x4)
  !-----------------------------------------------
  ! permutate global fields before writing
  ! used for re-ordering of ICON ensemble in LETKF
  !-----------------------------------------------
  type(t_grid) ,intent(in)              :: grid
  real(wp)     ,intent(inout) ,optional :: x  (:,:,:)
  real(wp)     ,intent(inout) ,optional :: x4 (:,:,:,:)

#if HAVE_F2008_CONTIGUOUS
  CONTIGUOUS :: x, x4
#endif

    if (grid% d_gme(1) >= 0 .and. associated (grid% marr)) then
      if (present (x )) x  (:,1,  1) = x ( grid% marr(2,:,1,1), 1,    1)
      if (present (x4)) then
!$omp parallel workshare
                        x4 (:,1,:,1) = x4( grid% marr(2,:,1,1), 1, :, 1)
!$omp end parallel workshare
      end if
    endif

  end subroutine permut_write
!==============================================================================
  subroutine construct_ctr_slot (ctr, nslots, nens, det)
    type(t_ctr_slot)   ,intent(out) ,allocatable :: ctr(:)   ! container
    integer            ,intent(in)  ,optional    :: nslots   ! time slots
    integer            ,intent(in)  ,optional    :: nens     ! ensemble size
    logical            ,intent(in)  ,optional    :: det      ! determinsitic?

    integer                              :: n_slots     ! local copy of nslots
    integer                              :: n_ens       ! local copy of nens
    logical                              :: ldet        ! local copy of det
    integer                              :: i_det       ! index of deterministic
    integer                              :: size        ! n_ens + det
    integer                              :: i           ! loop index

    n_slots = 1;       if (present(nslots)) n_slots = nslots
    n_ens   = 0;       if (present(nens))   n_ens   = nens
    ldet    = .false.; if (present(det))    ldet    = det

    i_det = -1;    if (ldet) i_det = n_ens + 1
    size  = n_ens; if (ldet) size  = n_ens + 1

    ! construct outer slot container
    allocate (ctr(n_slots))
    ctr(:)% slot   = (/ (i, i = 1, n_slots) /)
    ctr(:)% nslots = n_slots
    ctr(:)% nens   = n_ens

    ! construct inner ensemble container
    do i = 1, n_slots
      allocate (ctr(i)% c(size))
      if (ldet) ctr(i)% c(i_det)% det = .true.
    end do

  end subroutine construct_ctr_slot
!------------------------------------------------------------------------------
  subroutine destruct_ctr_slot (ctr)
    type(t_ctr_slot)    ,intent(inout) ,allocatable :: ctr(:)  ! container

    integer    :: ns
    integer    :: i, j

    if (allocated(ctr)) then
      ns = size(ctr)
      do i = 1, ns
        do j = i + 1, ns
          if (associated (ctr(i)% c, ctr(j)% c)) nullify(ctr(j)% c)
        end do
      end do

      do i = 1, ns
        if (associated(ctr(i)% c)) then
          do j = 1, size(ctr(i)% c)
            if (associated (ctr(i)% c(j)% invt)) deallocate (ctr(i)% c(j)% invt)
            nullify (ctr(i)% c(j)% invt)
          end do
          deallocate (ctr(i)% c)
        end if
      end do
      deallocate (ctr)
    end if

  end subroutine destruct_ctr_slot
!------------------------------------------------------------------------------
  subroutine read_inv_ctr_slot (ctr, file, nens, suffix, members, ifile)
  !-------------------------------------------------
  ! read inventories of ensemble members in parallel
  !-------------------------------------------------
  type (t_ctr_slot)            ,intent(inout)  :: ctr          ! inventory container
  character (len=*)            ,intent(in)     :: file         ! GRIB file name
  integer           ,optional  ,intent(in)     :: nens         ! ensemble size
  character (len=*) ,optional  ,intent(in)     :: suffix       ! suffix
  integer           ,optional  ,intent(in)     :: members(:)   ! ensemble members to read
  integer           ,optional  ,intent(in)     :: ifile        ! file index

    call read_inv_ctr (ctr% c, file, nens, suffix, members, ifile)

  end subroutine read_inv_ctr_slot
!------------------------------------------------------------------------------
  subroutine read_inv_ctr (ctr, file, nens, suffix, members, ifile)
  !---------------------------------------------------
  ! read inventories of ensemble members in parallel
  ! inventories are keeped on the reading pe's
  ! deterministic inventory is broadcasted to all pe's
  !---------------------------------------------------
  type (t_ctr_inv)             ,intent(inout)  :: ctr(:)       ! inventory container
  character (len=*)            ,intent(in)     :: file         ! GRIB file name
  integer           ,optional  ,intent(in)     :: nens         ! ensemble size
  character (len=*) ,optional  ,intent(in)     :: suffix       ! suffix
  integer           ,optional  ,intent(in)     :: members(:)   ! ensemble members to read
  integer           ,optional  ,intent(in)     :: ifile        ! file index

    !----------------
    ! local variables
    !----------------
    character (len=256) ,allocatable  :: fname(:)      ! grib file name
    integer             ,allocatable  :: pio(:)        ! PEs to use for reading
    integer                           :: nchunk        ! number of chunks to read
    integer                           :: mcsize        ! maximum chunk size
    integer                           :: i             ! start index of chunk to read
    integer                           :: nm            ! size        of chunk to read
    integer                           :: ie            ! ensemble member index in chunk
    integer                           :: ke            ! ensemble member index

    !-----------------------
    ! read deterministic run
    !-----------------------
    if (.not.present(nens)) then
      allocate (fname(1))
      fname(1) = path_file ('',file,suffix=suffix)
      do i = 1, size(ctr)
        if (ctr(i)% det) then
          ke = i
          exit
        end if
      end do
      call get_inventory (ctr(ke)% invt, fname(1), ifile=ifile)
      ctr(ke)% name   = fname(1)
      ctr(ke)% global = .true.
      ctr(ke)% pe     = dace% pio
      ctr(ke)% size   = size(ctr(ke)% invt)
      return
    end if

    !----------------------
    ! read ensemble members
    !----------------------
    if (present(nens)) then
      !-----------------------------------------------
      ! Set I/O processor pattern for ensemble members
      !-----------------------------------------------
      allocate (pio(nens))
      call set_iope_ens (pio, nchunk, mcsize, stride=-1)

      allocate (fname(nens))
      do i = 1, nens
        ke = i; if (present(members)) ke = members(i)
        fname(i) = path_file ('',file, iens=ke, suffix=suffix)
      end do

      !-----------------------------------------------------------------
      ! distribute read of inventories, 1 member per chosen PE at a time
      !-----------------------------------------------------------------
      if (dace% lpio) write (6,'(a,/)') repeat('-',79)
      do i = 1, nens, mcsize              ! loop over bunches to read
        nm = min (nens-i+1, mcsize)       ! number of members to read
        if (dace% lpio) &
          write(6,'(a,i3,a,i3,/)') ' reading inventory of members',i,' to',i+nm-1
        ie = -1
        ke = -1
        if (any (pio(i:i+nm-1) == dace% pe)) then
          ie = minloc (abs (pio(i:i+nm-1) - dace% pe),1)  ! index within chunk
          ke = i - 1 + ie                                 ! member to read on this PE
          !print *, "### get_inventory_read: p_pe,ie,ke =", dace% pe, ie, ke
          call get_inventory (ctr(ke)% invt, fname(ke), pio=pio(ke),  &
                              comm=self% comm, ifile=ifile            )
        end if
      end do

      !----------------------------
      ! broadcast result to all PEs
      !----------------------------
      do i = 1, nens
!        call p_bcast (ctr(i)% invt, pio(i))
        if (.not.associated(ctr(i)% invt)) allocate (ctr(i)% invt(0))
        ctr(i)% name    = fname(i)
        ctr(i)% global = .false.
        ctr(i)% pe      = pio(i)
        ctr(i)% size    = size(ctr(i)% invt)
      end do
    end if

  end subroutine read_inv_ctr
!------------------------------------------------------------------------------
  subroutine read_inv_files (ctr, files, ifile)
  !-----------------------------------------
  ! read inventories of files given by files
  ! inventories are broadcasted to all pe's
  !-----------------------------------------
  type (t_ctr_inv)             ,intent(inout)  :: ctr(:)       ! inventory container
  character (len=*)            ,intent(in)     :: files(:)     ! GRIB file names
  integer           ,optional  ,intent(in)     :: ifile        ! file index

    !----------------
    ! local variables
    !----------------
    integer                           :: nfiles        ! number of files
    integer             ,allocatable  :: pio(:)        ! PEs to use for reading
    integer                           :: nchunk        ! number of chunks to read
    integer                           :: mcsize        ! maximum chunk size
    integer                           :: i             ! start index of chunk to read
    integer                           :: nm            ! size        of chunk to read
    integer                           :: ie            ! ensemble member index in chunk
    integer                           :: ke            ! ensemble member index
    logical                           :: ex            ! file existence

    if (size(ctr) /= size(files)) call finish ('read_inv_files','size(ctr) /= size(files)')
    nfiles = size(files)

    !-------------------------
    ! check existence of files
    !-------------------------
    if (dace% lpio) then
      do i = 1, nfiles
        inquire (file=files(i), exist=ex)
        if (.not.ex) call finish ('read_inv_files','missing file: '//trim(files(i)))
      end do
    end if

    !--------------------------------------------
    ! Set I/O processor pattern for reading files
    !--------------------------------------------
    allocate (pio(nfiles))
    call set_iope_ens (pio, nchunk, mcsize, stride=-1)

    !---------------------------------------------------------------
    ! distribute read of inventories, 1 file per chosen PE at a time
    !---------------------------------------------------------------
    if (dace% lpio) write (6,'(a,/)') repeat('-',79)
    do i = 1, nfiles, mcsize              ! loop over bunches to read
      nm = min (nfiles-i+1, mcsize)       ! number of files to read
      if (dace% lpio) &
        write(6,'(a,i3,a,i3,/)') ' reading inventory of files',i,' to',i+nm-1
      ie = -1
      ke = -1
      if (any (pio(i:i+nm-1) == dace% pe)) then
        ie = minloc (abs (pio(i:i+nm-1) - dace% pe),1)  ! index within chunk
        ke = i - 1 + ie                                 ! file to read on this PE
        !print *, "### get_inventory_read: p_pe,ie,ke =", dace% pe, ie, ke
        call get_inventory (ctr(ke)% invt, files(ke), pio=pio(ke),  &
                            comm=self% comm, ifile=ifile            )
      end if
    end do

    !----------------------------
    ! broadcast result to all PEs
    !----------------------------
    do i = 1, size(files)
      call p_bcast (ctr(i)% invt, pio(i))
      if (.not.associated(ctr(i)% invt)) allocate (ctr(i)% invt(0))
      ctr(i)% name   = files(i)
      ctr(i)% global = .true.
      ctr(i)% pe     = pio(i)
      ctr(i)% size   = size(ctr(i)% invt)
    end do

  end subroutine read_inv_files
!==============================================================================
end module mo_grib
