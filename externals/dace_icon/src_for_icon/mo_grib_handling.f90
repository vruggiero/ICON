!
!+ GRIB handling library:  get inventory,  set PDS/GDS entries, ...
!
MODULE mo_grib_handling
!
! Description:
!   GRIB handling library.
!   Facilities:
!   - get inventory of a GRIB-File
!   - get verification or reference time from a GRIB record
!   - set PDS or GDS entries before writing a GRIB record
!   - compare (components of) a GRIB record
!   - write GRIB file
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
! V1_4         2009/03/26 Andreas Rhodin
!  changes for COSMO LETKF
! V1_5         2009/05/25 Andreas Rhodin
!  Extend GRIB file inventory
! V1_6         2009/06/10 Andreas Rhodin
!  interprete GRIB PDS time_range=13 (DWD special value for nudging)
! V1_7         2009/08/24 Andreas Rhodin
!  Changes for COSMO/LETKF; fix century in DWD database-entry date
! V1_8         2009/12/09 Andreas Rhodin
!  get_inventory: return zero size inventory if file is not present
!  use DWD5_IFS_FC = 14 (DWD specific time range indicator for IFS fields)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  get_inventory: correctly derive no.points (for reduced grids)
! V1_13        2011/11/01 Andreas Rhodin
!  changes for ensembles, GRIB2 API
! V1_19        2012-04-16 Andreas Rhodin
!  print_inventory: write zen (Zusatzelementnummer)
!  set Zusatzelementnummer for GME W_SO to 255 for GME LETKF
! V1_20        2012-06-18 Andreas Rhodin
!  change comment lines
! V1_22        2013-02-13 Harald Anlauf
!  changes for GRIB API / CDI
! V1_23        2013-03-26 Andreas Rhodin
!  LETKF: don't write local ensemble extension (253) for deterministic run
! V1_26        2013/06/27 Harald Anlauf
!  Changes for GRIB2/GRIB_API and ICON.
!  FTRACE instrumentation. Optimize write performance when using GRIB_API.
!  disable uuid for __ibm__
! V1_27        2013-11-08 Harald Anlauf
!  Fixes for generalized vertical coordinate, GRIB_API, centre=Meteoswiss
! V1_28        2014/02/26 Harald Anlauf
!  Fixes for GRIB_API, GRIB1 and GRIB2 handling; clean up for ICON, IAU scheme
! V1_37        2014-12-23 Harald Anlauf
!  bugfix/workaround for GRIB_API 1.12.3 tables and for ECMWF
! V1_42        2015-06-08 Harald Anlauf
!  ICON local patch
! V1_43        2015-08-19 Harald Anlauf
!  changes for GRIB2 and FLAKE
! V1_44        2015-09-30 Andreas Rhodin
!  LETKF determ.first guess: re-label fields from SST-'analysis' to 'forecast'
! V1_45        2015-12-15 Harald Anlauf
!  Cosmetic changes
! V1_46        2016-02-05 Harald Anlauf
!  GRIB2: use tables 4.0, 4.7; DWD local definiton 252, encoding of mean, spread
!  Fix timeRangeIndicator for KENDA, GRIB_API/GRIB1;
!  properly set number of bits for encoded values, default ensemble id
!  relabel_sstana: include w_sno,t_snow,w_i,freshsnw,rho_snow,h_snow
! V1_47        2016-06-06 Harald Anlauf
!  Enable GRIB encoding of selected analysis fields with 24 bits
! V1_48        2016-10-06 Andreas Rhodin
!  set_dwd_defaults: make all arguments optional
! V1_50        2017-01-09 Andreas Rhodin
!  change GRIB inventory printout (5 columns for time range)
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Luis Kornblueh  MPIfM  2001  original code
! Andreas Rhodin  DWD    2002  extended, debugged
!                        2003  now interfaces EMOS instead of wgrib
!                              writing of GRIB-files implemented
!                        2003  use DWD GRIB tables
! Harald Anlauf   DWD    2007  bugfixes, ECMWF ensemble, large file support
!-----------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------

  !=============
  ! Modules used
  !=============
  USE mo_kind,       ONLY: wp               ! working precision kind parameter
  USE mo_exception,  ONLY: finish           ! abort on error condition
  USE mo_mpi_dace,   ONLY: dace,           &! DACE default communicator info
                           self,           &! SELF         communicator info
                           p_bcast          ! overloaded broadcast routine
  USE mo_emos_grib1, ONLY: t_grib1,        &! GRIB record data type
!                          local_dwd,      &! DWD   local extension present ?
!                          local_ecmwf,    &! ECMWF local extension present ?
                           pbopen,         &! Open a GRIB file
                           pbclose,        &! Close a GRIB file
                           pbsetraw,       &! Set record header (dwd or raw)
!                          pbgrib,         &! Read next GRIB product
                           pbwrite,        &! Stream output
                           gribex,         &! Encode/Decode GRIB record
                           reallocate_data,&! reallocate internal data buffer
                           reallocate_buffer,&! allocate internal GRIB buffer
                           destruct,       &! deallocate buffers
                           INVALID_HANDLE   ! invalid GRIB2-API handle
  USE mo_wmo_tables
  USE mo_time,       ONLY: t_time,         &! date+time data type
!                          init_time,      &! initialisation routine
                           cdate,          &! derive date character string
                           ctime,          &! derive time character string
                           ihhmm,          &! derive hours/minutes (integer)
                           ihhhmm,         &! derive hours/minutes (integer)
                           ihhmmss,        &! derive hours/min/sec (integer)
                           iyyyy,          &! derive years         (integer)
                           imm,            &! derive months        (integer)
                           idd,            &! derive days          (integer)
                           ihh,            &! derive hours         (integer)
                           imi,            &! derive minutes       (integer)
                           days,           &! derive days          (real)
                           hours,          &! derive hours         (real)
                           minutes,        &! derive minutes       (real)
                           seconds,        &! derive seconds       (real)
                           OPERATOR(-),    &! subtract times
                           OPERATOR(+),    &! add      times
                           OPERATOR(==),   &! compare  times
                           zero_time        ! time variable with value zero
  USE mo_version,    ONLY: major, minor     ! 3dvar version number
  USE mo_dace_string,ONLY: split            ! split words into character array
  USE mo_grib_invt,  ONLY: t_inventory,    &! Inventory Data Type
                           t_ds,           &! characteristics of the data set
                           t_ct,           &! identification of the center
                           t_pa,           &! characteristics of the parameter
                           t_ti,           &! time information
                           t_lv,           &! level information
                           t_gr,           &! horizontal grid
                           t_en,           &! ensemble data
                           OPERATOR(==),   &! compare  inventory entries
!                          ver_time,       &! get verification time from grib1
                           ref_time,       &! get reference time from grib1
!                          db_time,        &! get data-base time from grib1
                           inventory_entry,&! derive table entry from GRIB
                           grib_library,   &! GRIB API to use
                           ENS_MEAN,       &! ensemble mean
                           ENS_SPREAD       ! ensemble spread
  USE mo_grib12,     ONLY: read_gribex,    &! read + decode GRIB 1 or 2
                           grib1togrib2     ! convert GRIB 1 -> GRIB 2
  USE mo_grib12_dwd, ONLY: search_grib1,   &! Search parameter in GRIB1 table
                           search_grib2,   &! Search parameter in GRIB2 table
                           t_par_grib1,    &! Derived type for GRIB1 parameters
                           t_par_grib2      ! Derived type for GRIB2 parameters
!                          setup_tables_grib12 ! Set up GRIB1/2 tables
#ifdef GRIB_API
  USE grib_api,      ONLY: grib_get_message_size,&
                           grib_copy_message,    &
                           grib_new_from_samples,&
                           grib_clone,           &
                           grib_release,         &
                           grib_get,             &
                           grib_set,             &
                           grib_set_missing,     &
                           grib_gribex_mode_on,  &
                           grib_open_file,       &
                           grib_close_file,      &
                           grib_write,           &
                           kindOfSize
  USE mo_t_grib_api, ONLY: GRIB_4_0_ANA_FC      ,&! Analysis or forecast
                           GRIB_4_0_FC_ENS      ,&! Individual ensemble forecast
                           GRIB_4_0_FC_ENS_ALL  ,&! Derived forecast
                           GRIB_4_0_ANA_FC_CHEM ,&! Analysis or forecast for chemical constituents
                           GRIB_4_0_CHEM_ENS    ,&! Individual ensemble forecast for chemical constituents
!                          GRIB_4_0_AERO        ,&! Aerosol analysis or forecast
!                          GRIB_4_0_AERO_ENS    ,&! Aerosol individual ensemble forecast
                           GRIB_4_3_ANA,         &! Tab.4.3: analysis
                           GRIB_4_3_FC,          &! Tab.4.3: forecast
                           GRIB_4_3_ANA_ERR,     &! Tab.4.3: analysis error
                           GRIB_4_3_ANA_INC     ,&! Tab.4.3: analysis increment
                           GRIB_4_7_MEAN        ,&! Unweighted mean of all members
                           GRIB_4_7_SPREAD        ! Spread of all members
#endif
  IMPLICIT NONE

  !================
  ! Public entities
  !================
  PRIVATE
  !------------------------
  ! predefined ensemble ids
  !------------------------
  PUBLIC :: EID_DETERM       ! deterministic run, no ensemble
  PUBLIC :: EID_ROUTINE      ! COSMO-DE routine ensemble
  PUBLIC :: EID_SREPS        ! SREPS boundary conditions
  PUBLIC :: EID_BCEPS        ! BCEPS boundary conditions
  PUBLIC :: EID_GME          ! GME   boundary conditions
  !----------------------------
  ! ensemble derived quantities
  !----------------------------
  PUBLIC :: ENS_MEAN         ! ensemble mean
  PUBLIC :: ENS_SPREAD       ! ensemble spread
  !-------------------------------------------------
  ! data type definition: t_inventory and components
  !-------------------------------------------------
  PUBLIC :: t_inventory      ! Inventory Data Type
  PUBLIC ::   t_ds           !   characteristics of the data set
  PUBLIC ::   t_ct           !   identification of the center
  PUBLIC ::   t_pa           !   characteristics of the parameter
  PUBLIC ::   t_ti           !   time information
  PUBLIC ::   t_lv           !   level information
  PUBLIC ::   t_gr           !   horizontal grid
  PUBLIC ::   t_en           !   ensemble data
  !--------------------------------------------------
  ! routines to derive the inventory from a GRIB file
  !--------------------------------------------------
  PUBLIC :: get_inventory    ! read inventory table
  PUBLIC :: print_inventory  ! print inventory
  PUBLIC :: print_legend     ! print legend
  PUBLIC :: relabel_fields   ! relabel 'analysis' to 'forecast'
  PUBLIC :: p_bcast          ! broadcast inventory via MPI
  !------------------
  ! write a GRIB file
  !------------------
  PUBLIC :: open_gribfile    ! open a file, set default grib sections
  PUBLIC :: write_gribrecord ! encode a record and write
  PUBLIC :: close_gribfile   ! close the file
  !---------------------------------
  ! prepare a GRIB record for output
  !---------------------------------
  PUBLIC :: set_defaults     ! specify default values
  PUBLIC :: set_dwd_defaults ! specify dwd local part of PDB
  PUBLIC :: set_dwd_def_ens  ! specify dwd local PDB for ensemble prediction
  PUBLIC :: set_center       ! specify center, subcenter, process
  PUBLIC :: set_dwd_local    ! specify dwd local part of PDB
  PUBLIC :: set_dwd_ens      ! specify dwd local PDB for ensemble prediction
  PUBLIC :: set_ecmwf_def_ens! specify defaults for ECMWF local extension
  PUBLIC :: set_ecmwf_ens    ! specify ECMWF local extension for ensemble fc.
  PUBLIC :: set_ref_time     ! specify reference time
  PUBLIC :: set_ver_time     ! specify verification time
  PUBLIC :: set_grid_vert    ! specify leveltype, vertical coordinate params.
  PUBLIC :: set_grid_latlon  ! specify lat/lon grid
  PUBLIC :: set_grid_gauss   ! specify Gaussian grid
  PUBLIC :: set_grid_tri     ! specify triangular grid
  PUBLIC :: set_grid_icon    ! specify ICON grid
  PUBLIC :: set_code         ! specify code, table, bits
  PUBLIC :: set_level        ! specify level
  PUBLIC :: set_data         ! transfer data
  PUBLIC :: set_grib         ! select GRIB edition and library
  PUBLIC :: grib_edition     ! GRIB edition to write
! PUBLIC :: grib_library     ! GRIB library version (GRIBEX or GRIB API)
  PUBLIC :: grib_api_version ! GRIB API library version numbers
  PUBLIC :: gribtabversion   ! GRIB2 minimum tablesVersion
  PUBLIC :: grib_packing     ! GRIB packing or compression
  !===========
  ! Interfaces
  !===========
  !---------------------------------
  ! subroutine print_inventory (inv)
  !---------------------------------
  INTERFACE print_inventory
    MODULE PROCEDURE print_inv_0 ! type (t_inventory) :: inv
    MODULE PROCEDURE print_inv_1 ! type (t_inventory) :: inv (:)
  END INTERFACE print_inventory

  !-------------------------
  ! subroutine p_bcast (inv)
  !-------------------------
  INTERFACE p_bcast
    MODULE PROCEDURE p_bcast_invt
  END INTERFACE p_bcast

  !---------------------
  ! subroutine set_level
  !---------------------
  INTERFACE set_level
!   MODULE PROCEDURE set_level1 ! set levels using GRIB1 conventions
    MODULE PROCEDURE set_level2 ! set levels using GRIB2 conventions
  END INTERFACE set_level

  !--------------------
  ! subroutine set_data
  !--------------------
  INTERFACE set_data
    MODULE PROCEDURE set_data_1 ! pass rank 1 array
    MODULE PROCEDURE set_data_2 ! pass rank 2 array
    MODULE PROCEDURE set_data_3 ! pass rank 3 array
  END INTERFACE set_data

#ifdef GRIB_API
  !----------------------------------
  ! retrieve GRIB_API version triplet
  !----------------------------------
  interface
     function grib_get_api_version () bind(C,NAME="grib_get_api_version")
       use, intrinsic  :: iso_c_binding, only: C_LONG
       integer(C_LONG) :: grib_get_api_version
     end function grib_get_api_version
  end interface
#endif

  !=================
  ! module variables
  !=================
  !----------------------------------
  ! Some default values are kept here
  !----------------------------------
  integer ,protected :: grib_edition   =  1   ! GRIB edition to write (-1,1,2)
  integer            :: gribtabversion = -1   ! GRIB2 minimum tablesVersion
  character(len=64)  :: grib_packing   = ""   ! GRIB packing (default:none)
  !------------------------
  ! predefined ensemble ids
  !------------------------
  INTEGER, PARAMETER :: EID_DETERM  =   0     ! deterministic run, no ensemble
  INTEGER, PARAMETER :: EID_ROUTINE =   1     !       routine ensemble
  INTEGER, PARAMETER :: EID_SREPS   = 201     ! SREPS boundary conditions
  INTEGER, PARAMETER :: EID_BCEPS   = 202     ! BCEPS boundary conditions
  INTEGER, PARAMETER :: EID_GME     = 203     ! GME   boundary conditions

  INTEGER            :: d_c     = 255         ! Identification of centre.
  INTEGER            :: d_s     = 255         ! Sub-centre identifier.
  INTEGER            :: d_p     = 255         ! Generating process id. number.
  INTEGER            :: d_ptyp  = 255         ! Generating process type.
  INTEGER            :: d_b     =  16         ! Number of bits used for encoding
  INTEGER            :: d_r     =   1         ! raw flag: 0=dwd format, 1=raw
  INTEGER            :: d_l     =   0         ! local flag
  INTEGER            :: d_run   =   3         ! run type
  INTEGER            :: d_ran   =   0         ! time range indicator for analysis
  INTEGER            :: d_e     =   0         ! experiment number
  TYPE(t_time) ,SAVE :: d_t                   ! run time
  TYPE(t_time) ,SAVE :: d_rt                  ! reference time
  LOGICAL            :: d_wlev  = .false.     ! write levels in GRIB section 2
  INTEGER            :: d_n_en  =   0         ! ensemble size
  INTEGER            :: d_i_en  =   0         ! ensemble member number
  INTEGER            :: d_id_en = EID_ROUTINE ! ensemble id
  character(len=4)   :: d_expid = "e???"      ! (ECMWF) Experiment id.

  INTEGER ,PARAMETER :: d_edition = 1         ! fallback GRIB edition to write
  logical            :: d_pbio    = .false.   ! Prefer PBIO over grib_api/write

CONTAINS
!==============================================================================
  SUBROUTINE set_grib (edition, library)
    INTEGER ,OPTIONAL ,INTENT(in) :: edition    ! GRIB edition to write
    INTEGER ,OPTIONAL ,INTENT(in) :: library    ! GRIB API to use

    if (present (edition)) grib_edition = edition
    if (present (library)) grib_library = library
  END SUBROUTINE set_grib
!------------------------------------------------------------------------------
  SUBROUTINE get_inventory (inventory, file, pio, comm, ifile, append)
  !-----------------------------------------------------------------------
  ! get the inventory of a file
  !
  ! The inventory pointer is allocated with the appropriate size.
  ! The pointer status must be defined (nullified or allocated) on entry.
  ! If the file is not present, a zero size inventory is returned.
  !
  ! By default (pio not present) the file is scanned on PE dace% pio.
  ! By default (comm not present) the inventory is broadcasted to all PEs.
  ! (pass comm=self% com to keep inventory on PE pio only).
  !-----------------------------------------------------------------------
  TYPE (t_inventory) ,POINTER    :: inventory(:) ! GRIB file inventory
  CHARACTER (len=*)  ,INTENT(in) :: file         ! GRIB file name
  INTEGER ,OPTIONAL  ,INTENT(in) :: pio          ! I/O processor index
  INTEGER ,OPTIONAL  ,INTENT(in) :: comm         ! communicator
  INTEGER ,OPTIONAL  ,INTENT(in) :: ifile        ! file index
  LOGICAL ,OPTIONAL  ,INTENT(in) :: append       ! append to inventory
    !----------------
    ! local variables
    !----------------
    INTEGER                     :: ierror           ! error return code
    TYPE (t_inventory) ,POINTER :: tmp1(:), tmp2(:) ! temporary storage
    TYPE (t_inventory)          :: gi
    INTEGER                     :: lpio             ! local copy of pio
    INTEGER                     :: n                ! record count
    INTEGER                     :: lfile            ! file index
    LOGICAL                     :: lappend          ! append to inventory
    INTEGER          ,PARAMETER :: defaultsize = 16384

    TYPE (t_grib1) :: grib
    !--------------------
    ! read on one PE only
    !--------------------
    lpio     = dace% pio ;if (present(pio   )) lpio    = pio
    lfile    = 0         ;if (present(ifile )) lfile   = ifile
    lappend  = .false.   ;if (present(append)) lappend = append
    IF (dace% pe == lpio) THEN
      !---------------------------------------------------
      ! reuse memory if 'inventory' pointer is associated,
      ! else allocate temporary
      !---------------------------------------------------
      n = 0
      IF (ASSOCIATED (inventory)) THEN
        tmp1           =>     inventory
        if (lappend) n = size(inventory)
      ELSE
        ALLOCATE (tmp1(defaultsize))
      ENDIF
      !----------
      ! open file
      !----------
      CALL open_gribfile (grib, TRIM(file), 'r', kret=ierror)
      IF (ierror /= 0) goto 99
      !-----------------------------------------
      ! read inventory, one grib block at a time
      !-----------------------------------------
      DO
        CALL read_gribex (grib, 'J', kret=ierror)
!       CALL pbgrib (grib, ierror)
!       IF (ierror /= 0) EXIT
!       CALL gribex (grib, 'J', ierror)
        IF (ierror /= 0) EXIT
        n = n + 1
        !-----------------------------------------
        ! derive inventory entry from GRIB1 record
        !-----------------------------------------
        gi = inventory_entry (grib)

        !-----------------------------------
        ! enlarge pointer array if too small
        !-----------------------------------
        IF (n > SIZE(tmp1)) THEN
          tmp2 => tmp1
          ALLOCATE (tmp1 (max (defaultsize, 2*SIZE (tmp2))))
          tmp1 (1:SIZE(tmp2)) = tmp2
          DEALLOCATE (tmp2)
        ENDIF
        tmp1 (n)        = gi
        tmp1 (n)% ifile = lfile
      END DO
      !-----------
      ! close file
      !-----------
      CALL close_gribfile (grib)
      !-------------------------------------------
      ! allocate return variable with correct size
      !-------------------------------------------
99    CONTINUE
      IF (SIZE(tmp1) /= n) THEN
        ALLOCATE (inventory (n))
        inventory = tmp1 (1:n)
        DEALLOCATE (tmp1)
      ELSE
        inventory => tmp1
      ENDIF

      call match_entries (inventory)

    ENDIF

    !----------------------------
    ! broadcast result to all PEs
    !----------------------------
    CALL p_bcast (inventory, lpio, comm)

  END SUBROUTINE get_inventory
!------------------------------------------------------------------------------
  subroutine match_entries (inventory)
  TYPE (t_inventory) ,intent(inout):: inventory(:) ! GRIB file inventory

    LOGICAL        ,ALLOCATABLE :: mask (:)         ! mask for entries
    INTEGER                     :: n                ! record count
    INTEGER                     :: i,j,l
    INTEGER                     :: nvar, ntime, nlev, nmem
    INTEGER                     :: ivar

    n = size (inventory)
    ALLOCATE (mask (n))

    !
    ! same variable
    !
    inventory% nvar = 0
    DO i=1,n
      IF (inventory(i)% nvar /= 0) CYCLE
      mask = (inventory(i)% pa == inventory% pa)
      nvar = COUNT (mask)
      WHERE (mask) inventory% nvar = nvar
    END DO

    !
    ! same time
    !
    inventory% ntime = 0
    DO i=1,n
      IF (inventory(i)% ntime /= 0) CYCLE
      mask = (inventory(i)% pa == inventory% pa &
        .AND. inventory(i)% lv == inventory% lv &
        .AND. inventory(i)% en == inventory% en )
      ntime = COUNT (mask)
      l = 0
!NEC$ ivdep
      DO j = 1,n
        IF (mask(j)) THEN
          l = l+1
          inventory(j)% ntime = ntime
          inventory(j)% jtime = l
        END IF
      END DO
    END DO

    !
    ! Same level.  Discriminate between levels and layers.
    !
    inventory% nlev = 0
    !-------
    ! Levels
    !-------
    ivar = 0
    DO i=1,n
      IF (inventory(i)% nlev /= 0) CYCLE
      if (inventory(i)% lv% levels(1) /= inventory(i)% lv% levels(2)) cycle ! Layer?
      ivar = ivar + 1
      mask = (inventory(i)% pa            == inventory% pa &
        .AND. inventory(i)% ti            == inventory% ti &
        .AND. inventory(i)% en            == inventory% en &
        .AND. inventory(i)% lv% leveltype == inventory% lv% leveltype)
      mask = mask .and. (inventory% lv% levels(1) == inventory% lv% levels(2))
      nlev = COUNT (mask)
      l = 0
!NEC$ ivdep
      DO j = 1,n
        IF (mask(j)) THEN
          l = l+1
          inventory(j)% nlev = nlev
          inventory(j)% jlev = l
          inventory(j)% ivar = ivar
        END IF
      END DO
    END DO
    !-------
    ! Layers
    !-------
    ivar = 0
    DO i=1,n
      IF (inventory(i)% nlev /= 0) CYCLE
      if (inventory(i)% lv% levels(1) == inventory(i)% lv% levels(2)) cycle ! Level?
      ivar = ivar + 1
      mask = (inventory(i)% pa            == inventory% pa &
        .AND. inventory(i)% ti            == inventory% ti &
        .AND. inventory(i)% en            == inventory% en &
        .AND. inventory(i)% lv% leveltype == inventory% lv% leveltype)
      mask = mask .and. (inventory% lv% levels(1) /= inventory% lv% levels(2))
      nlev = COUNT (mask)
      l = 0
!NEC$ ivdep
      DO j = 1,n
        IF (mask(j)) THEN
          l = l+1
          inventory(j)% nlev = nlev
          inventory(j)% jlev = l
          inventory(j)% ivar = ivar
        END IF
      END DO
    END DO

    !
    ! same ensemble member
    !
    inventory% nmem = 0
    inventory% jmem = 0
    DO i=1,n
      IF (inventory(i)% nmem /= 0) CYCLE
      IF (inventory(i)% en% size == 0 .AND. &
          inventory(i)% en% no   == 0 .AND. &
          inventory(i)% en% id   == 0 ) CYCLE
      mask = (inventory(i)% pa == inventory% pa &
        .AND. inventory(i)% ti == inventory% ti &
        .AND. inventory(i)% lv == inventory% lv )
      nmem = COUNT (mask)
      l = 0
!NEC$ ivdep
      DO j = 1,n
        IF (mask(j)) THEN
          l = l+1
          inventory(j)% nmem = nmem
          inventory(j)% jmem = l
        END IF
      END DO
    END DO

    DEALLOCATE (mask)

  end subroutine match_entries
!------------------------------------------------------------------------------
  subroutine relabel_fields (invt, par)
  type(t_inventory) ,intent(inout) :: invt (:) ! inventory
  character(len=*)  ,intent(in)    :: par      ! parameters to relabel
  !----------------------------------------------------------
  ! relable SST analysis fields (t_so) and first guess fields
  ! so that the analysis fields are preferred for reading
  !----------------------------------------------------------
    integer          :: j      ! variable index
    integer          :: i1, i2 ! inventory entry index
    logical          :: change (size(invt))
    character(len=8) :: var(20)
    integer          :: n
    !-----------------------------
    ! mark entries to be relabeled
    !-----------------------------
    call split (var, par, n)
    change = .false.
    do j = 1, n
      do i1 = 1, size (invt)
        if (invt(i1)% pa% iname   == var(j) .and. &
            invt(i1)% pa% runtype == 'analysis'   ) then
          change(i1) = .true.
          do i2 = 1, size (invt)
            if (invt(i2)% lv% levels(1) == invt(i1)% lv% levels(1) .and. &
                invt(i2)% ti% ver_time  == invt(i1)% ti% ver_time  .and. &
                invt(i2)% pa% iname     == var(j)                  .and. &
                invt(i2)% pa% runtype   == 'forecast'                    ) then
              change(i2) = .true.
              exit
            endif
          end do
          exit
        endif
      end do
    end do
    !-------------------
    ! print out original
    !-------------------
    if (dace% lpio) then
      write(6,*)
      write(6,'(a)') repeat('-',79)
      write(6,*)
      write(6,*)' relabel_fields: original inventory:'
      write(6,*)
    endif
    call print_inventory (invt, first=.true., mask=change)
    !---------
    ! re-label
    !---------
    where (change .and. invt(:)% pa% runtype == 'forecast') &
                        invt(:)% pa% runtype =  'unused'
    where (change .and. invt(:)% pa% runtype == 'analysis') &
                        invt(:)% pa% runtype =  'forecast'
    call match_entries   (invt)
    !------------------
    ! print out changes
    !------------------
    if (dace% lpio) then
      write(6,*)
      write(6,*)' relabel_fields: modified inventory:'
      write(6,*)
    endif
    call print_inventory (invt, first=.true., mask=change)
    if (dace% lpio) write(6,*)

  end  subroutine relabel_fields
!------------------------------------------------------------------------------
  SUBROUTINE p_bcast_invt (buffer, source, comm)
  !----------------------------------------------------
  ! Broadcast GRIB file inventory.
  ! by default ('comm' missing) dace% com is used.
  !----------------------------------------------------
  TYPE(t_inventory) ,POINTER    :: buffer (:) ! inventory to broadcast
  INTEGER           ,INTENT(in) :: source     ! source processor index
  INTEGER ,OPTIONAL ,INTENT(in) :: comm       ! communicator

    TYPE (t_inventory) :: dummy
    INTEGER            :: size_invt
    INTEGER            :: count
    INTEGER            :: lcomm

    IF (dace% npe == 1) RETURN

    size_invt = SIZE (TRANSFER ( dummy,(/' '/)))
    lcomm = dace% comm; IF(PRESENT(comm)) lcomm = comm
    if (lcomm == self% comm) return

    IF (dace% pe == source) count = SIZE (buffer)
    CALL p_bcast (count, source, lcomm)
    IF (ASSOCIATED (buffer)) THEN
      IF (SIZE (buffer) /= count) DEALLOCATE (buffer)
    ENDIF
    IF (.NOT.ASSOCIATED (buffer)) ALLOCATE (buffer(count))
    CALL p_bcast_derivedtype (buffer, count*size_invt, source, lcomm)

  END SUBROUTINE p_bcast_invt
!==============================================================================
  SUBROUTINE print_inv_0 (i, liname, ldbtime, llname)
  TYPE(t_inventory) ,INTENT(in) :: i
  LOGICAL ,OPTIONAL ,INTENT(in) :: liname  ! print internal name
  LOGICAL ,OPTIONAL ,INTENT(in) :: ldbtime ! print data base time,not obs.time
  LOGICAL ,OPTIONAL ,INTENT(in) :: llname  ! print longname
    character(len=20)  :: name
    character(len=128) :: longname
    logical            :: li, ldb, ll
    integer            :: vv, rg
    character(len=4)   :: runclass
    type(t_time)       :: time     ! verification or data-base time to print
    character(len=8)   :: triple   ! discipline.category.parameter# (len=8)
    character(len=11)  :: buffer   ! auxiliary buffer
    character(len=2)   :: ed       ! edition (annotated for packing/compression)
    type(t_par_grib1)  :: t1       ! GRIB 1 table
    type(t_par_grib2)  :: t2       ! GRIB 2 table

    li   = .false.;         if (present(liname )) li   = liname
    ldb  = .false.;         if (present(ldbtime)) ldb  = ldbtime
    ll   = .false.;         if (present(llname )) ll   = llname
    name = i%pa% shortname; if (li)               name = i%pa% iname
    time = i%ti% ver_time;  if (ldb)              time = i%ds% db_time
    vv   = ihhhmm (i%ti% ver_time - i%ti% ref_time)
    rg   = ihhhmm (i%ti% rng_time)

    select case (i%pa% runclass)
    case (0)
      runclass = 'main'
    case (1)
      runclass = 'pre'
    case (2)
      runclass = 'ass'
    case default
      write (runclass,'(i4)') i%pa% runclass
    end select

    write(ed,'(i2)') i%ds% edition
    if (i%ds% packing /= 0) then
      ed(1:1) = "*"                     ! mark compressed records in inventory
    end if

    if (i%ds% edition == 1) then
      write(triple,'(i3,"|",i3)')        i%pa% table, i%pa%code
      triple(8:8) = " "
    else
      write(buffer,'(i0,".",i0,".",i0)') i%pa% discipline, i%pa% category, &
                                         i%pa% number
      triple = buffer(1:8)
      if (buffer(9:9) /= " ") triple = "********"
    end if

    if (ll) then
      longname = ''
      select case (i%ds% edition)
      case (1)
        t1 = search_grib1 (name=name)
        longname = t1% longname
      case (2)
        t2 = search_grib2 (name=name)
        longname = t2% longname
      end select
      write (6,'(i5,"|",i12,"|",a,"|",a,"| ",a)')    &
        i%ds% rec, i%ds% ptr, ed, name, trim(longname)
    else
      select case (i%lv% leveltype)
      case (100,103,105,111)
        write (6,'(i5,"|",i12,"|",a,"|",a,"|",a10,"|",a8,"|",a8,5("|",i3),"|",i7,&
             &"|",a,"|",i5,"|",i5,"|",i4,5("|",i3),"|",a4,"|",i5)')              &
          i%ds% rec, i%ds% ptr, ed, name, cdate(time), ctime(time),              &
          triple, i%ct% center, i%ct% subcenter, i%ct% process,                  &
          i%pa% zen, i%lv% leveltype, i%lv% levelvalue,                          &
          i%pa% runtype, vv, rg, i%gr% ni, i%gr% gridtype,                       &
          i%nlev,i%en%no, i%nmem, i%en%id, runclass, i%pa% expid
      case default
        write (6,'(i5,"|",i12,"|",a,"|",a,"|",a10,"|",a8,"|",a8,7("|",i3),&
             &"|",a,"|",i5,"|",i5,"|",i4,5("|",i3),"|",a4,"|",i5)')       &
          i%ds% rec, i%ds% ptr, ed, name, cdate(time), ctime(time),       &
          triple, i%ct% center, i%ct% subcenter, i%ct% process,           &
          i%pa% zen, i%lv% leveltype, max (i%lv% levels(:), 0),           &
          i%pa% runtype, vv, rg, i%gr% ni, i%gr% gridtype,                &
          i%nlev,i%en%no, i%nmem, i%en%id, runclass, i%pa% expid
      end select
    endif

  END SUBROUTINE print_inv_0
!------------------------------------------------------------------------------
  SUBROUTINE print_inv_1 (inv, first, mask, liname, ldbtime, llname)
  !=============================================================
  ! Write the inventory table on the I/O PE on stdout (unit 6).
  !
  ! If 'first' is present and true only one field of a sequence
  ! (e.g. one level, one ensemble member, etc.) is printed.
  !=============================================================
  TYPE (t_inventory) ,INTENT(in) :: inv (:) ! inventory table to print
  LOGICAL  ,OPTIONAL ,INTENT(in) :: first   ! print first/last code entry
  LOGICAL  ,OPTIONAL ,INTENT(in) :: mask(:) ! mask to select table entries
  LOGICAL  ,OPTIONAL ,INTENT(in) :: liname  ! print internal name
  LOGICAL  ,OPTIONAL ,INTENT(in) :: ldbtime ! print data base time
  LOGICAL  ,OPTIONAL ,INTENT(in) :: llname  ! print longname
    !----------------
    ! local variables
    !----------------
    INTEGER :: i
    LOGICAL :: fl
    LOGICAL :: ll
    !----------------
    ! write on I/O PE
    !----------------
    IF (dace% lpio) THEN
      !-------------
      ! write header
      !-------------
      ll = .false.; if (present(llname)) ll = llname
      if (ll) then
        write (6,'(a)') &
          "  rec      address ed name                 description"
      else if (all (inv(:)%ds % edition == 2)) then
        write (6,'(a)') &
          "  rec      address ed name                 valid date     time triple   cnt sub prc &
          &zen lvt levelvl  runtype  vvmm range   ni grd nlv mem nme eid clas expid"
      else
        write (6,'(a)') &
          "  rec      address ed name                 valid date     time tab cde  cnt sub prc &
          &zen lvt levelvl  runtype  vvmm range   ni grd nlv mem nme eid clas expid"
      end if
      !------------
      ! write lines
      !------------
      fl = .FALSE.; IF (PRESENT(first)) fl = first
      DO i=1,SIZE(inv)
!       IF(fl .AND. ANY(inv(:i-1)%pa% code == inv(i)%pa% code) &
!             .AND. ANY(inv(i+1:)%pa% code == inv(i)%pa% code) ) CYCLE
        IF (fl .AND. inv(i)% jlev /= 1) CYCLE
        IF (fl .AND. inv(i)% jmem >  1) CYCLE
        IF (PRESENT (mask)) THEN
          IF (.NOT.mask(i)) CYCLE
        ENDIF
        CALL print_inventory (inv (i), liname, ldbtime, llname)
      END DO

    ENDIF

  END SUBROUTINE print_inv_1
!------------------------------------------------------------------------------
  subroutine print_legend
    if (dace% lpio) then
      write (6,*)
      write (6,*) '  Legend:'
      write (6,*)
      write (6,*) '    rec     : record number in file'
      write (6,*) '    address : address in file'
      write (6,*) '    ed      : GRIB edition (1 or 2)'
      write (6,*) '    date    : valid date'
      write (6,*) '    time    : valid time'
      write (6,*) '    tab     : GRIB table number (GRIB 1)'
      write (6,*) '    cde     : GRIB code number  (GRIB 1)'
      write (6,*) '    triple  : GRIB triple       (GRIB 2)'
      write (6,*) '    cnt     : center'
      write (6,*) '    sub     : subcenter'
      write (6,*) '    prc     : generating process'
      write (6,*) '    zen     : DWD Zusatzelementnummer'
      write (6,*) '    lvt     : level type (vertical grid)'
      write (6,*) '    levelvl : level value(s)'
      write (6,*) '    runtype : run type (forecast, ass, ..)'
      write (6,*) '    vvmm    : forecast lead time (hhhmm)'
      write (6,*) '    range   : time range (hhmm) for accumulated/averaged fields'
      write (6,*) '    ni      : grid resolution parameter'
      write (6,*) '    grd     : grid type (horizontal grid)'
      write (6,*) '    nlv     : number of levels'
      write (6,*) '    mem     : ensemble member id'
      write (6,*) '    nme     : number of members in file'
      write (6,*) '    eid     : ensemble id'
      write (6,*) '    clas    : run class (main, ass, pre)'
      write (6,*) '    expid   : experiment id (1=routine, 50=para, 51=para1)'
      write (6,*)
    endif
  end subroutine print_legend
!==============================================================================
!    SUBROUTINE decode_bds_float_2 (x)
!    REAL(sp), INTENT(out) :: x(:,:) ! data array to read
!      INTEGER :: nxny
!      nxny = size(x)
!      call decode_bds_float (nxny, x)
!    END SUBROUTINE decode_bds_float_2
!------------------------------------------------------------------------------
!    SUBROUTINE decode_bds_double_2 (x)
!    REAL(dp), INTENT(out) :: x(:,:) ! data array to read
!      INTEGER :: nxny
!      nxny = size(x)
!      call decode_bds_double (nxny, x)
!    END SUBROUTINE decode_bds_double_2
!------------------------------------------------------------------------------
!    SUBROUTINE decode_bds_float_3 (x)
!    REAL(sp), INTENT(out) :: x(:,:,:) ! data array to read
!      INTEGER :: nxny
!      nxny = size(x)
!      call decode_bds_float (nxny, x)
!    END SUBROUTINE decode_bds_float_3
!------------------------------------------------------------------------------
!    SUBROUTINE decode_bds_double_3 (x)
!    REAL(dp), INTENT(out) :: x(:,:,:) ! data array to read
!      INTEGER :: nxny
!      nxny = size(x)
!      call decode_bds_double (nxny, x)
!    END SUBROUTINE decode_bds_double_3
!==============================================================================
  !=========================================
  ! routines to set a GRIB record for output
  !=========================================
!------------------------------------------------------------------------------
  !--------------------------------------------------------
  ! Set some default values valid for all output GRIB files
  !--------------------------------------------------------
  SUBROUTINE set_defaults (center, subcenter, process, process_type, &
                           bits, raw, ref_time, wlev)
  INTEGER ,INTENT(in) ,OPTIONAL :: center        ! Identification of centre
  INTEGER ,INTENT(in) ,OPTIONAL :: subcenter     ! Sub-centre identifier
  INTEGER ,INTENT(in) ,OPTIONAL :: process       ! Generating process id
  INTEGER ,INTENT(in) ,OPTIONAL :: process_type  ! Generating process type
  INTEGER ,INTENT(in) ,OPTIONAL :: bits          ! Bits used for encoding
  INTEGER ,INTENT(in) ,OPTIONAL :: raw           ! 0=dwd format, 1=raw
  LOGICAL ,INTENT(in) ,OPTIONAL :: wlev          ! write levels to sec.2
  TYPE(t_time) ,INTENT(in) ,OPTIONAL :: ref_time ! Reference time

    IF (PRESENT (center)      ) d_c    = center
    IF (PRESENT (subcenter)   ) d_s    = subcenter
    IF (PRESENT (process)     ) d_p    = process
    IF (PRESENT (process_type)) then
                                d_ptyp = process_type
    else
                                d_ptyp = d_p     ! Same as process id.
    end if
    IF (PRESENT (bits)        ) d_b    = bits
    IF (PRESENT (raw)         ) d_r    = raw
    IF (PRESENT (ref_time)    ) d_rt   = ref_time
    IF (PRESENT (wlev)        ) d_wlev = wlev

  END SUBROUTINE set_defaults
!------------------------------------------------------------------------------
  !-------------------------------------------
  ! specify defaults for DWD local part of PDB
  !-------------------------------------------
  SUBROUTINE set_dwd_defaults (run_time, run_type, nex, range_ana)
  TYPE(t_time),INTENT(in) ,optional :: run_time  ! time of assimilation run
  INTEGER     ,INTENT(in) ,optional :: run_type  ! 0=haupt, 1=pre, 2=ass, 3=test
  INTEGER     ,INTENT(in) ,optional :: nex       ! experiment number
  INTEGER     ,INTENT(in) ,optional :: range_ana ! analysis time range indicator
    if (present (run_time))  d_t   = run_time
    if (present (run_type))  d_run = run_type
    if (present (nex))       d_e   = nex
    if (present (range_ana)) d_ran = range_ana
    select case (d_c)
    case (WMO0_DWD,WMO0_COSMO,WMO0_MSWISS,WMO0_COMET)
       d_l = 254                          ! set local_dwd   if center==DWD
    case (WMO0_ECMWF)
       d_l = 1                            ! set local_ecmwf if center==ECMWF
    end select
  end SUBROUTINE set_dwd_defaults
!------------------------------------------------------------------------------
  SUBROUTINE set_dwd_def_ens (n_ens, i_ens, id_ens)
  INTEGER        ,INTENT(in)    :: n_ens  ! ensemble size
  INTEGER        ,INTENT(in)    :: i_ens  ! ensemble member number
  INTEGER        ,INTENT(in)    :: id_ens ! ensemble id
    d_n_en  = n_ens
    d_i_en  = i_ens
    d_id_en = id_ens
    if (d_c == WMO0_DWD   .or. &
        d_c == WMO0_COSMO .or. &
        d_c == WMO0_MSWISS.or. &
        d_c == WMO0_COMET      ) then
       if (n_ens < 255) then
          d_l = 253                       ! local_dwd, "small" ensemble
       else
          d_l = 153                       ! local_dwd, "big" ensemble
       end if
    end if
  END SUBROUTINE set_dwd_def_ens
!------------------------------------------------------------------------------
  SUBROUTINE set_ecmwf_def_ens (n_ens, i_ens, expid)
  INTEGER         ,INTENT(in)           :: n_ens  ! ensemble size
  INTEGER         ,INTENT(in)           :: i_ens  ! ensemble member number
  CHARACTER(len=*),INTENT(in), optional :: expid
    d_n_en = n_ens
    d_i_en = i_ens
    if (present (expid)) d_expid = expid
    if (d_c == WMO0_ECMWF) d_l   = 1      ! set local_ecmwf if center==ECMWF
  END SUBROUTINE set_ecmwf_def_ens
!------------------------------------------------------------------------------
  !-----------------------------------
  ! specify center, subcenter, process
  !-----------------------------------
  SUBROUTINE set_center (grib, center, subcenter, process)
  TYPE (t_grib1) ,INTENT(inout) :: grib
  INTEGER        ,INTENT(in)    :: center
  INTEGER        ,INTENT(in)    :: subcenter
  INTEGER        ,INTENT(in)    :: process

    grib% isec1% center     = center
    grib% isec1% sub_center = subcenter
    grib% isec1% process    = process

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       call grib_set (grib% handle,'centre'                     ,center)
       call grib_set (grib% handle,'subCentre'                  ,subcenter)
       call grib_set (grib% handle,'generatingProcessIdentifier',process)
    end if
#endif
  END SUBROUTINE set_center
!------------------------------------------------------------------------------
  !------------------------------
  ! specify DWD local part of PDB
  !------------------------------
  SUBROUTINE set_dwd_local (grib, run_time, run_type, nex)
  TYPE (t_grib1) ,INTENT(inout) :: grib     ! grib buffer
  TYPE (t_time)  ,INTENT(in)    :: run_time ! time of assimilation run
  INTEGER        ,INTENT(in)    :: run_type ! 0=haupt, 2=ass, 3=test
  INTEGER        ,INTENT(in)    :: nex      ! experiment number

#ifdef GRIB_API
    integer :: h        ! handle
    integer :: vers     ! local version number
#endif

    grib% isec1 % local_flag    =   1
    grib% isec1 % local_ident   = 254
    grib% s1_dwd% local_ident   = 254
    grib% s1_dwd% day_number    =   0 ! not used
    grib% s1_dwd% record_number =   0 ! not used
    grib% s1_dwd% decoding      =   0 ! not used
    grib% s1_dwd% element_no    =   0 ! not used
    if (run_time == zero_time) then
     grib%s1_dwd% year          = 255
     grib%s1_dwd% month         = 255
     grib%s1_dwd% day           = 255
     grib%s1_dwd% hour          = 255
     grib%s1_dwd% minute        = 255
    else
     grib%s1_dwd% year          =      iyyyy (run_time) - 1900
     grib%s1_dwd% month         =      imm   (run_time)
     grib%s1_dwd% day           =      idd   (run_time)
     grib%s1_dwd% hour          =      ihh   (run_time)
     grib%s1_dwd% minute        =      imi   (run_time)
     if (grib% s1_dwd% year < 0 .or. grib% s1_dwd% year > 255) &
      grib%s1_dwd%year          = 255 ! Mark out-of-bounds values
    end if
    grib% s1_dwd% exp           = nex
    grib% s1_dwd% run_type      = run_type
    grib% s1_dwd% user_id       = 0
    grib% s1_dwd% experiment_id = 0
    grib% s1_dwd% ensemble_id   = d_id_en       ! default: EID_ROUTINE
    grib% s1_dwd% ensemble_size = 0
    grib% s1_dwd% ensemble_no   = 0
    grib% s1_dwd% major_version = major()
    grib% s1_dwd% minor_version = minor()

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       h = grib% handle
       select case (grib% isec0% edition)
       case (1)
          call grib_set (h,'setLocalDefinition'   ,1)
          call grib_set (h,'localDefinitionNumber'  ,grib% isec1% local_ident)
          call grib_set (h,'localElementNumber'     ,grib% s1_dwd% element_no)
          call grib_set (h,'localVersionNumber'     ,nex + 2**14 * run_type)
          call grib_set (h,'localDecodeDateYear'    ,grib% s1_dwd% year)
          call grib_set (h,'localDecodeDateMonth'   ,grib% s1_dwd% month)
          call grib_set (h,'localDecodeDateDay'     ,grib% s1_dwd% day)
          call grib_set (h,'localDecodeDateHour'    ,grib% s1_dwd% hour)
          call grib_set (h,'localDecodeDateMinute'  ,grib% s1_dwd% minute)
       case (2)
          call grib_set (h,'grib2LocalSectionPresent',1)
          call grib_set (h,'localDefinitionNumber'  ,grib% isec1% local_ident)
!         call grib_set_missing (h,'localHostIdentifier')
          call grib_set (h,'localHostIdentifier'    ,0) ! 1=oper.,2=backup,???
          call grib_set (h,'localCreationDateYear'  ,grib% s1_dwd% year + 1900)
          call grib_set (h,'localCreationDateMonth' ,grib% s1_dwd% month)
          call grib_set (h,'localCreationDateDay'   ,grib% s1_dwd% day)
          call grib_set (h,'localCreationDateHour'  ,grib% s1_dwd% hour)
          call grib_set (h,'localCreationDateMinute',grib% s1_dwd% minute)
          call grib_set (h,'localCreationDateSecond',0)
          call grib_set (h,'localNumberOfExperiment',nex)
          call grib_set (h,'localInformationNumber' ,grib% s1_dwd% element_no)
!         vers = grib_get_api_version ()
          vers = 10000 * major() + 100 * minor()    ! currently no patchlevel
          call grib_set (h,'localVersionNumber'     ,vers)
          call grib_set (h,'backgroundProcess'      ,run_type)
       end select
    end if
#endif
  end SUBROUTINE set_dwd_local
!------------------------------------------------------------------------------
  !------------------------------
  ! specify DWD local part of PDB
  !------------------------------
  SUBROUTINE set_dwd_ens (grib, n_ens, i_ens, id_ens)
  TYPE (t_grib1) ,INTENT(inout) :: grib     ! grib buffer
  INTEGER        ,INTENT(in)    :: n_ens    ! ensemble size
  INTEGER        ,INTENT(in)    :: i_ens    ! ensemble member number
  INTEGER        ,INTENT(in)    :: id_ens   ! ensemble id

    integer :: en_id    ! Actual ensemble id
#ifdef GRIB_API
    integer :: h        ! handle
    integer :: der      ! derived forecast ident.
#endif
    logical :: bigens   ! DWD "big" ensemble

    en_id = id_ens; if (id_ens < 0) en_id = d_id_en
    bigens = n_ens > 254

    if (n_ens > 0 .and. i_ens > 0) then
      if (bigens) then
        grib% isec1 % local_ident = 153     ! DWD "big" ensemble
        grib% s1_dwd% local_ident = 153
      else
        grib% isec1 % local_ident = 253     ! DWD "small" ensemble
        grib% s1_dwd% local_ident = 253
      end if
      grib% s1_dwd% ensemble_size = n_ens
      grib% s1_dwd% ensemble_no   = i_ens
      grib% s1_dwd% ensemble_id   = en_id
    else if (n_ens > 0 .and. any (i_ens == [ENS_MEAN,ENS_SPREAD])) then
      if (grib% isec0% edition == 1) then
        grib% isec1 % local_ident = 253     ! DWD GRIB_API:
        grib% s1_dwd% local_ident = 253     ! only 253 for grib1 available
      else if (bigens) then
        grib% isec1 % local_ident = 152
        grib% s1_dwd% local_ident = 152
      else
        grib% isec1 % local_ident = 252
        grib% s1_dwd% local_ident = 252
      end if
      grib% s1_dwd% ensemble_size = n_ens
      grib% s1_dwd% ensemble_no   = -1
      grib% s1_dwd% ensemble_id   = en_id
    else
      grib% isec1 % local_ident   = 254
      grib% s1_dwd% local_ident   = 254
    endif

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       h = grib% handle
       call grib_set (h,'localDefinitionNumber',grib% isec1% local_ident)
       select case (grib% isec0% edition)
       case (1)
         if (grib% isec1% local_ident == 253) then
           call grib_set (h,'localModelMajorVersionNumber',grib% s1_dwd% major_version)
           call grib_set (h,'localModelMinorVersionNumber',grib% s1_dwd% minor_version)
         end if
         if (n_ens > 0 .and. i_ens > 0) then
           call grib_set (h,'localEnsembleIdentification'      ,en_id)
           call grib_set (h,'localNumberOfEnsembleMembers'     ,n_ens)
           call grib_set (h,'localActualNumberOfEnsembleNumber',i_ens)
         else if (n_ens > 0 .and. any (i_ens == [ENS_MEAN,ENS_SPREAD])) then
           call grib_set (h,'localEnsembleIdentification'      ,en_id)
           call grib_set (h,'localNumberOfEnsembleMembers'     ,n_ens)
         end if
       case (2)
         if (n_ens > 0 .and. i_ens > 0) then
           grib% pdtn = GRIB_4_0_FC_ENS
           call grib_set (h,'typeOfProcessedData'            ,5)            ! Tab.1.4
           call grib_set (h,'productDefinitionTemplateNumber',grib% pdtn)   ! Ens.
           call grib_set (h,'typeOfEnsembleForecast'         ,192)          ! Tab.4.6
           call grib_set (h,'localTypeOfEnsembleForecast'    ,en_id)
           if (bigens) then
             call grib_set (h,'localNumberOfForecastsInEnsemble',n_ens)     ! DWD
             call grib_set (h,'localPerturbationNumber'         ,i_ens)
           else
             call grib_set (h,'numberOfForecastsInEnsemble'  ,n_ens)
             call grib_set (h,'perturbationNumber'           ,i_ens)
           end if
         else if (n_ens > 0 .and. any (i_ens == [ENS_MEAN,ENS_SPREAD])) then
           der        = GRIB_4_7_MEAN; if (i_ens == ENS_SPREAD) der = GRIB_4_7_SPREAD
           grib% pdtn = GRIB_4_0_FC_ENS_ALL
           call grib_set (h,'typeOfProcessedData'            ,5)            ! Tab.1.4
           call grib_set (h,'productDefinitionTemplateNumber',grib% pdtn)
           call grib_set (h,'localTypeOfEnsembleForecast'    ,en_id)
           call grib_set (h,'derivedForecast'                ,der)
           if (bigens) then
             call grib_set (h,'localNumberOfForecastsInEnsemble',n_ens)     ! DWD
           else
             call grib_set (h,'numberOfForecastsInEnsemble'     ,n_ens)
           end if
         else
           grib% pdtn = GRIB_4_0_ANA_FC
           call grib_set (h,'productDefinitionTemplateNumber',grib% pdtn)   ! Det.
         end if
       end select
    end if
#endif
  END SUBROUTINE set_dwd_ens
!------------------------------------------------------------------------------
  !--------------------------------
  ! specify ECMWF local part of PDS
  ! (some entries preliminarily hard-coded for testing)
  !--------------------------------
  SUBROUTINE set_ecmwf_ens (grib, n_ens, i_ens, expid)
  TYPE(t_grib1)   ,INTENT(inout)       :: grib  ! grib buffer
  INTEGER         ,INTENT(in)          :: n_ens ! ensemble size (incl. control)
  INTEGER         ,INTENT(in)          :: i_ens ! ensemble member number
  CHARACTER(len=*),INTENT(in),optional :: expid ! Experiment identification

#ifdef GRIB_API
    integer :: h        ! handle
#endif

    grib% isec1 % local_flag    = 1     ! Indicates ECMWF local extensions
    grib% isec1 % local_ident   = 1     ! local GRIB use definition identifier
    grib% s1_ecmwf% local_ident = 1     ! Mars labelling or ensemble fc data
    grib% s1_ecmwf% class       = 2     ! =ECMWF research dep.
!   grib% s1_ecmwf% class       = 106   ! =de (Germany)
    if (n_ens > 1) then
      if (i_ens > 0) then
        grib% s1_ecmwf% type    = 11    ! =perturbed forecast (ensembles)
      else
        grib% s1_ecmwf% type    = 10    ! =control   forecast (ensembles)
      end if
      grib%   s1_ecmwf% stream  = 1035  ! =enfo
      grib%   s1_ecmwf% i_en    = i_ens
      grib%   s1_ecmwf% n_en    = n_ens
    else
      if (d_run == 0) then
        grib% s1_ecmwf% type    = 9     ! =forecast
      else
        grib% s1_ecmwf% type    = 2     ! =analysis
      end if
      grib%   s1_ecmwf% stream  = 1025  ! =oper
      grib%   s1_ecmwf% i_en    = 0
      grib%   s1_ecmwf% n_en    = 1
    end if
    grib% s1_ecmwf% expid       = d_expid
    if (present (expid)) grib% s1_ecmwf% expid = expid

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       h = grib% handle
       call grib_set (h,'localDefinitionNumber',grib% isec1% local_ident)
       if (n_ens > 1) then
         write(0,*) "Warning: ECMWF ensemble not yet supported with GRIB_API"
         CALL finish ('set_ecmwf_ens','GRIB_API not yet supported')
       else
         write(0,*) "Warning: ECMWF not yet fully supported with GRIB_API"
       end if
    end if
#endif
  END SUBROUTINE set_ecmwf_ens
!------------------------------------------------------------------------------
  !-----------------------
  ! specify reference time
  !-----------------------
  SUBROUTINE set_ref_time (grib, ref_time)
  TYPE (t_grib1) ,INTENT(inout) :: grib
  TYPE (t_time)  ,INTENT(in)    :: ref_time

    INTEGER :: yyyy
    INTEGER :: hhmm
#ifdef GRIB_API
    integer :: yyyymmdd
#endif

    hhmm                 = ihhmm (ref_time)
    yyyy                 = iyyyy (ref_time)
    grib% isec1% day     = idd   (ref_time)
    grib% isec1% month   = imm   (ref_time)
    grib% isec1% hour    = hhmm       /100
    grib% isec1% minute  = MOD (hhmm  ,100)
    grib% isec1% year    = MOD (yyyy-1,100) + 1
    grib% isec1% century = (yyyy-grib% isec1% year)/100  + 1
    !---------------------------------------------------------------
    ! by default, set the time range to forecast
    !             set verification time to reference time (analysis)
    !---------------------------------------------------------------
    grib% isec1% time_range = WMO5_FORECAST
    grib% isec1% p1         = 0
    grib% isec1% p2         = 0
    grib% isec1% time_unit  = 0

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       yyyymmdd = yyyy * 10000 + grib% isec1% month * 100 + grib% isec1% day
       if (grib% isec0% edition == 2) then
          ! grib2/tables/4/1.2.table  :: 0=ana,1=fc_start,2=veri_time
          call grib_set (grib% handle,'significanceOfReferenceTime', 0)
       end if
       call grib_set (grib% handle,'dataDate',yyyymmdd)
       call grib_set (grib% handle,'dataTime',hhmm)
       call grib_set (grib% handle,'stepType','instant')
    end if
#endif
  END SUBROUTINE set_ref_time
!------------------------------------------------------------------------------
  !--------------------------------------
  ! specify verification time, time range
  !--------------------------------------
  SUBROUTINE set_ver_time (grib, ver_time, range, unit)
  TYPE (t_grib1) ,INTENT(inout)        :: grib
  TYPE (t_time)  ,INTENT(in)           :: ver_time
  INTEGER        ,INTENT(in) ,OPTIONAL :: range
  INTEGER        ,INTENT(in) ,OPTIONAL :: unit

    INTEGER       :: lrange
    INTEGER       :: lunit
    INTEGER       :: hhmmss
    INTEGER       :: p1, p2
    TYPE (t_time) :: r_time, p1_time
    !----------------------------------------------
    ! set default/initial values of local variables
    !----------------------------------------------
    lrange = WMO5_FORECAST; IF (PRESENT(range)) lrange = range
    p1     = 0
    p2     = 0
    lunit  = 0; IF (PRESENT(unit)) lunit = unit
    IF (grib% isec1% year == 255) THEN
      CALL set_ref_time (grib, ver_time)
      r_time = ver_time
    ELSE
      r_time = ref_time (grib)
    ENDIF
    !----------------------------------
    ! calculate quantities for forecast
    !----------------------------------
    SELECT CASE (lrange)
    CASE (WMO5_FORECAST)
      p1_time = ver_time - r_time
      !-----------------------------
      ! set time unit if not yet set
      !-----------------------------
      hhmmss     = ihhmmss (p1_time)
      IF (.NOT.PRESENT (unit)) THEN
        IF (         hhmmss       ==0) THEN
!         lunit = WMO4_DAY currently always hours !!!
          lunit = WMO4_HOUR
        ELSE IF (MOD(hhmmss,10000)==0) THEN
          lunit = WMO4_HOUR
        ELSE IF (MOD(hhmmss,  100)==0) THEN
          lunit = WMO4_MINUTE
        ELSE
          lunit = WMO4_SECOND
        ENDIF
      ENDIF
      !-------------------------
      ! set p1 in units required
      !-------------------------
      SELECT CASE (lunit)
      CASE (WMO4_DAY)
        p1 = int (days    (p1_time))
      CASE (WMO4_HOUR)
        p1 = int (hours   (p1_time))
      CASE (WMO4_MINUTE)
        p1 = int (minutes (p1_time))
      CASE (WMO4_SECOND)
        p1 = int (seconds (p1_time))
      CASE default
        CALL finish ('set_ver_time','time unit not implemented')
      END SELECT
      if (p1_time == zero_time .and. .not. present(range)) lrange = d_ran
    CASE default
      CALL finish ('set_ver_time','time range not implemented')
    END SELECT
    !--------------------------------------
    ! finally set entries in GRIB section 1
    !--------------------------------------
    grib% isec1% time_range = lrange
    grib% isec1% time_unit  = lunit
    grib% isec1% p1         = p1
    grib% isec1% p2         = p2

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       if (grib% isec0% edition == 2) then
          ! grib2/tables/4/1.2.table  :: 0=ana,1=fc_start,2=veri_time
          if (p1 == 0 .and. p2 == 0) then
             call grib_set (grib% handle,'significanceOfReferenceTime', 0)
          else
             call grib_set (grib% handle,'significanceOfReferenceTime', 1)
          end if
          ! Section 1, grib2/1.4.table :: 0=ana,1=fc,4:5=pert.fc,...
          if (any (grib% isec1 % local_ident == [152,153,252,253]) .and. &
                   grib% s1_dwd% ensemble_size > 0                       ) then
             call grib_set (grib% handle,'typeOfProcessedData', 5)
          else if (p1 == 0 .and. p2 == 0) then
             call grib_set (grib% handle,'typeOfProcessedData', 0)
          else
             call grib_set (grib% handle,'typeOfProcessedData', 1)
          end if
       end if

       select case (lunit)
       case (WMO4_SECOND)
         lunit = 0
       case (DWD4_15_MINUTES)
         call finish('set_ver_time','no GRIB-2 equivalent for DWD4_15_MINUTES')
       end select

       call grib_set   (grib% handle,'stepType' ,'instant')
       if (grib% isec0% edition == 1) then
         call grib_set (grib% handle,'timeRangeIndicator',lrange)
       end if
       call grib_set   (grib% handle,'stepUnits',lunit)
       call grib_set   (grib% handle,'startStep',p1)
       call grib_set   (grib% handle,'endStep'  ,p1+p2)
    end if
#endif
  END SUBROUTINE set_ver_time
!------------------------------------------------------------------------------
  !--------------------------------------------------
  ! specify leveltype, vertical coordinate parameters
  !--------------------------------------------------
  SUBROUTINE set_grid_vert (grib, leveltype, a, b, wlev, nlev, gridnum, uuid)
  TYPE (t_grib1)  ,INTENT(inout)        :: grib
  INTEGER         ,INTENT(in)           :: leveltype
  REAL(wp)        ,INTENT(in) ,OPTIONAL :: a (:)
  REAL(wp)        ,INTENT(in) ,OPTIONAL :: b (:)
  LOGICAL         ,INTENT(in) ,OPTIONAL :: wlev      ! write levels
  INTEGER         ,INTENT(in) ,OPTIONAL :: nlev      ! no.levels
  INTEGER         ,INTENT(in) ,OPTIONAL :: gridnum   ! numberOfVGridUsed
  CHARACTER(len=1),INTENT(in) ,OPTIONAL :: uuid(:)   ! uuidOfVGrid
    INTEGER :: nvcp, na
    LOGICAL :: wl
#ifdef GRIB_API
    integer :: h, h2    ! handles
    integer :: ffs, sfs ! Type of first/second fixed surface
    logical :: lgenv    ! Generalized vertical coord.
    integer :: oldffs   ! Old type of first fixed surface
!   integer :: oldnv    ! Old number of vertical coord. parameters
#endif

    wl = d_wlev; if (present(wlev)) wl = wlev
    nvcp = 0
    IF (wl) THEN
      IF (PRESENT (a)) THEN
        nvcp = nvcp + SIZE (a)
        grib% rsec2% vcp (1:nvcp) = a
      ENDIF
      IF (PRESENT (b)) THEN
        na   = nvcp
        nvcp = nvcp + SIZE (b)
        grib% rsec2% vcp (na+1:nvcp) = b
      ENDIF
    ENDIF
    grib% isec1% level_type = leveltype
    grib% latlon% nvcp      = nvcp
    grib% gauss%  nvcp      = nvcp
    grib% tri%    nvcp      = nvcp
    grib% sph%    nvcp      = nvcp

#ifdef GRIB_API
    lgenv = .false.
    if (grib% handle /= INVALID_HANDLE) then
       h = grib% handle

       !-------------------------------------------------------------
       ! Explicit cloning to enforce/maintain consistency of metadata
       !-------------------------------------------------------------
       call grib_clone   (h ,h2)
       call grib_release (h)
       grib% handle = h2
       h = h2

       if (grib% isec0% edition == 1) then
          call grib_set (h,'indicatorOfTypeOfLevel',leveltype)
       else
          ffs = 255
          sfs = 255
          select case (leveltype)
          case (WMO3_SURFACE)
             ffs = 1
          case (WMO3_ISOBARIC)
             ffs = 100
          case (WMO3_ABOVESUR)
             ffs = 103
!            sfs = 181  ! for t2m_land, rh2m_land
          case (WMO3_BELOWSUR)
             ffs = 106
!            sfs = 106  ! for w_so, w_so_ice
          case (WMO3_SURF_HORZ)
             ffs = WMO3_SURF_HORZ
          case (WMO3_HYBRID)
             ffs = 105
          case (WMO3_HYBRIDB, WMO3_HHYBRID)
             ffs = 105
             sfs = 105
          case (WMO3_GENV)
             if (.not. present (nlev)) &
               call finish ("set_grid_vert","nlev required for GENV")
             if (.not. present (gridnum)) &
               call finish ("set_grid_vert","gridnum required for GENV")
             if (.not. present (uuid)) &
               call finish ("set_grid_vert","uuid required for GENV")
             ffs = 150
             sfs = 150
             lgenv = .true.
             wl    = .false.
          case default
             write(0,*) "leveltype =", leveltype
             call finish ("set_grid_vert","unsupported leveltype")
          end select
          call grib_get (h, 'typeOfFirstFixedSurface', oldffs)
!         print *, "oldffs =", oldffs
!         call grib_get (h, 'NV', oldnv)
!         print *, "NV old =", oldnv
       end if
       if (wl .and. nvcp > 0) then
          if (grib% isec0% edition == 2) then
             call grib_set (h,'typeOfFirstFixedSurface' ,ffs)
             call grib_set (h,'typeOfSecondFixedSurface',sfs)
          end if
          call grib_set (h,'PVPresent',1)               ! Force PVPresent!
          call grib_set (h,'pv'       ,grib% rsec2% vcp(1:nvcp))
       else
          if (grib% isec0% edition == 1) then
             call grib_set (h,'PVPresent',0)            ! Remove PV
             call grib_set (h,'NV'       ,0)
          else
          !-----------------------------------------------------
          ! GRIB_API: must set NV before typeOf*FixedSurface !!!
          !-----------------------------------------------------
             if (lgenv) then
               call grib_set   (h,'genVertHeightCoords'     ,1   ) ! for 1.12.3
               call grib_set   (h,'NV'                      ,6   )
               call grib_set   (h,'typeOfFirstFixedSurface' ,ffs )
               call grib_set   (h,'typeOfSecondFixedSurface',sfs )
!              call grib_set   (h,'genVertHeightCoords'     ,1   )
               call grib_set   (h,'nlev'                    ,nlev)
               call grib_set   (h,'numberOfVGridUsed'       ,gridnum)
#if !defined(__ibm__)
!!! This requires GRIB_API version >= 1.11.0
               call grib_set   (h,'uuidOfVGrid'             ,uuid)
#endif
             else
               call grib_set   (h,'genVertHeightCoords'     ,0   ) ! for 1.12.3
               call grib_set   (h,'typeOfFirstFixedSurface' ,ffs )
               call grib_set   (h,'typeOfSecondFixedSurface',sfs )
               if (oldffs == 150) then
                 call grib_set (h,'NV'                      ,0   ) ! Beware!
               end if

               !-------------------------------------------------------------
               ! Explicit cloning to enforce/maintain consistency of metadata
               !-------------------------------------------------------------
               call grib_clone   (h ,h2)
               call grib_release (h)
               grib% handle = h2
               h = h2

               if (oldffs == 105) then
                 call grib_set (h,'PVPresent'               ,0   ) ! Remove PV
                 call grib_set (h,'NV'                      ,0   )
               end if
             end if
          end if
       end if
    end if
#endif
  END SUBROUTINE set_grid_vert
!------------------------------------------------------------------------------
  !---------------------
  ! specify lat/lon grid
  !---------------------
  SUBROUTINE set_grid_latlon (grib, ni, nj, dlatf, dlonf, di, dj, &
                              scan, dlatr, dlonr)
  TYPE (t_grib1) ,INTENT(inout)        :: grib
  INTEGER        ,INTENT(in)           :: ni    ! Number of points in i
  INTEGER        ,INTENT(in)           :: nj    ! Number of points in j
  REAL(wp)       ,INTENT(in)           :: dlatf ! first latitude       [degree]
  REAL(wp)       ,INTENT(in)           :: dlonf ! first longitude      [degree]
  REAL(wp)       ,INTENT(in)           :: di    ! i direction increment[degree]
  REAL(wp)       ,INTENT(in)           :: dj    ! j direction increment[degree]
  INTEGER        ,INTENT(in) ,OPTIONAL :: scan  ! Scanning mode flags (table 8)
  REAL(wp)       ,INTENT(in) ,OPTIONAL :: dlatr ! latitude  of rotated S-pole
  REAL(wp)       ,INTENT(in) ,OPTIONAL :: dlonr ! longitute of rotated S-pole
    INTEGER :: sc
#ifdef GRIB_API
    integer  :: h          ! handle
    real(wp) :: lon1, lon2 ! Longitudes modulo 360
#endif
    integer  :: rfac       ! resolution factor for grib1/grib2

    if (grib% isec0% edition == 1) then
       rfac = 1000         ! resolution: 1/1000°
    else
       rfac = 1000000      ! resolution: 1/1000000°
    end if
    !--------------------------------
    ! set scan mode (cf. WMO table 8)
    !--------------------------------
    sc = 0
    IF (di < 0) sc = ior (sc, WMO8_I_NEGATIVE)  ! 128
    IF (dj > 0) sc = ior (sc, WMO8_J_POSITIVE)  !  64
    IF(PRESENT(scan)) sc = scan
    !--------------------------------------------
    ! set parameters in GRIB section 2 (defaults)
    !--------------------------------------------
!   grib% latlon% nvcp             ! Number of vert.coord.param. set elsewhere
    grib% latlon% lat_rot     =  0 ! Latitude of southern pole of rotation.
    grib% latlon% lon_rot     =  0 ! Longitude of southern pole of rotation.
    grib% latlon% lat_strech  =  0 ! Latitude of the pole of stretching.
    grib% latlon% lon_strech  =  0 ! Longitude of the pole of stretching.
    grib% latlon% reduced     =  0 ! 0: Regular grid, 1: reduced grid.
    grib% latlon% earth       =  0 ! 0: spherical r=6367.47km
    grib% latlon% components  =  0 ! 0: u and v components in E/N
    grib% latlon% reserved    =  0
    !---------------------------------
    ! set parameters in GRIB section 2
    !---------------------------------
    grib% isec2   (1)         = WMO6_LATLON
    grib% latlon% repr        = WMO6_LATLON  ! Data representation type
    grib% latlon% ni          = ni           ! Number of points along parallel.
    grib% latlon% nj          = nj           ! Number of points along meridian.
    grib% latlon% lat_first   = nint (dlatf * rfac) ! Latitude  1st grid point
    grib% latlon% lon_first   = nint (dlonf * rfac) ! Longitude 1st grid point
    grib% latlon% increments  = 128          ! Direction increments given.
    grib% latlon% di          = nint (di    * rfac) ! i direction increment
    grib% latlon% dj          = nint (dj    * rfac) ! j direction increment
    grib% latlon% lat_last    = nint (((nj-1)*dj + dlatf) * rfac) ! last Latit.
    grib% latlon% lon_last    = nint (((ni-1)*di + dlonf) * rfac) ! last Long.
    grib% latlon% scan_mode   = sc           ! Scanning mode flags(WMO table 8)
    !-------------------------------------------------------
    ! don't give di, dj if integer representation is inexact
    !-------------------------------------------------------
    if (grib% latlon% dj * (nj-1) /=                          &
        grib% latlon% lat_last - grib% latlon% lat_first .or. &
        grib% latlon% di * (ni-1) /=                          &
        grib% latlon% lon_last - grib% latlon% lon_first      ) then
      grib% latlon% increments  = 0          ! Direction increments not given.
      grib% latlon% di          = 65535      ! Set all bits to ones
      grib% latlon% dj          = 65535      ! (required by GRIB1)
    endif
    !-------------------------------------------
    ! optionally set parameters for rotated grid
    !-------------------------------------------
    if (present (dlatr)) then
      grib% isec2   (1)       = WMO6_ROTLL
      grib% latlon% repr      = WMO6_ROTLL          ! Data representation type
      grib% latlon% lat_rot   = nint (dlatr * rfac) ! Latitude  S-pole
      grib% latlon% lon_rot   = nint (dlonr * rfac) ! Longitude S-pole
      !---------------------
      ! u,v relative to grid
      !---------------------
      grib% latlon% increments= grib% latlon% increments + 8
    endif
    !--------------------------------------
    ! derive number of gridpoints per level
    !--------------------------------------
    grib% isec4% n_data   = ni * nj
    grib% isec4% non_miss = ni * nj

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       h = grib% handle
       call grib_set (h,'latitudeOfFirstGridPointInDegrees' ,dlatf)
       call grib_set (h,'latitudeOfLastGridPointInDegrees'  ,dlatf+(nj-1)*dj)
       lon1 = dlonf
       lon2 = dlonf + (ni-1)*di
       if (grib% isec0% edition == 2 .and. lon1 < 0) lon1 = lon1 + 360._wp
       if (grib% isec0% edition == 2 .and. lon2 < 0) lon2 = lon2 + 360._wp
       call grib_set (h,'longitudeOfFirstGridPointInDegrees', lon1)
       call grib_set (h,'longitudeOfLastGridPointInDegrees' , lon2)
       if (iand (grib% latlon% increments, 128) == 0) then
          call grib_set (h,'ijDirectionIncrementGiven'    ,0)
       else
          call grib_set (h,'iDirectionIncrementInDegrees' ,di)
          call grib_set (h,'jDirectionIncrementInDegrees' ,dj)
       end if
       select case (grib% latlon% repr)
       case (WMO6_LATLON)
          call grib_set (h,'gridType','regular_ll')
       case (WMO6_ROTLL)
          call grib_set (h,'gridType','rotated_ll')
          call grib_set (h,'uvRelativeToGrid'                ,1)
          call grib_set (h,'latitudeOfSouthernPoleInDegrees' ,dlatr)
          call grib_set (h,'longitudeOfSouthernPoleInDegrees',dlonr)
          call grib_set (h,'angleOfRotationInDegrees'        ,0)
       end select
       call grib_set (h,'Ni',ni)
       call grib_set (h,'Nj',nj)
       call grib_set (h,'iScansNegatively', &
                      merge (1, 0, iand (sc, WMO8_I_NEGATIVE) /= 0))
       call grib_set (h,'jScansPositively', &
                      merge (1, 0, iand (sc, WMO8_J_POSITIVE) /= 0))
    end if
#endif
  END SUBROUTINE set_grid_latlon
!------------------------------------------------------------------------------
  !----------------------
  ! specify Gaussian grid
  !----------------------
  SUBROUTINE set_grid_gauss (grib, ni, nj, dlatf, dlonf, di, scan)
  TYPE (t_grib1) ,INTENT(inout)        :: grib
  INTEGER        ,INTENT(in) ,OPTIONAL :: ni    ! Number of points in i
  INTEGER        ,INTENT(in)           :: nj    ! Number of points in j
  REAL(wp)       ,INTENT(in)           :: dlatf ! first latitude       [degree]
  REAL(wp)       ,INTENT(in) ,OPTIONAL :: dlonf ! first longitude      [degree]
  REAL(wp)       ,INTENT(in) ,OPTIONAL :: di    ! i direction increment[degree]
  INTEGER        ,INTENT(in) ,OPTIONAL :: scan  ! Scanning mode flags (table 8)
    INTEGER  :: sc, lni
    REAL(wp) :: lonf, ldi
#ifdef GRIB_API
    integer  :: h        ! handle
#endif
    !------------------------
    ! set optional parameters
    !------------------------
    lni  = nj * 2      ;if (present(ni))    lni  = ni
    lonf = 0           ;if (present(dlonf)) lonf = dlonf
    ldi  = 360._wp/lni ;if (present(di))    ldi  = di
    !--------------------------------
    ! set scan mode (cf. WMO table 8)
    !--------------------------------
    sc = 0                                      ! default: -j
    IF (ldi < 0) sc = ior (sc, WMO8_I_NEGATIVE) ! 128
!   IF (dj > 0) sc = ior (sc, WMO8_J_POSITIVE)  !  64
    IF(PRESENT(scan)) sc = scan
    !---------------------------------
    ! set parameters in GRIB section 2
    !---------------------------------
    grib% isec2   (1)        = WMO6_GAUSSIAN
    grib% gauss% repr        = WMO6_GAUSSIAN ! Data representation type (192)
    grib% gauss% ni          = lni           ! Number of points along parallel.
    grib% gauss% nj          = nj            ! Number of points along meridian.
    grib% gauss% lat_first   = nint (dlatf * 1000) ! Latitude  1st grid point
    grib% gauss% lon_first   = nint (lonf  * 1000) ! Longitude 1st grid point
    grib% gauss% increments  = 128           ! Direction increments given.
    grib% gauss% di          = nint (ldi   * 1000) ! i direction increment
    grib% gauss% nglh        = nj / 2        ! parallels betw. pole and equator
    grib% gauss% lat_last    = - grib% gauss% lat_first
    grib% gauss% lon_last    = grib% gauss% lon_first &
                             + nint (ldi * (ni-1) * 1000)
    grib% gauss% scan_mode   = sc            ! Scanning mode flags(WMO table 8)
    !-----------------------------------------------
    ! remaining parameters are set to default values
    !-----------------------------------------------
!   grib% gauss% nvcp             ! Number of vert.coord.param. set elsewhere
    grib% gauss% lat_rot     =  0 ! Latitude of southern pole of rotation.
    grib% gauss% lon_rot     =  0 ! Longitude of southern pole of rotation.
    grib% gauss% lat_strech  =  0 ! Latitude of the pole of stretching.
    grib% gauss% lon_strech  =  0 ! Longitude of the pole of stretching.
    grib% gauss% reduced     =  0 ! 0: Regular grid, 1: reduced grid.
    grib% gauss% earth       =  0 ! 0: spherical r=6367.47km
    grib% gauss% components  =  0 ! 0: u and v components in E/N
    grib% gauss% reserved    =  0
    !--------------------------------------
    ! derive number of gridpoints per level
    !--------------------------------------
    grib% isec4% n_data   = ni * nj
    grib% isec4% non_miss = ni * nj

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       h = grib% handle
       call grib_set (h,'latitudeOfFirstGridPointInDegrees' , dlatf)
       call grib_set (h,'latitudeOfLastGridPointInDegrees'  ,-dlatf)
       call grib_set (h,'longitudeOfFirstGridPointInDegrees',  lonf)
       call grib_set (h,'longitudeOfLastGridPointInDegrees' ,  lonf+(ni-1)*ldi)
       call grib_set (h,'iDirectionIncrementInDegrees'      ,ldi)
       call grib_set (h,'gridType','regular_gg')
       call grib_set (h,'Ni',lni)
       call grib_set (h,'Nj', nj)
       call grib_set (h,'N' , nj/2)
!      call grib_set (h,'global',1) ! Beware: this flips first/last latitude!
       call grib_set (h,'iScansNegatively', &
                      merge (1, 0, iand (sc, WMO8_I_NEGATIVE) /= 0))
       call grib_set (h,'jScansPositively', &
                      merge (1, 0, iand (sc, WMO8_J_POSITIVE) /= 0))
    end if
#endif
  END SUBROUTINE set_grid_gauss
!------------------------------------------------------------------------------
  !------------------------
  ! specify triangular grid
  !------------------------
  SUBROUTINE set_grid_tri (grib, ni, nd)
  TYPE (t_grib1) ,INTENT(inout)        :: grib
  INTEGER        ,INTENT(in)           :: ni   ! Number of subdivisions
  INTEGER        ,INTENT(in)           :: nd   ! Number of diamonds

    INTEGER :: ni2, ni3
#ifdef GRIB_API
    integer :: h        ! handle
    integer :: stat     ! error status for handling of key renaming
#endif

    CALL factorize_ni (ni, ni2, ni3)
    !---------------------------------
    ! set parameters in GRIB section 2
    !---------------------------------
    grib% isec2 (1)        = DWD6_ICOSAHEDRON
    grib% tri% repr        = DWD6_ICOSAHEDRON ! Data representation type (192).
    grib% tri% ni          = ni    ! Number of triangular subdivisions.
    grib% tri% ni2         = ni2   ! Times of factor 2 in factorisation of NI.
    grib% tri% ni3         = ni3   ! Times of factor 3 in factorisation of NI.
    grib% tri% nd          = nd    ! Number of diamonds.
    !-----------------------------------------------
    ! remaining parameters are set to default values
    !-----------------------------------------------
    grib% tri% orient      =   128 ! Flag for orientation of diamonds.
    grib% tri% lat_pole    = 90000 ! Latitude of pole point.
    grib% tri% lon_pole    =     0 ! Longitude of pole point.
    grib% tri% lon_dia1    =     0 ! Longitude of the first diamond.
    grib% tri% scan_mode   =     0 ! Flag for storage sequence.
    grib% tri% reserved_11 =     0 !
!   grib% tri% nvcp                ! Number of vert.coord.param. set elsewhere
    grib% tri% reserved    =     0 !
    !-----------------------------------------------------
    ! take care that this derived type is shared with ICON
    !-----------------------------------------------------
    grib% tri% npts        =     0 !
    grib% tri% grid_num    =     0 !
    grib% tri% grid_ref    =     0 !
    grib% tri% uuid        = achar (0)
    !--------------------------------------
    ! derive number of gridpoints per level
    !--------------------------------------
    grib% isec4% n_data    = (ni+1)**2 * nd
    grib% isec4% non_miss  = grib% isec4% n_data

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       h = grib% handle
       call grib_set (h,'gridType','triangular_grid')
       call grib_set (h,'n2',ni2)
       call grib_set (h,'n3',ni3)
       call grib_set (h,'nd',nd)
       call grib_set (h,'Ni',ni)
       call grib_set (h,'numberingOrderOfDiamonds' ,grib% tri% orient)
       call grib_set (h,'scanningModeForOneDiamond',grib% tri% scan_mode)
       select case (grib% isec0% edition)
       case (1)
          call grib_set (h,'latitudeOfIcosahedronPoleInDegrees'        ,90)
          call grib_set (h,'longitudeOfIcosahedronPoleInDegrees'       ,0)
          call grib_set (h,'longitudeOfFirstDiamondCenterLineInDegrees',0,stat)
          if (stat /= 0) &
          call grib_set (h,'longitudeOfFirstDiamondCentreLineInDegrees',0)
       case (2)
          if (grib% isec1% center /= WMO0_ECMWF) then   ! Bug in ECMWF tables
           call grib_set(h,'latitudeOfThePolePointInDegrees'  ,90)
          end if
          call grib_set (h,'longitudeOfThePolePointInDegrees' ,0)
          call grib_set (h,'longitudeOfFirstDiamondCenterLine',0,stat)
          if (stat /= 0) &
          call grib_set (h,'longitudeOfFirstDiamondCentreLine',0)
          ! grib2/tables/4/3.8.table :: 0 = triangle vertices
          call grib_set (h,'gridPointPosition',0)
          call grib_set (h,'totalNumberOfGridPoints'          ,(ni+1)**2 * nd)
       end select
    end if
#endif
  END SUBROUTINE set_grid_tri
!------------------------------------------------------------------------------
  !-----------------------------
  ! specify ICON triangular grid
  !-----------------------------
  SUBROUTINE set_grid_icon (grib, npts, grid_num, uuid)
  TYPE (t_grib1)  ,INTENT(inout)        :: grib
  INTEGER         ,INTENT(in)           :: npts     ! numberOfDataPoints
  INTEGER         ,INTENT(in), optional :: grid_num ! numberOfGridUsed
  CHARACTER(len=1),INTENT(in), optional :: uuid(:)  ! Unique ID of hor. grid

    INTEGER :: ni, ni2, ni3, ierr
#ifdef GRIB_API
    integer :: h        ! handle
#endif

    if (grib% isec0% edition /= 2) then
       call finish ("set_grid_icon","ICON requires GRIB edition 2!")
    end if
    ni = nint (sqrt (npts / 20._wp))
    CALL factorize_ni (ni, ni2, ni3, ierr)
    !---------------------------------
    ! set parameters in GRIB section 2
    !---------------------------------
    grib% isec2 (1)        = DWD6_ICON
    grib% tri% repr        = DWD6_ICON ! Data representation type (193).
    grib% tri% ni          = ni    ! Number of triangular subdivisions.
    grib% tri% ni2         = ni2   ! Times of factor 2 in factorisation of NI.
    grib% tri% ni3         = ni3   ! Times of factor 3 in factorisation of NI.
    grib% tri% nd          = 1     ! No diamonds.
    if (present (uuid)) then
       grib% tri% uuid     = uuid  ! ICON unique horizontal grid ID.
    else
       grib% tri% uuid     = achar (0)
    end if
    !-----------------------------------------------
    ! remaining parameters are set to default values
    !-----------------------------------------------
    grib% tri% orient      =     0 ! Flag for orientation of diamonds.
    grib% tri% lat_pole    =     0 ! Latitude of pole point.
    grib% tri% lon_pole    =     0 ! Longitude of pole point.
    grib% tri% lon_dia1    =     0 ! Longitude of the first diamond.
    grib% tri% scan_mode   =     0 ! Flag for storage sequence.
    grib% tri% reserved_11 =     0 !
!   grib% tri% nvcp                ! Number of vert.coord.param. set elsewhere
    grib% tri% reserved    =     0 !
    grib% tri% npts        = npts  ! numberOfDataPoints
    grib% tri% grid_num    =     0 ! numberOfGridUsed
    if (present (grid_num)) &
      grib% tri% grid_num  = grid_num
    grib% tri% grid_ref    =     1 ! numberOfGridInReference (1=Triangle centers)
    !--------------------------------------
    ! derive number of gridpoints per level
    !--------------------------------------
    grib% isec4% n_data    = npts
    grib% isec4% non_miss  = npts

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       h = grib% handle
       ! grib2/tables/4/3.1.table
!      call grib_set (h,'gridType','unstructured_grid')     ! GRIB_API 1.11+
       call grib_set (h,'gridDefinitionTemplateNumber',101)
       ! grib2/tables/4/3.2.table :: R = 6,371,229.0 m
       call grib_set (h,'shapeOfTheEarth'             ,6)
       call grib_set (h,'numberOfGridUsed'            ,grib% tri% grid_num)
       call grib_set (h,'numberOfGridInReference'     ,grib% tri% grid_ref)
#if !defined(__ibm__)
!!! This requires GRIB_API version >= 1.11.0
       call grib_set (h,'uuidOfHGrid'                 ,grib% tri% uuid)
#endif
       call grib_set (h,'numberOfDataPoints'          ,grib% isec4% n_data)
    end if
#endif
  END SUBROUTINE set_grid_icon
!------------------------------------------------------------------------------
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE factorize_ni (kni, kni2, kni3, ierr)
  INTEGER, INTENT(in)            :: kni
  INTEGER, INTENT(out)           :: kni2, kni3
  INTEGER, INTENT(out), optional :: ierr

    INTEGER  :: mx

    mx    = kni
    kni2  = 0
    kni3  = 0

    if (present (ierr)) ierr = 0
    DO WHILE (mx > 1)
      IF (MOD(mx,2) == 0) THEN
        kni2  = kni2+1
        mx    = mx/2
      ELSE IF (MOD(mx,3) == 0) THEN
        kni3  = kni3+1
        mx    = mx/3
      ELSE if (present (ierr)) then
        ierr = mx
        exit
      ELSE
        CALL finish('factorize_ni','error in factorisation')
      ENDIF
    ENDDO

  END SUBROUTINE factorize_ni
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!------------------------------------------------------------------------------
  !--------------------------
  ! specify code, table, bits
  !--------------------------
  SUBROUTINE set_code (grib, code, table, dis, cat, num, element_no,ctyp, bits)
  TYPE (t_grib1)    ,INTENT(inout) :: grib
  INTEGER           ,INTENT(in)    :: code       ! grib code number
  INTEGER           ,INTENT(in)    :: table      ! grib table number
  INTEGER ,OPTIONAL ,INTENT(in)    :: dis        ! GRIB2 parameter discipline
  INTEGER ,OPTIONAL ,INTENT(in)    :: cat        ! GRIB2 parameter category
  INTEGER ,OPTIONAL ,INTENT(in)    :: num        ! GRIB2 parameter number
  INTEGER ,OPTIONAL ,INTENT(in)    :: element_no ! Zusatzelementnummer
  INTEGER ,OPTIONAL ,INTENT(in)    :: cTyp       ! constituentType
  INTEGER ,OPTIONAL ,INTENT(in)    :: bits       ! number of bits used for encoding

#ifdef GRIB_API
    integer :: h        ! handle
#endif

    integer :: pdtn ! temporary for Product Definition Template Number

    grib% isec1%  code       = code
    grib% isec1%  table      = table
    grib% s1_dwd% element_no = 0
    grib% ctyp               = -1
    IF (PRESENT(element_no)) grib% s1_dwd% element_no = element_no
    IF (PRESENT(dis))        grib% dis                = dis
    IF (PRESENT(cat))        grib% cat                = cat
    IF (PRESENT(num))        grib% num                = num
    IF (PRESENT(ctyp))       grib% ctyp               = ctyp
    IF (PRESENT(bits)) then
                             grib% isec4% bits        = d_b
       if (bits > 0)         grib% isec4% bits        = bits
    end IF
#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       h = grib% handle
       select case (grib% isec0% edition)
       case (1)
          call grib_set (h,'table2Version'       ,table)
          call grib_set (h,'indicatorOfParameter',code)
          select case (grib% isec1% center)
          case (WMO0_DWD,WMO0_COSMO,WMO0_MSWISS,WMO0_COMET)
            if (grib% s1_dwd% local_ident > 0) &
              call grib_set (h,'localElementNumber',grib% s1_dwd% element_no)
          end select
       case (2)
          call grib_set   (h,'discipline'       ,grib% dis )
          call grib_set   (h,'parameterCategory',grib% cat )
          call grib_set   (h,'parameterNumber'  ,grib% num )
          if (grib% ctyp > 0) then
            !-------------------------------------------------------------
            ! for constituents: set proper productDefinitionTemplateNumber
            !                         and constituentType
            !-------------------------------------------------------------
            select case (grib% pdtn)
            case (GRIB_4_0_ANA_FC)
              call grib_set (h,'productDefinitionTemplateNumber',GRIB_4_0_ANA_FC_CHEM)
            case (GRIB_4_0_FC_ENS)
              call grib_set (h,'productDefinitionTemplateNumber',GRIB_4_0_CHEM_ENS)
            case default
              write(0,*) ' default productDefinitionTemplateNumber : ',grib% pdtn
              call finish('set_code','cannot handle productDefinitionTemplateNumber for constituents')
            end select
            call grib_set (h,'constituentType'  ,grib% ctyp)
          else
            !-------------------------------------------------------
            ! no constituents: reset productDefinitionTemplateNumber
            !-------------------------------------------------------
            call grib_get (h,'productDefinitionTemplateNumber',pdtn)
            if (grib% pdtn /= pdtn) &
              call grib_set (h,'productDefinitionTemplateNumber',grib% pdtn)
          endif
          select case (grib% isec1% center)
          case (WMO0_DWD,WMO0_COSMO,WMO0_MSWISS,WMO0_COMET)
            if (grib% s1_dwd% local_ident > 0) &
              call grib_set (h,'localInformationNumber',grib% s1_dwd% element_no)
          end select
       end select
    end if
#endif
  END SUBROUTINE set_code
!------------------------------------------------------------------------------
  SUBROUTINE set_level1 (grib, level, level_b)
  !----------------------------------------
  ! specify level (using GRIB1 conventions)
  !----------------------------------------
  TYPE (t_grib1) ,INTENT(inout)        :: grib
  INTEGER        ,INTENT(in) ,OPTIONAL :: level
  INTEGER        ,INTENT(in) ,OPTIONAL :: level_b

#ifdef GRIB_API
    integer :: h        ! handle
#endif

    if (grib% isec0% edition == 2) then
       CALL finish ('set_level1','must not be called for GRIB edition 2!')
    end if

    SELECT CASE (grib% isec1% level_type)
    !-----------------------------
    ! levels need not be specified
    !-----------------------------
    CASE (WMO3_SURFACE)
      grib% isec1% level_st = 0
      grib% isec1% level_b  = 0
    !--------------------------
    ! single level is specified
    !--------------------------
    CASE (WMO3_HYBRID, WMO3_ISOBARIC, WMO3_ABOVESUR, WMO3_BELOWSUR)
      IF(.NOT.PRESENT(level)) &
        CALL finish ('set_level','level must be present')
      grib% isec1% level_st = level
      grib% isec1% level_b  = 0
    !------------------------------------
    ! top and bottom levels are specified
    !------------------------------------
    CASE (WMO3_HYBRIDB)
      IF(.NOT.PRESENT(level)) &
        CALL finish ('set_level','level must be present')
      grib% isec1% level_st = level
      IF(PRESENT(level_b)) THEN
        grib% isec1% level_b  = level_b
      ELSE
        grib% isec1% level_b  = level + 1
      ENDIF
    CASE default
      WRITE(0,*) 'set_level: level type not implemented:',grib% isec1% level_type
      CALL finish ('set_level','level type not implemented')
    END SELECT

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       h = grib% handle
       call grib_set (h,'indicatorOfTypeOfLevel' ,grib% isec1% level_type)
       SELECT CASE (grib% isec1% level_type)
       CASE (WMO3_SURFACE, WMO3_ISOBARIC, WMO3_ABOVESUR, WMO3_BELOWSUR)
          call grib_set (h,'level' ,grib% isec1% level_st)
       CASE (WMO3_HYBRID)
          call grib_set (h,'level' ,grib% isec1% level_st)
       CASE (WMO3_HYBRIDB)
          call grib_set (h,'topLevel'   ,grib% isec1% level_st)
          call grib_set (h,'bottomLevel',grib% isec1% level_b)
       case default
          write(0,*) "Level type =", grib% isec1% level_type
          CALL finish ('set_level1','unsupported level type')
       END SELECT
    end if
#endif
  END SUBROUTINE set_level1
!------------------------------------------------------------------------------
  SUBROUTINE set_level2 (grib, level, levtyp, level2, levtyp2)
    !----------------------------------------
    ! specify level (using GRIB2 conventions)
    !----------------------------------------
    TYPE(t_grib1) ,INTENT(inout)        :: grib
    real(wp)      ,INTENT(in)           :: level     ! Level value 1
    integer       ,INTENT(in), OPTIONAL :: levtyp    ! Type of 1st surf.
    real(wp)      ,INTENT(in), OPTIONAL :: level2    ! Level value 2
    integer       ,INTENT(in), OPTIONAL :: levtyp2   ! Type of 2nd surf.

    integer :: scf1, scv1   ! Scale factor, scaled value, level 1
    integer :: scf2, scv2   ! Scale factor, scaled value, level 2
    integer :: ltyp         ! Level type
    integer :: ltyp1        ! Level type (1st surface, for grib2)
    integer :: ltyp2        ! Level type (2nd surface)
    integer :: ffs          ! Type of 1st surface
    integer :: sfs          ! Type of 2nd surface
#ifdef GRIB_API
    integer :: h            ! handle
#endif

    ltyp  = WMO3_SURFACE; if (present (levtyp))  ltyp  = levtyp
    ltyp2 = 255;          if (present (levtyp2)) ltyp2 = levtyp2

    !-----------------------------------------------
    ! Special handling for generalized vertical grid
    !-----------------------------------------------
    ltyp1 = ltyp
    if (grib% isec1% level_type == WMO3_GENV) then
                                  ltyp1 = WMO3_GENV
       if (ltyp2 == WMO3_HYBRIDB) ltyp2 = WMO3_GENV
    end if

    !-------------------------
    ! Emulate GRIB1 convention
    !-------------------------
    if (grib% isec0% edition == 1) then
       grib% isec1% level_type = ltyp
       select case (ltyp)
       case (WMO3_SURFACE)
          call set_level1 (grib, 0)
       case (WMO3_ABOVESUR)
          call set_level1 (grib, nint (level))           ! nearest m
       case (WMO3_BELOWSUR)
          call set_level1 (grib, nint (level * 100._wp)) ! m  -> cm
       case (WMO3_ISOBARIC)
          call set_level1 (grib, nint (level / 100._wp)) ! Pa -> hPa
       case (WMO3_HYBRID, WMO3_HYBRIDB)
          if (present (level2)) then
             call set_level1 (grib, nint (level), nint (level2))
          else
             call set_level1 (grib, nint (level))
          end if
       end select
       return
    end if

!print *,"Enter set_level2",grib% isec0% edition,grib% isec1% level_type,level

    scf1 = 255; scv1 = -1; scf2 = 255; scv2 = -1
    if (ltyp /= 1)        call normalize_value (level,  scf1, scv1)
    if (present (level2)) call normalize_value (level2, scf2, scv2)

    ffs = levtyp_1to2 (ltyp1)
    sfs = levtyp_1to2 (ltyp2)

    select case (ltyp)
    case (WMO3_SURFACE)
!      scf1 = 0         ! Force level to 0 (default: 255 = missing)
!      scv1 = 0         ! Scaled value: -1 = missing!
    case (WMO3_ISOBARIC)
    case (WMO3_SEALEVEL)
    case (WMO3_HEIGHT)
    case (WMO3_ABOVESUR)
    case (WMO3_BELOWSUR)
    case (WMO3_LAKE_BOT,WMO3_SEDIM_BOT,WMO3_MIX_LAYER)
    case (WMO3_SURF_HORZ)
    case (WMO3_HYBRID)
    case (WMO3_HYBRIDB) !,WMO3_HHYBRID,WMO3_GENV
       if (ltyp2 == 255) then
          ltyp2 = ltyp1
          sfs   = ffs
       end if
    case default
       write(0,*) "leveltype =", ltyp
       call finish ("set_level2","leveltype not implemented")
    end select

#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) then
       h = grib% handle
       call grib_set (h,'typeOfFirstFixedSurface' ,ffs)
       call grib_set (h,'typeOfSecondFixedSurface',sfs)
       select case (ltyp)
       case (WMO3_SURFACE,WMO3_LAKE_BOT,WMO3_SEDIM_BOT,WMO3_MIX_LAYER)
          call grib_set_missing (h,'scaleFactorOfFirstFixedSurface')
          call grib_set_missing (h,'scaledValueOfFirstFixedSurface')
          call grib_set_missing (h,'scaleFactorOfSecondFixedSurface')
          call grib_set_missing (h,'scaledValueOfSecondFixedSurface')
       case (WMO3_SEALEVEL)
          call grib_set         (h,'scaleFactorOfFirstFixedSurface' ,0)
          call grib_set         (h,'scaledValueOfFirstFixedSurface' ,0)
          call grib_set_missing (h,'scaleFactorOfSecondFixedSurface')
          call grib_set_missing (h,'scaledValueOfSecondFixedSurface')
       case (WMO3_ISOBARIC,WMO3_HEIGHT)
          call grib_set         (h,'scaleFactorOfFirstFixedSurface' ,scf1)
          call grib_set         (h,'scaledValueOfFirstFixedSurface' ,scv1)
          call grib_set_missing (h,'scaleFactorOfSecondFixedSurface')
          call grib_set_missing (h,'scaledValueOfSecondFixedSurface')
       case (WMO3_ABOVESUR,WMO3_BELOWSUR)
          call grib_set         (h,'scaleFactorOfFirstFixedSurface' ,scf1)
          call grib_set         (h,'scaledValueOfFirstFixedSurface' ,scv1)
          if (scv2 >= 0) then ! Layer below surface
           call grib_set        (h,'scaleFactorOfSecondFixedSurface',scf2)
           call grib_set        (h,'scaledValueOfSecondFixedSurface',scv2)
          else                ! Level below surface
           call grib_set_missing(h,'scaleFactorOfSecondFixedSurface')
           call grib_set_missing(h,'scaledValueOfSecondFixedSurface')
          end if
       case (WMO3_HYBRID,WMO3_HYBRIDB,WMO3_HHYBRID,WMO3_GENV)
          select case (ltyp2)
          case (255)
             call grib_set (h,'level'      ,scv1)
          case (WMO3_SEALEVEL)  ! Beware of special case of HHL
             call grib_set (h,'topLevel'   ,scv1)
             call grib_set (h,'bottomLevel',0   )
          case default
             call grib_set (h,'topLevel'   ,scv1)
             call grib_set (h,'bottomLevel',scv1 + 1)
          end select
       case default
          write(0,*) "leveltype =", ltyp
          call finish ("set_level2","leveltype not implemented for GRIB_API")
       end select
    end if
#endif
  END SUBROUTINE set_level2
!------------------------------------------------------------------------------
  integer function levtyp_1to2 (levtyp, strict) result (levtyp2)
    integer, intent(in)           :: levtyp
    logical, intent(in), optional :: strict
    !---------------------------------------------------------
    ! Convert from WMO grib1, table 3 to grib2, code table 4.5
    !---------------------------------------------------------
    logical :: lstrict
    select case (levtyp)
    case (WMO3_SURFACE)
       levtyp2 = 1
    case (WMO3_ISOBARIC)
       levtyp2 = 100
    case (WMO3_SEALEVEL)
       levtyp2 = 101
    case (WMO3_HEIGHT)
       levtyp2 = 102
    case (WMO3_ABOVESUR)
       levtyp2 = 103
    case (WMO3_HYBRID,WMO3_HYBRIDB)
       levtyp2 = 105
    case (WMO3_BELOWSUR)
       levtyp2 = 106
    case (WMO3_HHYBRID)
       levtyp2 = 105 ! Should be 118, but need backwards compatibility
    case (WMO3_GENV)
       levtyp2 = 150
    case (WMO3_LAKE_BOT,WMO3_SEDIM_BOT,WMO3_MIX_LAYER)
       levtyp2 = levtyp
    case (WMO3_TILE_LAND)
       levtyp2 = levtyp
    case (WMO3_SURF_HORZ)
       levtyp2 = levtyp
    case (255)
       levtyp2 = 255 ! Keep missing value
    case default
       lstrict = .false.; if (present (strict)) lstrict = strict
       if (lstrict) then
          write(0,*) "levtyp=", levtyp
          call finish('levtyp_1to2','invalid level type')
       else
          levtyp2 = 255 ! Missing value
       end if
    end select
  end function levtyp_1to2
!------------------------------------------------------------------------------
  pure subroutine normalize_value (x, scale, value)
    real(wp), intent(in)  :: x      ! Value to normalize
    integer,  intent(out) :: scale  ! Decimal scale
    integer,  intent(out) :: value  ! Scaled value

    real(wp)            :: val, delta
    real(wp), parameter :: eps   = 1.e-5_wp  ! Threshold for approximation
    integer,  parameter :: maxit = 9         ! Max.scale for safe termination

    val   = x
    scale = 0
    do
       value = nint (val)             ! Approximate scaled value
       delta = abs  (val - value)
       if (delta <  eps  ) exit     ! Is it good enough?
       if (scale >= maxit) exit
       scale = scale + 1
       val   = val * 10._wp
    end do
  end subroutine normalize_value
!------------------------------------------------------------------------------
  !--------------
  ! transfer data
  !--------------
  SUBROUTINE set_data_1  (grib, DATA)
  TYPE (t_grib1) ,INTENT(inout) :: grib
  REAL(wp)       ,INTENT(in)    :: DATA (:)

    INTEGER :: n
    n = SIZE(DATA)

    IF (n /= grib% isec4% n_data) THEN
      WRITE(0,*)'set_data_1, invalid data size:, n, isec4 =',&
                n,grib% isec4% n_data
      CALL finish ('set_data_1','invalid data size')
    ENDIF
    CALL reallocate_data (grib, n)
    grib% rsec4 (1:n) = DATA
  END SUBROUTINE set_data_1
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE set_data_2  (grib, DATA)
  TYPE (t_grib1) ,INTENT(inout) :: grib
  REAL(wp)       ,INTENT(in)    :: DATA (:,:)

    INTEGER :: n
    n = SIZE(DATA)

    IF (n /= grib% isec4% n_data) THEN
      WRITE(0,*)'set_data_2, invalid data size:, nx,ny,n,isec4 =',&
        shape(data), n, grib% isec4% n_data
      CALL finish ('set_data_2','invalid DATA size')
    ENDIF
    CALL reallocate_data (grib, n)
    grib% rsec4 (1:n) = RESHAPE (DATA,(/n/))
  END SUBROUTINE set_data_2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE set_data_3  (grib, DATA)
  TYPE (t_grib1) ,INTENT(inout) :: grib
  REAL(wp)       ,INTENT(in)    :: DATA (:,:,:)

    INTEGER :: n
    n = SIZE(DATA)

    IF (n /= grib% isec4% n_data) THEN
      WRITE(0,*)'set_data_3, invalid data size:, n1,n2,n3,n,isec4 =',&
        shape(data), n, grib% isec4% n_data
      CALL finish ('set_data','invalid data size')
    ENDIF
    CALL reallocate_data (grib, n)
#ifndef __SX__
    grib% rsec4 (1:n) = RESHAPE (DATA,(/n/))
#else
    if (size(DATA,2) == 1 .and. size(DATA,3) == 1) then
       grib% rsec4 (1:n) = DATA(1:n, lbound(data,2), lbound(data,3))
    else
       grib% rsec4 (1:n) = RESHAPE (DATA,(/n/))
    end if
#endif
  END SUBROUTINE set_data_3
!==============================================================================
  SUBROUTINE open_gribfile (grib, file, mode, edition, kret)
  !---------------------------------------
  ! open a gribfile
  !
  !   for writing: set some default values
  !---------------------------------------
  TYPE (t_grib1)    ,INTENT(inout)          :: grib
  CHARACTER (len=*) ,INTENT (in)            :: file
  CHARACTER (len=*) ,INTENT (in)            :: mode     ! 'r','r+','w','a'
  integer           ,INTENT (in)  ,OPTIONAL :: edition  ! GRIB edition
  INTEGER           ,INTENT (out) ,OPTIONAL :: kret

    integer :: ed       ! GRIB edition to use
    logical :: lw       ! write mode?

    ed = grib_edition
    if (present (edition)) ed = edition
    if (ed < 1 .or. ed > 2) &
         call finish ("open_gribfile","supported editions are 1 and 2")
    lw = .false.
    !---------------------------
    ! call the low level routine
    !---------------------------
    select case (mode)
    case ('w','a')
      if (ed > 1 .and. grib_library == 1) &
        call finish ("open_gribfile","CGRIBEX can only write edition 1")
      if (grib_library == 1 .or. d_pbio) then
        CALL pbopen (grib, file, mode, kret)
        if (ed == 1) then
          CALL pbsetraw (grib, d_r)
        else
          CALL pbsetraw (grib, 1)       ! Always raw mode for edition=2
        end if
      else
#ifdef GRIB_API
        call grib_open_file (grib% file, file, mode)
        grib% pbio = .false.            ! Don't use PBIO for writing
#else
        call finish ("open_gribfile","GRIB_API not linked!")
#endif
      end if
      lw = .true.
    case default
      CALL pbopen (grib, file, mode, kret)
    end select
    !---------------
    ! error handling
    !---------------
    IF (PRESENT(kret)) THEN
      IF (kret /= 0) RETURN
    ENDIF
    !-------------
    ! set defaults
    !-------------
    grib% isec0% n_octets    =   0
    grib% isec0% edition     =  ed
    !--------------------------------------
    ! Section 1: Product Definition Section
    !--------------------------------------
    grib% isec1% center      = d_c ! Identification of centre.
    grib% isec1% sub_center  = d_s ! Sub-centre identifier.
    grib% isec1% process     = d_p ! Generating process identification.
    grib% isec1% code        = 255 ! Parameter indicator.
    grib% isec1% table       = 255 ! Version number of code table.
    grib% isec1% century     = 255 ! Reference time of data.
    grib% isec1% year        = 255 !                  "
    grib% isec1% day         = 255 !                  "
    grib% isec1% month       = 255 !                  "
    grib% isec1% hour        = 255 !                  "
    grib% isec1% minute      = 255 !                  "
    IF (d_rt% days /= 0) CALL set_ref_time (grib, d_rt) ! "
    grib% isec1% time_range  =   0 ! Time range indicator.
    grib% isec1% time_unit   = 255 ! Time unit indicator.
    grib% isec1% p1          =   0 ! Time period.
    grib% isec1% p2          =   0 ! Time period.
    grib% isec1% level_type  = 255 ! Type of level indicator.
    grib% isec1% level_st    =   0 ! Single level or top of layer.
    grib% isec1% level_b     =   0 ! Bottom of layer.
    !--------------------------------------------
    !    these values are not changed up to now :
    !--------------------------------------------
    grib% isec1% grid        = 255 ! Grid description follows in section 2.
    grib% isec1% present_2_3 = 128 ! 128: Section 2; 64: Section 3 included.
    grib% isec1% n_average   =   0 ! Number of products included in an average.
    grib% isec1% n_missing   =   0 ! Number of products missing from an average
    grib% isec1% factor      =   0 ! Decimal scale factor.
    grib% isec1% local_flag  =   0 ! 0: No local use, 1: Local use in Section 1
    grib% isec1% reserved    =   0 !
    grib% isec1% local_ident =   0 ! local GRIB use definition identifier.
    grib% isec1% local_use   =   0 ! used for ECMWF local extensions.
    !----------------------------
    ! Section 1, local part (DWD)
    !----------------------------
    grib% s1_dwd% year       = 255 ! Creation date (default)
    grib% s1_dwd% day        = 255 !       "
    grib% s1_dwd% month      = 255 !       "
    grib% s1_dwd% hour       = 255 !       "
    grib% s1_dwd% minute     = 255 !       "

!print *, "open_gribfile", d_c, d_l, d_n_en, d_i_en, d_e, d_expid
    select case (d_l)
    case (254)
       call set_dwd_local (grib, d_t, d_run, d_e)
    case (152,153,252,253)
       call set_dwd_local (grib, d_t, d_run, d_e)
       call set_dwd_ens   (grib, d_n_en, d_i_en, d_id_en)
    case (1)
       if (d_c == WMO0_ECMWF) then
          if (d_e >= 0) write (d_expid,'(i4.4)') mod (d_e, 10000)
          call set_ecmwf_ens (grib, d_n_en, d_i_en, d_expid)
       end if
    end select
    if (grib% s1_dwd% year < 0 .or. grib% s1_dwd% year > 255) &
      grib%s1_dwd%year       = 255 ! Mark out-of-bounds values
    !-----------------------------------
    ! Section 2: Grid Definition Section
    !-----------------------------------
    grib% isec2              = 255
    grib% latlon% nvcp       =   0
    grib% gauss%  nvcp       =   0
    grib% tri%    nvcp       =   0
    grib% sph%    nvcp       =   0

    grib% rsec2%  rot_angle  =   0._wp
    grib% rsec2%  str_factor =   1._wp
    grib% rsec2%  reserved   =   0._wp
    grib% rsec2%  vcp        =   0._wp
    !-------------------------------------
    ! Section 3: Bitmap Definition Section
    !-------------------------------------
    grib% isec3%  bitmap     =   0
    grib% isec3%  missing    =   0

    grib% rsec3%  ignored    =   0._wp
    grib% rsec3%  missing    =   0._wp
    !------------------------
    ! Section 4: Data Section
    !------------------------
    grib% isec4%  n_data     =   0 ! Number of data values in array PSEC4.
    grib% isec4%  bits       = d_b ! Number of bits used for encoded values.
    grib% isec4%  grid_type  =   0 ! 0: Grid point data, 128: Sph.harm.coeff.
    grib% isec4%  packing    =   0 ! 0: Simple packing,   64: Complex packing.
    grib% isec4%  data_repr  =   0 ! 0: Floating point,   32: Integer data.
    grib% isec4%  flags      =   0 ! 0: no                16: Additional flags.
    grib% isec4%  reserved   =   0 ! Reserved. Set to 0.
    grib% isec4%  matrix     =   0 ! 0: single datum      64: Matrix
    grib% isec4%  bitmap2    =   0 ! 0: no Secondary bit. 32: S.bitmaps present
    grib% isec4%  width      =   0 ! 0: constant width    16: different width
    grib% isec4%  bits2      =   0 ! Number of bits for second order values.
    grib% isec4%  wmo        =   0 ! Reserved for WMO reserved flag fields.
    grib% isec4%  start_cplx =   0 ! For complex packing, start of packed data.
    grib% isec4%  scale_cplx =   0 ! For complex packing, the scaling factor P.
    grib% isec4%  j_cplx     =   0 ! pentagonal resolution parameter J.
    grib% isec4%  k_cplx     =   0 ! pentagonal resolution parameter K.
    grib% isec4%  m_cplx     =   0 ! pentagonal resolution parameter M.
    grib% isec4%  non_miss   =   0 ! The number of non-missing values.
    grib% isec4%  reserved2  =   0 ! Reserved. Set to 0.
    grib% isec4%  o_coded    =   0 ! offset (returned by 'G', 'I' or 'J')
    grib% isec4%  remaining  =   0 !

#ifdef GRIB_API
    !----------------------------------------------------------------
    ! Set up basic template for GRIB_API, depending on chosen edition
    !----------------------------------------------------------------
    if (lw .and. grib_library == 2) then
       select case (ed)
       case (1)
          call setup_grib1_template (grib)
       case (2)
          call setup_grib2_template (grib)
       end select
    end if
#endif
  END SUBROUTINE open_gribfile
!-----------------------------------------------------------------------------
#ifdef GRIB_API
  subroutine setup_grib1_template (grib)
    type(t_grib1), intent(inout) :: grib
    !--------------------------------
    ! Set up basic template for GRIB1
    !--------------------------------
    integer :: h        ! handle

    call grib_new_from_samples (h, "GRIB1")
    grib% handle = h

    !--------------------------------------
    ! Section 1: Product Definition Section
    !--------------------------------------
    call grib_set (h,'centre'                     ,grib% isec1% center)
    call grib_set (h,'subCentre'                  ,grib% isec1% sub_center)
    call grib_set (h,'generatingProcessIdentifier',grib% isec1% process)
    call grib_set (h,'table2Version'              ,grib% isec1% table)
    call grib_set (h,'indicatorOfParameter'       ,grib% isec1% code)
    call grib_set (h,'indicatorOfTypeOfLevel'     ,grib% isec1% level_type)
    call grib_set (h,'level'                      ,0)
!   call grib_set (h,'dataDate',date1)
!   call grib_set (h,'dataTime',time1)
    call grib_set (h,'stepType'                   ,'instant')
    call grib_set (h,'timeRangeIndicator'         ,d_ran    )

    if (grib% isec1% local_flag > 0) then
       call grib_set (h,'setLocalDefinition'   ,1)
       call grib_set (h,'localDefinitionNumber',grib% isec1% local_ident)
    end if
    select case (grib% isec1% center)
    case (WMO0_DWD,WMO0_COSMO,WMO0_MSWISS,WMO0_COMET)
       if (grib% isec1% local_flag > 0) then
          call grib_set (h,'localElementNumber'   ,255 ) !grib% s1_dwd% element_no
          call grib_set (h,'localVersionNumber'   ,d_e + 2**14 * grib% s1_dwd% run_type)

          call grib_set (h,'localDecodeDateYear'  ,grib% s1_dwd% year)
          call grib_set (h,'localDecodeDateMonth' ,grib% s1_dwd% month)
          call grib_set (h,'localDecodeDateDay'   ,grib% s1_dwd% day)
          call grib_set (h,'localDecodeDateHour'  ,grib% s1_dwd% hour)
          call grib_set (h,'localDecodeDateMinute',grib% s1_dwd% minute)
       end if
       select case (grib% isec1% local_flag)
       case (252,253)
          call grib_set (h,'localEnsembleIdentification'      ,grib% s1_dwd% ensemble_id)
          call grib_set (h,'localNumberOfEnsembleMembers'     ,grib% s1_dwd% ensemble_size)
          call grib_set (h,'localActualNumberOfEnsembleNumber',grib% s1_dwd% ensemble_no)
       end select
    case (WMO0_ECMWF)
    end select

    ! Section 2
    call grib_set (h,'numberOfVerticalCoordinateValues',0)
    ! Grid definition follows: default is regular_ll, S->N
    call grib_set (h,'dataRepresentationType',0)
    call grib_set (h,'latitudeOfFirstGridPointInDegrees',-90.0)
    call grib_set (h,'latitudeOfLastGridPointInDegrees' , 90.0)
    call grib_set (h,'jScansPositively', 1)

    ! Section 4
    call grib_set (h,'packingType' ,'grid_simple')
    call grib_set (h,'bitsPerValue', grib% isec4% bits)

    ! Force full coding of constant fields (no packing)
    call grib_gribex_mode_on ()

  end subroutine setup_grib1_template
!-----------------------------------------------------------------------------
  subroutine setup_grib2_template (grib)
    type(t_grib1), intent(inout) :: grib
    !--------------------------------
    ! Set up basic template for GRIB2
    !--------------------------------
    integer            :: h              ! handle
    integer            :: vers           ! local version number
    integer            :: tabvers        ! tablesVersion
    integer            :: mintabvers = 5 ! min. tablesVersion for GME
!   integer            :: mintabvers = 6 ! min. tablesVersion for ICON
    integer            :: togp           ! typeOfGeneratingProcess
    integer            :: gpid           ! generatingProcessIdentifier

    call grib_new_from_samples (h, "GRIB2")
    grib% handle = h

    !----------
    ! Section 1
    !----------
    ! grib1/0.table
    call grib_set (h,'centre'                     ,grib% isec1% center)
    call grib_set (h,'subCentre'                  ,grib% isec1% sub_center)

    call grib_get (h, 'tablesVersion'             ,tabvers)
    mintabvers = max (gribtabversion, mintabvers)
    if (tabvers < mintabvers) then
       call grib_set (h, 'tablesVersion'          ,mintabvers)
    end if

    ! grib2/tables/4/1.2.table :: 0=ana,1=fc_start,2=veri_time
    call grib_set (h,'significanceOfReferenceTime', 0)
!   call grib_set (h,'dataDate',date1)
!   call grib_set (h,'dataTime',time1)

    ! grib2/tables/5/1.3.table   :: 0=oper,1=test,2=exp
    select case (d_e)
    case (0:50)
       call grib_set (h,'productionStatusOfProcessedData', 0)
    case (51:99)
       call grib_set (h,'productionStatusOfProcessedData', 1)
    case default
       call grib_set (h,'productionStatusOfProcessedData', 2)
    end select
    ! grib2/1.4.table   :: 0=ana,1=fc,4:5=pert.fc,...
    call grib_set (h,'typeOfProcessedData', 0)
    call grib_set (h,'stepType','instant')

    !---------------------
    ! Section 2: local use
    !---------------------
    if (grib% isec1% local_flag > 0) then
       call grib_set (h,'grib2LocalSectionPresent', 1)
       call grib_set (h,'localDefinitionNumber',grib% isec1% local_ident)
       select case (grib% isec1% center)
       case (WMO0_DWD,WMO0_COSMO,WMO0_MSWISS,WMO0_COMET)
!         call grib_set_missing (h,'localHostIdentifier')
          call grib_set (h,'localHostIdentifier'    ,0) ! 1=oper.,2=backup,???
          call grib_set (h,'localCreationDateYear'  ,grib% s1_dwd% year + 1900)
          call grib_set (h,'localCreationDateMonth' ,grib% s1_dwd% month)
          call grib_set (h,'localCreationDateDay'   ,grib% s1_dwd% day)
          call grib_set (h,'localCreationDateHour'  ,grib% s1_dwd% hour)
          call grib_set (h,'localCreationDateMinute',grib% s1_dwd% minute)
          call grib_set (h,'localCreationDateSecond',0)

          call grib_set (h,'localNumberOfExperiment',d_e)
          call grib_set (h,'localInformationNumber' ,255)  ! ZEN
          select case (grib% isec1% local_ident)
          case (254)
!            vers = grib_get_api_version ()
             vers = 10000 * major() + 100 * minor() ! currently no patchlevel
             call grib_set (h,'localVersionNumber'  ,vers)
          case (152,153,252,253)
             call grib_set (h,'localTypeOfEnsembleForecast',0)
          end select
       end select
    end if
    !------------------------------
    ! Section 3 (GDS)
    ! template.3.*: grid definition
    !------------------------------
    call grib_set (h,'numberOfDataPoints',0)   ! will be set later
    ! grib2/3.11.table
    call grib_set (h,'interpretationOfNumberOfPoints',0)

    ! grib2/3.1.table
    call grib_set (h,'gridType','regular_ll')  ! Default
    ! grib2/3.2.table :: r=6371.229 km
    call grib_set (h,'shapeOfTheEarth',6)
    call grib_set (h,'latitudeOfFirstGridPointInDegrees',-90.0)
    call grib_set (h,'latitudeOfLastGridPointInDegrees' , 90.0)
    call grib_set (h,'Ni'      ,0)
    call grib_set (h,'Nj'      ,0)
    call grib_set (h,'jScansPositively', 1)

    !----------
    ! Section 4 (mostly)
    !----------
    ! template.4
    ! grib2/4.0.table :: 0=Deterministic, 1=Ensemble, 2=Mean/Spread,7=fg/anaerr...
    grib% pdtn = GRIB_4_0_ANA_FC
    call grib_set (h,'productDefinitionTemplateNumber', grib% pdtn)

    ! grib2/0.0.table
    call grib_set (h,'discipline'       ,0)
    ! grib2/4.1.*.table
    call grib_set (h,'parameterCategory',0)
    ! grib2/4.2.*.*.table
    call grib_set (h,'parameterNumber'  ,0)

    !------------------------------------------------------------------
    ! grib2/4.3.table :: 0=ana,1=ini,2=fc,4=ens.fc,7=anaerr,9=clim.,...
    ! Derive typeOfGeneratingProcess and generatingProcessIdentifier
    !------------------------------------------------------------------
    gpid = grib% isec1% process
    togp = grib% isec1% process
    if (gpid == d_p .and. d_p == d_ptyp) then
       ! Derive from known DWD-defined generating process identifiers:
       ! (see grib2/tables/local/edzw/1/generatingProcessIdentifier.table)
       select case (grib% isec1% process)
       case (26,27,189:196,210)
          togp = GRIB_4_3_ANA_ERR
       case (1:24,100,110,120,130,137,139,141,143,145,147,149,  &
             161,163,165,167,169,171,173,175,177,179,181,206,255)
          togp = GRIB_4_3_ANA
       case (GRIB_4_3_ANA_INC)  ! (201) Difference analysis - first guess
          togp = GRIB_4_3_ANA_INC
       case default
          togp = GRIB_4_3_FC
       end select
    else
       togp = d_ptyp
    end if
    call grib_set (h,'typeOfGeneratingProcess', togp)

    ! grib2/backgroundProcess.table :: 0=main,1=pre,2=ana,3=test
    if (grib% s1_dwd% local_ident > 0) then
       call grib_set (h,'backgroundProcess', grib% s1_dwd% run_type)
    else
       call grib_set (h,'backgroundProcess', 255)
    end if
    ! call grib_set (h,'backgroundGeneratingProcessIdentifier', 0)
    ! grib2/generatingProcessIdentifier.table
    call grib_set (h,'generatingProcessIdentifier', gpid)

    call grib_set_missing (h,'hoursAfterDataCutoff')
    call grib_set_missing (h,'minutesAfterDataCutoff')

    ! grib2/4.4.table
    call grib_set (h,'indicatorOfUnitOfTimeRange',1)
    call grib_set (h,'stepUnits',1)

    ! Requires productDefinitionTemplateNumber=1,11,...
    !call grib_set (h,'typeOfEnsembleForecast', 192)
    !call grib_set (h,'numberOfForecastsInEnsemble', 32)
    !call grib_set (h,'perturbationNumber', 17)

    call grib_set (h,'numberOfValues', 0)
    !----------
    ! Section 5
    !----------
    call grib_set (h,'packingType' ,'grid_simple')
    !----------
    ! Section 7
    !----------
    call grib_set (h,'bitsPerValue', grib% isec4% bits)

  end subroutine setup_grib2_template
#endif
!-----------------------------------------------------------------------------
  SUBROUTINE write_gribrecord (grib, edition)
  TYPE (t_grib1) ,INTENT(inout) :: grib     ! GRIB record to write
  INTEGER ,OPTIONAL             :: edition  ! GRIB edition to use

    INTEGER                :: ed
#ifdef GRIB_API
    TYPE (t_grib1)         ::  grib2    ! GRIB 2 record
    INTEGER                :: igrib2    ! GRIB 2 API handle
    INTEGER                :: int_size
    CHARACTER ,ALLOCATABLE :: message (:)
    INTEGER(kindOfSize)    :: byte_size ! size of message in byte
#endif
    !----------------------------
    ! check GRIB edition to write
    !----------------------------
    ed = grib_edition
    if (present (edition)) ed = edition
    select case (ed)
    case (1,2)
    case default
      ed = d_edition
    end select
#ifdef GRIB_API
    select case (grib_library)
    case (1)                            ! library=1
    !---------------------------------
    ! Use CGRIBEX for primary encoding
    !---------------------------------
#endif
    !-----------------------
    ! encode the GRIB record
    !-----------------------
    CALL gribex (grib, 'C')
    !---------------------
    ! finally write record
    !---------------------
#ifdef GRIB_API
    select case (ed)
    case (1)                                 ! library=1 edition=1
#endif
      CALL pbwrite (grib)
#ifdef GRIB_API
    case (2)                                 ! library=1 edition=2
      !-------------------------
      ! convert GRIB 1 -> GRIB 2
      !-------------------------
      call grib1togrib2 (grib, igrib2)
      !---------------------------------
      ! convert from GRIB API to CGRIBEX
      !---------------------------------
      call grib_get_message_size (igrib2, byte_size)
      int_size = (byte_size + 3) / 4
      allocate (message (int_size*4))
      call reallocate_buffer     ( grib2, int_size )
      call grib_copy_message     (igrib2, message  )
      call grib_release          (igrib2)
      grib2% kgrib = transfer (message, grib2% kgrib)
      deallocate (message)
      !---------------------
      ! finally write record
      !---------------------
      grib2% blen = byte_size
      grib2% file = grib% file
      CALL pbwrite  (grib2)
      call destruct (grib2)
    end select
#endif
#ifdef GRIB_API
    case (2)                                 ! library=2
      if (grib% handle == INVALID_HANDLE) &
           call finish ("write_gribrecord","invalid GRIB_API handle")
      !----------------------------------------
      ! Clone template before encoding the data
      !----------------------------------------
FTRACE_BEGIN("write_gribrecord:grib_clone")
      call grib_clone (grib% handle, igrib2)
FTRACE_END  ("write_gribrecord:grib_clone")
FTRACE_BEGIN("write_gribrecord:set_values")
      call grib_set   (igrib2,'bitsPerValue',grib% isec4% bits)
      call grib_set   (igrib2,'values',      grib% rsec4(1:grib% isec4% n_data))

      select case (grib_packing)
      case ("")                              ! No packing
      case default
        call grib_set (igrib2,'packingType' ,trim(grib_packing))
      case ("grid_second_order")
        !------------------------------------------------
        ! If grid_second_order packing requested, enforce
        ! reencoding of values (GRIB_API bug or feature?)
        !------------------------------------------------
        call grib_set (igrib2,'packingType' ,trim(grib_packing))
        call grib_set (igrib2,'values',      grib% rsec4(1:grib% isec4% n_data))
      end select
FTRACE_END  ("write_gribrecord:set_values")
      if (grib% pbio) then                   ! library=2 pbio=T
      !-------------------------------------
      ! Write to file using the PB* routines
      !-------------------------------------
        call grib_get_message_size (igrib2, byte_size)
        int_size = (byte_size + 3) / 4
        allocate (message (int_size*4))
        call reallocate_buffer     ( grib2, int_size )
FTRACE_BEGIN("write_gribrecord:copy_message")
        call grib_copy_message     (igrib2, message  )
FTRACE_END  ("write_gribrecord:copy_message")
FTRACE_BEGIN("write_gribrecord:grib_release")
        call grib_release          (igrib2)
FTRACE_END  ("write_gribrecord:grib_release")
FTRACE_BEGIN("write_gribrecord:transfer")
        grib2% kgrib = transfer (message, grib2% kgrib)
FTRACE_END  ("write_gribrecord:transfer")
        deallocate (message)
        !---------------------
        ! finally write record
        !---------------------
        grib2% blen = byte_size
        grib2% file = grib% file
FTRACE_BEGIN("write_gribrecord:pbwrite")
        CALL pbwrite  (grib2)
        call destruct (grib2)
FTRACE_END  ("write_gribrecord:pbwrite")
      else                                  ! library=2 pbio=F
FTRACE_BEGIN("write_gribrecord:grib_write")
         call grib_write   (igrib2, grib% file)
         call grib_release (igrib2)
FTRACE_END  ("write_gribrecord:grib_write")
      end if
    end select
#endif
  END SUBROUTINE write_gribrecord
!-----------------------------------------------------------------------------
  SUBROUTINE close_gribfile (grib)
  TYPE (t_grib1) ,INTENT(inout) :: grib
#ifdef GRIB_API
    if (grib% handle /= INVALID_HANDLE) call grib_release (grib% handle)
    grib% handle = INVALID_HANDLE
    if (grib% pbio) then
       CALL pbclose (grib)
    else
       call grib_close_file (grib% file)
       call destruct        (grib)
    end if
#else
    CALL pbclose (grib)
#endif
  END SUBROUTINE close_gribfile
!=============================================================================
  subroutine grib_api_version (major, minor, revision)
    ! GRIB API library version number (or -1 if not available)
    integer, intent(out)           :: major
    integer, intent(out)           :: minor
    integer, intent(out), optional :: revision
#ifdef GRIB_API
    integer :: triplet
    triplet = grib_get_api_version ()
    major                            =      triplet / 10000
    minor                            = mod (triplet / 100  , 100)
    if (present (revision)) revision = mod (triplet        , 100)
#else
    major = -1
    minor = -1
    if (present (revision)) revision = -1
#endif
  end subroutine grib_api_version
!=============================================================================
END MODULE mo_grib_handling
