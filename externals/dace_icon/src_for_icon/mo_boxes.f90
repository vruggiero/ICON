!
!+ derive distribution of observations over 'boxes' and processor elements
!
MODULE mo_boxes
!
! Description:
!   Derive the distribution of observations over boxes (used for
!   first-guess check and preconditioning) and processor elements.
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
! V1_3         2008/12/08 Andreas Rhodin
!  Work around NEC SX compiler bug (write format b31)
! V1_4         2009/03/26 Andreas Rhodin
!  fix for zero number of observations in a box
! V1_5         2009/05/25 Harald Anlauf
!  release_bx: optimize for SX-9
! V1_7         2009/08/24 Harald Anlauf
!  set_boxes: fix crash in box reordering if some boxes are empty
! V1_9         2010/04/20 Harald Anlauf
!  set_boxes: optimize distribution of obs->box and box->proc for GPSRO.
!             bugfix (avoid observation loss).
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  reorder boxes even if nproc_repro is set
! V1_14        2011/11/08 Harald Anlauf
!  Workaround for sxf90 bug with PACK(character array)
! V1_19        2012-04-16 Harald Anlauf
!  bugfix for possibly uninitialized nproc_use
! V1_20        2012-06-18 Harald Anlauf
!  Framework for arbitrary number of input feedback files (up to 999)
! V1_22        2013-02-13 Andreas Rhodin
!  changes for STD operator; treat OBS_ONLY reports;
!  fix double entries for DISSMISSed reports
! V1_29        2014/04/02 Andreas Rhodin
!  consistently use n_ot instead of n_obstype
! V1_31        2014-08-21 Andreas Rhodin
!  diagnose empty feedback-(fof-)files and skip reading them
! V1_42        2015-06-08 Andreas Rhodin
!  preparations for temporal interpolation for COSMO MEC
! V1_43        2015-08-19 Andreas Rhodin
!  set_cosmo_boxes: abort in case of error
! V1_44        2015-09-30 Harald Anlauf
!  move type t_bufr_inv to mo_t_obs to simplify dependencies
! V1_47        2016-06-06 Harald Anlauf
!  Set tentative 'cost' for GPSGB data
! V1_51        2017-02-24 Andreas Rhodin
!  accelerate subroutine set_boxes
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD 2005-2008  original source
! Harald Anlauf   DWD 2008       optimizations and bug fixes (IBM+NEC)
! Gerhard Paul    DWD 2008       changes for NetCDF input of observations
!==============================================================================

  !-------------
  ! Modules used
  !-------------
  use mo_kind,         only: wp, i8         ! kind parameters
  use mo_exception,    only: finish         ! error exit routine
  use slatec_module,   only: sort           ! sort
  use mo_obs_tables,   only: obstyp,       &! observation types, names (table)
                             write_report, &! write report to status file
                             write_pending,&! write pending reports
                             pass_stat,    &! passive  reports statistics table
                             reje_stat,    &! rejected reports statistics table
                             oonly_stat     ! observation only statistics table
  use mo_obs_bufr_dwd, only: read_obs_bufr  ! read BUFR files
  use mo_obs_netcdf,   only: scan_obs_netcdf! scan NetCDF observation files
  use mo_fdbk_tables,  only: n_ot,         &! number of observation types
                             OT_SYNOP,     &! observation type identifier
                             OT_AIREP,     &!   ..
                             OT_SATOB,     &!   ..
                             OT_DRIBU,     &!   ..
                             OT_TEMP,      &!   ..
                             OT_PILOT,     &!   ..
                             OT_SATEM,     &!   ..
                             OT_PAOB,      &!   ..
                             OT_SCATT,     &!   ..
                             OT_GPSRO,     &!   ..
                             OT_GPSGB,     &!   ..
                             OT_RAD         !   ..
  use mo_fdbk_in,      only: read_feedback  ! scan feedback file
  use mo_t_obs,        only: t_obs,        &! observation data type
                             t_icol,       &! column data type
                             t_spot,       &! report meta data type
                             t_box,        &! component of t_spot
                             set_xuv,      &! set components of t_icol datatype
                             set_inbx,     &! set spot% inbx
!                            m_source,     &! max no. of source files
                             n_source,     &!     no. of source files to read
                             source,       &! table of observation files
                             bufr_inv,     &! observation file inventory data
                             p_bcast,      &! generic MPI bcast routine
                             read_bufr,    &! flag to read BUFRs
                             read_netcdf,  &! flag to read NetCDF files
                             ft_name,      &! derive name of file type
                             netcdf_verb    ! verbosity of NetCDF decoding
  use mo_mpi_dace,     only: dace,         &! MPI group info
                             p_bcast,      &! generic MPI broadcast routine
                             p_sum,        &! generic MPI reduction routine
                             p_gather,     &! generic MPI gather routine
!                            p_min,        &! generic MPI minimum over PEs
                             p_max          ! generic MPI maximum over PEs
  use mo_run_params,   only: p_readbufr,   &! PE used to read BUFR data
                             npe_read_obs, &! number of PEs to read obsv
                             nproc_repro,  &! fict. PEs for repeatable runs
                             method         !
  use mo_p_output,     only: oline,        &! output line buffer
                             iol,          &! number of next line to write
                             nextline,     &! routine to increment line number
                             flush_buf      ! routine to write buffer
  use mo_t_use,        only: STAT_PASSIVE, &!
                             STAT_PAS_REJ, &!
                             STAT_OBS_ONLY,&!
                             STAT_REJECTED,&!
                             STAT_DISMISS   !
  use mo_physics,      only: rearth         ! earth radius
  use mo_cosmo_obs,    only: scan_cosmo_obs ! scan COSMO observation files
  use mo_atm_state,    only: t_atm          ! atmospheric state derived type
! use mo_cpu_time,     only: stop_time      ! breakdown of times in set_boxes
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: set_input_boxes ! associate observations with boxes for input
  public :: set_boxes       ! associate observations with boxes for PSAS
  public :: set_veri_boxes  ! associate observations with boxes for MEC


  !-----------------------------------------------------------------------
  ! Derived type definition:
  ! This structure is used to hold information concerning a report, relevant
  ! for the distribution of reports over boxes and processor elements
  !-----------------------------------------------------------------------
  type t_bx
    real(wp)    :: d1    ! distance from arbitrary initial point
    real(wp)    :: d2    ! distance
    type(t_box) :: box   ! box number (in FG scan)
    real(wp)    :: x(3)  ! coordinates (unit vector on the sphere)
    integer     :: count ! number of independent observations
    integer     :: otyp  ! observation type
    integer     :: cost  ! cost of observation processing
    integer     :: isrc  ! index of element in local array of type't_bx'
    logical     :: valid ! true if not yet distributed
  end type t_bx

  !----------------------------------------------------------
  ! Cost model for observation types depending on method/step
  !----------------------------------------------------------
  enum, bind(c)
     enumerator :: COST_OBS_READ = 0,  &! File read
                   COST_OBS_FG   = 1,  &! First-guess step
                   COST_OBS_ASS  = 2,  &! Variational assimilation
                   COST_OBS_VERI = 3    ! Verification step
  end enum
!------------------------------------------------------------------------------
contains
!==============================================================================
  subroutine set_input_boxes
  !----------------------------------------------------------------
  ! Decide which observational data is read on which PE. Data is
  ! redistributed later for the First Guess check and again for the
  ! analysis step.
  !
  ! Observation types are either read on one PE only (npe_read_obs < 1) or
  ! distributed round robin over processor elements (npe_read_obs >
  ! p_nprocs) or read on dedicated PEs.
  !----------------------------------------------------------------
    !================
    ! local variables
    !================
    integer          :: pe        ! processor index
    integer          :: ir        ! observation type index
    integer          :: if        ! index for obsv.types to handle first
    integer          :: is        ! file index
    integer          :: codetype  ! CMA codetype
    logical          :: first     ! true  for obsv.types to handle first
    integer          :: k_source  ! max. no. source files (for pretty printing)
    integer          :: subsets
    integer          :: subsetso (n_ot)
    integer          :: subsetl
    integer          :: subsetn
    integer          :: pe1, pel, pen
    integer          :: ier
    character(len=8) :: oname
    logical          :: l_dedic  (n_ot)             ! use dedicated pe?
    !----------------------------------------------------------------
    ! estimate cost of observation processing (reading and decoding):
    !----------------------------------------------------------------
    real(wp)         :: cost_obs     (n_ot) = 1._wp ! cost / 1 report
    real(wp)         :: cost_obs_tot (n_ot)         ! cost / all reports
    real(wp)         :: cost_pe      (0:dace% npe-1)! cost / processor
    real(wp)         :: cost_mean_pe                ! mean cost / PE
    real(wp)         :: cost_total
    real(wp)         :: cost_left
    real(wp)         :: cost_file    (n_source)
    integer          :: i_cost       (n_source)     ! sorted file ids.
    !----------------------------------------------
    ! list of observation types to process first
    ! and to read from dedicated processor elements
    !----------------------------------------------
    integer          :: use_first (7) = (/OT_TEMP,OT_PILOT,OT_DRIBU,OT_SYNOP,OT_GPSGB,OT_GPSRO,OT_RAD/)
    integer          :: use_dedic (7) = (/      0,       0,       1,       2,       0,      -1,    -1/)

    type t_dedic
       integer :: obstype     = -1   ! observation type (3D-Var module)
       integer :: codetype(2) = -1   ! CMA code type (range)
       integer :: pe          = -1   ! dedicated processor index (-1: none)
    end type t_dedic

    type(t_dedic), parameter :: obs_pe(14) = [ &
         t_dedic( OT_TEMP,        -1,   0 ), &! TEMP catchall
         t_dedic( OT_PILOT, [132,137],  5 ), &! PILOT wind profiler
         t_dedic( OT_PILOT,      157,   5 ), &! SCADA wind profiler
         t_dedic( OT_PILOT,       -1,   0 ), &! PILOT catchall
         t_dedic( OT_DRIBU,       -1,   1 ), &! SYNOP sea: BUOY
         t_dedic( OT_SYNOP, [ 21, 24],  1 ), &! SYNOP sea: SHIP (new)
!        t_dedic( OT_SYNOP, [ 21, 24],  2 ), &! SYNOP sea: SHIP (old)
         t_dedic( OT_SYNOP, [ 11, 14],  2 ), &! SYNOP land
         t_dedic( OT_SYNOP,      140,   2 ), &! SYNOP land: METAR
         t_dedic( OT_SYNOP,       20,   4 ), &! SYNOP land: CMAN
         t_dedic( OT_SYNOP,       16,   6 ), &! SYNOP land: SNOW
         t_dedic( OT_SYNOP,       -1,   2 ), &! SYNOP catchall
         t_dedic( OT_GPSGB,       -1,   3 ), &! GNSS ground-based (new)
!        t_dedic( OT_GPSGB,       -1,   0 ), &! GNSS ground-based (old)
         t_dedic( OT_GPSRO,       -1,  -1 ), &! GNSS RO
         t_dedic( OT_RAD,         -1,  -1 )  ]! RADiances
    !---------------------------------------------------------------
    ! scan BUFR and NetCDF files, store file inventory in 'bufr_inv'
    !---------------------------------------------------------------
    if (read_NetCDF) call scan_obs_netcdf
    if (read_bufr )  call read_obs_bufr
    call scan_cosmo_obs
    call read_feedback (pass=1)
    call p_bcast (bufr_inv,         dace% pio)
    call p_bcast (source% entries,  dace% pio)
    call p_bcast (source% filetype, dace% pio)

    !-----------------------------
    ! print file numbers and names
    !-----------------------------
    call flush_buf
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a,i0)') '  Number of observation input files: ',n_source
      write(6,'()')
      write (6,'(2x,2x,1x,a44,4a10)') &
        'File','filetype','obstype','entries','complete'
      do if = 1, n_source
        write (6,'(1x,i3,1x,a44,a10,2i10,9x,l1)') if, &
                      trim(source(if)% file),     &
                   ft_name(source(if)% filetype), &
                           source(if)% obstype,   &
                           source(if)% entries,   &
                           source(if)% complete
      end do
    endif
    !----------------
    ! print inventory
    !----------------
    call flush_buf
    if (dace% lpio) then
      k_source = max (((n_source+9)/10)*10, 50) ! Dynamic list format
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a,99(i5,:,5x))') ' obstype    records  subsets file:', &
           (       if    , if=1, k_source/10)
      write(6,'(a,999i1)'      ) '                             '     , &
           (modulo(if,10), if=1, k_source   )
      write(6,'()')
      do ir = 1, n_ot
        if (bufr_inv(ir)% nrec == 0 .and. bufr_inv(ir)% nsubset == 0) cycle
        write(6,'(2x,a,2i9,1x,999a1)') obstyp(ir)%name, &
          bufr_inv(ir)%nrec, bufr_inv(ir)%nsubset,      &
          merge ('*','.', bufr_inv(ir)% file(1:n_source))
      end do
    endif

    !-----------------------------------------
    ! estimated cost to read observation types
    !-----------------------------------------
    cost_obs           =    1._wp
    cost_obs(OT_SYNOP) =    4._wp
    cost_obs(OT_AIREP) =    3._wp
    cost_obs(OT_SATOB) =    2._wp
    cost_obs(OT_DRIBU) =    3._wp
    cost_obs(OT_TEMP ) =  500._wp
    cost_obs(OT_PILOT) =   30._wp
    cost_obs(OT_SATEM) =   20._wp
    cost_obs(OT_PAOB ) =    1._wp
    cost_obs(OT_SCATT) =    2._wp
    cost_obs(OT_GPSRO) = 1000._wp
    cost_obs(OT_GPSGB) =  100._wp

    cost_obs(OT_RAD )  =  100._wp

    cost_pe            = 0._wp
    cost_file          = 0._wp
    subsetso           = 0
    !--------------------
    ! derive distribution
    !--------------------
    if (npe_read_obs < 1) then
      !----------------------------------
      ! read all observations from one PE
      !----------------------------------
      if (dace% pe==p_readbufr) bufr_inv% subsetl = bufr_inv% nsubset
    else

      !===========================================
      ! process dedicated files (from BUFR2NetCDF)
      !===========================================
      do is=1,n_source
        ir = source(is)% obstype
        if (ir<1) cycle
        !----------------------------------
        ! derive cost for reading this file
        !----------------------------------
        if (is>1) then
          subsets     = bufr_inv(ir)% subseto(is-1)
        else
          subsets     = 0
        endif
        subsetl       = bufr_inv (ir)% subseto(is)
        cost_file(is) = cost_file(is) + cost_obs (ir) * (subsetl-subsets)
        if (.not.source(is)% complete) cycle
        subsetso (ir) = subsetl
        codetype      = source(is)% codetype
        !-------------------------------
        ! files assigned to PEs a priori
        !-------------------------------
        pe            = -1
        !------------
        ! old version
        !------------
!       do if = 1,size(use_first)
!         if (use_first(if)==ir) pe = min(use_dedic(if),dace% npe-1)
!       end do
        !------------
        ! new version
        !------------
        do if = 1, size (obs_pe)
          if (obs_pe(if)% obstype >  0 .and. &
              obs_pe(if)% obstype /= ir      ) cycle
          if (codetype > 0 .and. obs_pe(if)% codetype(1) > 0 .and.       &
                                (obs_pe(if)% codetype(1) > codetype .or. &
                                 obs_pe(if)% codetype(2) < codetype     )) cycle
          pe = min (obs_pe(if)% pe, dace% npe-1)
          exit
        end do
        source(is)% pe = pe
        if (pe >= 0) then
          cost_pe(pe)=cost_pe(pe)+cost_file(is)
        endif
      end do
      !----------------------------------
      ! process files not assigned so far
      !----------------------------------
      call sort (cost_file, i_cost, -1, ier)
      do if=1,n_source
        is = i_cost (if)
        if (.not.source(is)% complete) cycle
!       ir = source(is)% obstype
        pe = source(is)% pe
        if (pe < 0) then
          pe             = sum(minloc(cost_pe))-1
          source(is)% pe = pe
          cost_pe(pe)=cost_pe(pe)+cost_file(is)
        endif
      end do
      !-----------------------------------
      ! print distribution of NetCDF files
      !-----------------------------------
      if(dace% lpio) then
        write(6,'(a)') repeat('-',79)
        write(6,'()')
        write(6,'(a)') '  Distribution of NetCDF-files for input'
        write(6,'()')
      end if

      if(dace% lpio .and. count(source(:)% complete)>0) then
        write(6,'(a)') '  # file-name                                         &
                   &pe type     subsets(start,end,total)      cost     cost/pe'
        write(6,'(a)')
        do is=1,n_source
          if (.not.source(is)% complete) cycle
          pe = source(is)% pe
          ir = source(is)% obstype
          oname   = ''
          subsets = 0
          subsetl = source(is)% entries
          if (ir > 0) then
            oname = obstyp(ir)%name
            if (is>1) subsets   = bufr_inv(ir)% subseto(is-1)
            subsetl             = bufr_inv(ir)% subseto(is)
          endif
          write(6,'(i3,1x,a48,i4,1x,a8,3i8,2f12.0)')           &
            is, source(is)% file, pe, oname, subsets, subsetl, &
            subsetl-subsets, cost_file(is), cost_pe(pe)
        end do
        write(6,'()')
      endif

      if(dace% lpio .and. count(.not. source(:)% complete)>0) then
        write(6,'(a)') '  # file-name                                         &
                   &pe type     subsets(start,end,total)      cost     cost/pe'
        write(6,'(a)')
        do is=1,n_source
          if (source(is)% complete) cycle
          ir = source(is)% obstype
          oname   = ''
          subsets = 0
          subsetl = source(is)% entries
          if (ir > 0) then
            oname = obstyp(ir)%name
            if (is>1) subsets   = bufr_inv(ir)% subseto(is-1)
            subsetl             = bufr_inv(ir)% subseto(is)
          endif
          write(6,'(i3,1x,a48,3x,"*",1x,a8,3i8,2f12.0)')          &
            is, source(is)% file, oname, subsets, subsetl,        &
            subsetl-subsets, cost_file(is), cost_file(is)/dace% npe
        end do
        write(6,'()')
      end if

      cost_obs_tot = (bufr_inv% nsubset-subsetso) * cost_obs
      cost_total   = sum (cost_obs_tot)
      cost_pe      = 0._wp
!     cost_file    = 0._wp
      !----------------------------------------------
      ! loop over obstypes, first most expensive ones
      !----------------------------------------------
      l_dedic = .false.
      do ir = 1, size (obs_pe)
         if (obs_pe(ir)% obstype >  0 .and. &
             obs_pe(ir)% pe      > -1 ) l_dedic(obs_pe(ir)% obstype) = .true.
      end do

      if    = 1
      first = .true.
      do
        !----------------------------
        ! settings to read on all PEs
        !----------------------------
        pe1 = 0
        pel = dace% npe-1
        !---------------------------------------------------------
        ! special handling as requested in 'use_first','use-dedic'
        !---------------------------------------------------------
        if (first) then
          if (if>size(use_first)) then
            !----------------------------------------------
            ! set start values for round robin distribution
            !----------------------------------------------
            first = .false.
            ir    = 1
          else
            ir = use_first(if)
            if (use_dedic(if)>=0) then
              !-------------------------------------------------
              ! settings to read this report type on one PE only
              !-------------------------------------------------
              pe1 = min(use_dedic(if),dace% npe-1)
              pel = pe1
            endif
            if = if + 1
          endif
        else
          ir = ir + 1
        endif
        if (.not. first .and. any (use_first==ir)) cycle
        if (ir > n_ot) exit
        pen = max (pel - pe1 + 1, 1)
        if (npe_read_obs > dace% npe) then
          !-----------------------------
          ! use all PEs for all obstypes
          !-----------------------------
          cost_mean_pe = (sum(cost_pe) + cost_obs_tot(ir)) / pen
        else
          !-------------------------------
          ! use dedicated PEs for obstypes
          !-------------------------------
          cost_mean_pe = cost_total / pen
        endif
        subsets = subsetso (ir)
        do pe = pe1, pel
          cost_left = max (0._wp, cost_mean_pe - cost_pe(pe))
          subsetn = ceiling (cost_left / cost_obs(ir))
          subsetl = min (subsets+subsetn, bufr_inv(ir)% nsubset)
          if (pe == pel) subsetl = bufr_inv(ir)% nsubset
          subsetn = subsetl - subsets
          cost_pe(pe) = cost_pe(pe) + subsetn * cost_obs(ir)
          if (pe == dace% pe) then
            bufr_inv(ir)% subsets = subsets
            bufr_inv(ir)% subsetl = subsetl
          endif
          subsets = subsetl
        end do
      end do
    endif
    !-----------------------------
    ! determine files used on a PE
    !-----------------------------
    source% used = .false.
    do ir = 1, n_ot
      if (bufr_inv(ir)% subsetl > bufr_inv(ir)% subsets) then
        do if = 1, n_source
          source(if)% used = source(if)% used .or. bufr_inv(ir)% file(if)
        end do
      endif
    end do
    !=========
    ! printout
    !=========
    if (netcdf_verb > 0) then
     if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'(a)')
      write(6,'(a)') '  Distribution of BUFR reports for input'
      write(6,'(a)')
      write(6,'(a,10x,100a8)') '    pe  ',obstyp% name
      write(6,'(a)')
     endif
     if (any (bufr_inv% subsetl > bufr_inv% subsets)) then
      call nextline
      write(oline(iol),'(i6,a,100i8)') dace% pe,'  start : ', bufr_inv% subsets
      call nextline
      write(oline(iol),'(i6,a,100i8)') dace% pe,'  end   : ', bufr_inv% subsetl
      call nextline
      write(oline(iol),'(i6,a,100i8)') dace% pe,'  used  : ',&
                                   bufr_inv% subsetl - bufr_inv% subsets
      call nextline
      write(oline(iol),'(i6,a,100l8)') dace% pe,'        : ',&
                                   bufr_inv% subsetl > bufr_inv% subsets
      call nextline
     end if
     call flush_buf
    end if

    if (netcdf_verb > 0) then
     if (dace% lpio) then
      write(6,'(a,10x,100a8)') '  file  ',obstyp% name
      write(6,'()')
      write(6,'(i6,a,100i8)')   1 ,' subsets: ',bufr_inv% subseto (1)
      write(6,'()')
      do if = 2, n_source
        write(6,'(i6,a,100i8)') if,' offset : ',bufr_inv% subseto (if-1)
        write(6,'(i6,a,100i8)') if,' subsets: ',bufr_inv% subseto (if) - &
                                                bufr_inv% subseto (if-1)
        write(6,'()')
      end do
      write(6,'()')
      write(6,'(a,6x,99i10)') '  file  ', (       if    , if=1, k_source/10)
      write(6,'(a,6x,999i1)') '        ', (modulo(if,10), if=1, k_source   )
      write(6,'("  processor")')
     endif
     call nextline
     write(oline(iol),'(i6,a,999a1)') dace% pe,'      : ', &
                                      merge ('*','.', source(1:n_source)% used)
     call flush_buf
    end if

  end subroutine set_input_boxes
!------------------------------------------------------------------------------
  subroutine set_cost_model (mode, cost_obs)
    !--------------------------------------------
    ! Set relative cost of observation processing
    ! for load balancing over processors.
    !--------------------------------------------
    integer, intent(in)  :: mode
    integer, intent(out) :: cost_obs(:)

    integer :: cost_model

    if (size (cost_obs) < n_ot) &
         call finish ("set_cost_model","cost_obs too small")

    cost_model = mode

    select case (method)
    case default
       select case (mode)
       case (0)
          cost_model = COST_OBS_READ
       case (1)
          cost_model = COST_OBS_FG
       case (2)
          cost_model = COST_OBS_ASS
       case default
          write(0,*) "set_cost_model: mode, method =", mode, trim (method)
          call finish ("set_cost_model","unsupported cost model")
       end select
    case ("GMESTAT","VERI_ENS")
       cost_model = COST_OBS_VERI
    end select

    cost_obs(:) = 1

    select case (cost_model)
    case (COST_OBS_READ)
       !-----------------------------------------
       ! estimated cost to read observation types
       !-----------------------------------------
       cost_obs(OT_SYNOP) =    4
       cost_obs(OT_AIREP) =    3
       cost_obs(OT_SATOB) =    2
       cost_obs(OT_DRIBU) =    3
       cost_obs(OT_TEMP ) =  500
       cost_obs(OT_PILOT) =   30
       cost_obs(OT_SATEM) =   20
       cost_obs(OT_PAOB ) =    1
       cost_obs(OT_SCATT) =    2
       cost_obs(OT_GPSRO) = 1000
       cost_obs(OT_GPSGB) =  100
    case (COST_OBS_FG)
       !-----------------------------------------------
       ! estimated cost to run tangent linear operators
       !-----------------------------------------------
       cost_obs(OT_RAD  ) =    3
       cost_obs(OT_GPSRO) =  100
       cost_obs(OT_GPSGB) =  100
    case (COST_OBS_ASS)
       !-----------------------------------------------
       ! estimated cost of observations in assimilation
       !-----------------------------------------------
       cost_obs(OT_RAD  ) =    2
       cost_obs(OT_GPSRO) =   20
       cost_obs(OT_GPSGB) =   20
    case (COST_OBS_VERI)
       !--------------------------------------------
       ! estimated cost to run observation operators
       !--------------------------------------------
       cost_obs(OT_RAD  ) =    3
       cost_obs(OT_GPSRO) = 5000
       cost_obs(OT_GPSGB) =  100
    end select

  end subroutine set_cost_model
!------------------------------------------------------------------------------
  subroutine set_veri_boxes (obs, pes, n_slot, mode, cosmope, atm)
  !=========================================================
  ! associate observations with boxes for verification (MEC)
  !=========================================================
  type (t_obs) ,intent(inout)        :: obs (:) ! observation data type
  integer      ,pointer              :: pes (:) ! processor / boxes
  integer      ,intent(in)           :: n_slot  ! number of time slots
  integer      ,intent(in)           :: mode    !
  logical      ,intent(in)           :: cosmope ! glue COSMO observations to PE
  type (t_atm) ,intent(in) ,optional :: atm     ! atmospheric reference state

    integer          :: ntri
    integer          :: nobs_fg ! approx.number obs/box in fg-scan
    integer          :: nctr    ! tot.number of cntr.params.
!   integer          :: m_slot
!   integer ,pointer :: pes_ (:)
!   integer          :: i, ib, is
    !------------------------------------------
    ! no specific handling, just call set_boxes
    !------------------------------------------
    nullify  (pes)
    ntri       =    0
    nobs_fg    = 5000
    obs(1)% pe = dace% pe
    call set_inbx  (obs)
    if (cosmope) then
      call set_cosmo_boxes (obs, ntri, nctr, pes, atm)
    else
      call set_boxes (obs, ntri, nctr, pes, .false., nobs_fg, 1, .false., &
                      STAT_OBS_ONLY)
    endif
    select case (mode)
    case (1:)
      !---------------------------
      ! single slot, nothing to do
      !---------------------------
!   case (-1:0)
!     !------------------------------------------------------
!     ! setup for fgat, choose next time slice or interpolate
!     !------------------------------------------------------
!     if (ntri/=dace% npe) call finish('set_veri_boxes','ntri/=dace% npe')
!     m_slot = n_slot + mode
!     allocate (pes_ (ntri*m_slot))
!     do i = 0, m_slot - 1
!       pes_(i*ntri+1:(i+1)*ntri) = pes
!     end do
!     do ib=1,size(obs)
!       if (obs(ib)% pe /= dace% pe) cycle
!       do is = 1, obs(ib)%n_spot
!         i = obs(ib)% spot(is)% i_time
!         if (i<1.or.i>m_slot) call finish('set_veri_boxes','i<1.or.i>m_slot')
!         obs(ib)%   spot(is)% fgbx% box =            &
!           obs(ib)% spot(is)% fgbx% box + (i-1) * ntri
!       end do
!     end do
!     ntri = ntri * m_slot
    case default
!     call finish ('set_veri_boxes','invalid mode')
    end select

  end subroutine set_veri_boxes
!------------------------------------------------------------------------------
  subroutine set_cosmo_boxes (obs, nboxes, nctr, pes, atm)
  !=================================================================
  ! distribute observations over boxes and processors
  ! for COSMO observations (with pre-assigned PEs, 1 box per PE)
  ! Output:
  ! obs% spot% fgbx% pe  =  processor_element_number (0..dace% npe-1)
  ! obs% spot% fgbx% box =  box_number               (1..
  !=================================================================
  type (t_obs) ,intent(inout)        :: obs(:)  ! observation data type
  integer      ,intent(inout)        :: nboxes  ! number of boxes
  integer      ,intent(out)          :: nctr    ! number of control variables
  integer      ,pointer              :: pes (:) ! processor / boxes
  type (t_atm) ,intent(in) ,optional :: atm     ! atmospheric reference state

    integer          :: ib  ! index variable (box)
    integer          :: is  ! index variable (spot)
    integer          :: pe  ! index variable (PE)
    integer ,pointer :: ijdp (:)
    !--------------------
    ! simply 1 box per PE
    !--------------------
    nboxes = dace% npe
    allocate (pes (nboxes))
    pes = (/(pe,pe=0,dace% npe-1)/)
    nctr = 0
    do ib=1,size(obs)
      if (obs(ib)% pe /= dace% pe) cycle
      if (.not.associated(obs(ib)%spot)) allocate (obs(ib)%spot(0))
      do is = 1, obs(ib)%n_spot
        ijdp => obs(ib)% spot(is)% col% h% ijdp
        if (present (atm)) then
          pe = atm% grid% marr (1,ijdp(1),ijdp(2),ijdp(3))
        else
          pe = ijdp(4)
        endif
        obs(ib)% spot(is)% fgbx% pe  = pe
        obs(ib)% spot(is)% fgbx% box = pe + 1
        nctr = nctr + 1
        if(pe < 0) call finish('set_cosmo_boxes','pe < 0')
      end do
    end do

  end subroutine set_cosmo_boxes
!------------------------------------------------------------------------------
  subroutine set_boxes (obs, nboxes, nctr, pes, freeio, nobs_box, iscan, &
                        lsort, stat_ok)
  !==================================================
  ! distribute observations over boxes and processors
  ! associate boxes with 'spots'
  !
  ! Output:
  ! obs% spot% pre1% pe  =  processor_element_number (0..dace% npe-1)
  ! obs% spot% pre1% box =  box_number               (1..
  !==================================================
  type (t_obs)    ,intent(inout) :: obs(:)      ! observation data type
  integer         ,intent(inout) :: nboxes      ! number of boxes
  integer         ,intent(out)   :: nctr        ! number of control variables
  integer         ,pointer       :: pes (:)     ! processor / boxes
  logical         ,intent(in)    :: freeio      ! no computations on I/O PE
  integer         ,intent(in)    :: nobs_box    ! approx.no. of observ./box
  integer         ,intent(in)    :: iscan       ! sort for 1:fg 2:psas
  logical         ,intent(in)    :: lsort       ! redistribute boxes over PEs
  integer         ,intent(in)    :: stat_ok     ! lowest status to keep

    !----------------
    ! local variables
    !----------------
    integer               :: ib      ! index variable (box, in FG scan)
    integer               :: is      ! index variable (spot in FG scan)
    integer               :: i1, in  ! index range
    integer               :: i       ! index in list of t_bx info
    integer               :: j       ! index in list of t_bx info (sorted)
    integer               :: nbox    ! index variable (box, in PSAS scan)
    integer               :: npe     ! index variable (PE,  in PSAS scan)
    integer               :: nspota  ! number of spots accepted
    integer(i8)           :: nctr2   ! matrix size (long integer)
    integer               :: ndgftot ! total no. observations
    integer               :: ndgf1   ! no. observations distributed so far
    integer               :: ndgf2   ! no. observations left for current box
    integer               :: nprocs  ! number of processors to be used
    integer               :: pio     ! index of processor not to use
    type (t_icol)         :: col     ! 'column' type to hold coordinates
    type (t_bx)  ,pointer :: ibx(:)  ! local list of t_bx info
    type (t_bx)  ,pointer :: obx(:)  ! global list of t_bx info (sorted)
    type (t_bx)  ,pointer :: bx      ! pointer to current element in list
    integer               :: nibx    ! no. reports left in local list
    integer               :: nobx    ! no. reports in global list
    integer               :: ns      ! no. reports in current box
    type(t_spot) ,pointer :: si      ! pointer to current report
    integer               :: nibxs   ! global sum of nibox
!   real(wp)              :: d2_prev ! Distance of prev. obs to ref. point
    real(wp),allocatable  :: box_x (:,:)
    integer ,allocatable  :: ip    (:)
    integer ,allocatable  :: newbox(:)
    integer               :: bpp       ! Boxes / proc.
    integer               :: nboxr     ! # of remaining boxes
    integer               :: ier       ! Error flag
    integer               :: np1,ib1
!   integer               :: nprocr,ib2,i2(2)
!   real(wp),allocatable  :: d     (:), d2(:,:)
    logical               :: lastd     ! last observation with same distance
    integer               :: mboxes    ! Actual number of non-empty boxes
    integer               :: nproc_use ! # of PEs used to determine # of boxes
    integer               :: cost_tot  ! Total cost of observation processing
    integer               :: cost_box  ! Avg. cost per box (goal)
    integer               :: cost      ! Accumulated cost
    integer               :: cost0     ! Accumulated cost (previous boxes)
    integer               :: cost2     ! Cost goal for current box
    integer               :: cost_obs (n_ot)       ! Observation cost (norm.)
    integer               :: cost_pe(0:dace% npe-1)! Total cost for each pe
    integer               :: idx_pe (dace% npe)    ! Index vector for sort
    integer ,allocatable  :: cost_bx(:)            ! Total cost for each box
    integer ,allocatable  :: nobs_bx(:,:)          ! Obs. count per type,box
    logical               :: crit (2)
    !---------------------------------------
    ! no processing on I/O processor element
    !---------------------------------------
    if (freeio .and. dace% npe > 1) then
      nprocs = dace% npe - 1
      pio    = dace% pio
    else
      nprocs = dace% npe
      pio    = dace% npe + 1
    endif
    !-------------------------------
    ! fill in data for parallel sort
    !-------------------------------
    nullify (ibx)
    nibx = 0
    do ib=1,size(obs)
      if (obs(ib)% pe /= dace% pe) cycle
      ns = obs(ib)%n_spot
      if (.not.associated(obs(ib)%spot)) allocate (obs(ib)%spot(0))
      call fill_bx (ibx, nibx, obs(ib), ib, &
                    mask=obs(ib)%spot(1:ns)% use% state >= stat_ok)
      !-----------------------------
      ! reject observations not used
      !-----------------------------
      obs(ib)% spot% pre1% box = -1
      obs(ib)% spot% pre1% pe  = -1
      if (iscan==1)                    cycle
      do i=1,ns
        si => obs(ib)%spot(i)
        if (si% use% state >= stat_ok) cycle
        select case (int(si% use% state))
        case (:STAT_DISMISS)
          ! don't write_report
        case (STAT_OBS_ONLY)
          call write_report (si)
          oonly_stat   (si%use% check, si%hd% obstype) = &
            oonly_stat (si%use% check, si%hd% obstype) + 1
        case (STAT_REJECTED)
          call write_report (si)
          reje_stat   (si%use% check, si%hd% obstype) = &
            reje_stat (si%use% check, si%hd% obstype) + 1
        case (STAT_PASSIVE, STAT_PAS_REJ)
          call write_report (si)
          pass_stat   (si%use% check, si%hd% obstype) = &
            pass_stat (si%use% check, si%hd% obstype) + 1
        case default
          call write_report (si)
        end select
      end do
    end do
    if (.not.associated (ibx)) allocate (ibx(0))
    !---------------------------------------------
    ! Set relative cost for observation processing
    ! (tunable for load balancing over processors)
    !---------------------------------------------
    call set_cost_model (iscan, cost_obs)
    !------------------------------------------------
    ! sort according to distance from arbitrary point
    !------------------------------------------------
    col% c% dlon = 50._wp
    col% c% dlat =  0._wp
    call set_xuv (col)
    nctr = 0
    cost = 0
    do i=1,nibx
      is =  ibx(i)% box% spot
      ib =  ibx(i)% box% box
      si => obs(ib)% spot(is)
      ibx(i)% x     = si% col% c% x
      ibx(i)% d1    = rearth * sum((ibx(i)%x(:) - col% c% x)**2)
      ibx(i)% d2    = ibx(i)% d1
      select case (iscan)
      case (1)
        !-----------------------
        ! FG: cost = no. columns
        !-----------------------
        ibx(i)% count = 3
        if (associated (si% imcol)) ibx(i)% count = size(si% imcol)
!       !------------------------------
!       ! FG: cost = no. arguments to H
!       !------------------------------
!       ibx(i)% count = si% i% n
      case (2)
        !------------------------------
        ! PSAS: cost = no. observations
        !------------------------------
!       ibx(i)% count = si% o% n
        i1 = si% o% i + 1
        in = si% o% i + si% o% n
        ibx(i)% count = count (obs(ib)% body(i1:in)% use% state >= stat_ok)
!       if (si% hd% obstype == OT_RAD) ibx(i)% count = 2 * si% o% n
      end select
      ibx(i)% otyp = si% hd% obstype
      ibx(i)% cost = ibx(i)% count * cost_obs(si% hd% obstype)
      cost = cost  + ibx(i)% cost
      nctr = nctr  + ibx(i)% count
    end do
    nspota   = p_sum(nibx)
    nctr     = p_sum(nctr)
    cost_tot = p_sum(cost)

    allocate (ip(0:dace% npe-1))
    call p_gather (nibx, ip, root=dace% pio)
    if (dace% lpio) then
       write (6,'()')
       write (6,*) "set_boxes : nibx vs. processors:"
       if (dace% npe > 8) then
          write (6,'(12i8)') ip(:)
       else
          write (6,'(16i8)') ip(:)
       end if
       write (6,'()')
       write (6,*) "set_boxes : min, max, avg, rms =", &
            minval (ip), maxval (ip), &
            nint (      sum (real (ip,wp)   )/dace% npe), &
            nint (sqrt (sum (real (ip,wp)**2)/dace% npe))
       write (6,'()')
    end if
    deallocate (ip)

    if (nproc_repro > 0) then
      !--------------------------------------------------------
      ! Use given # of fictitious PEs for repeatable # of boxes
      !--------------------------------------------------------
      nproc_use = nproc_repro
    else
      nproc_use = nprocs
    end if
    !------------------------------------
    ! modify number of boxes if necessary
    !------------------------------------
    if (nboxes <=0) then
      if (nobs_box == 0) then
        nboxes = nproc_use
      else
        nboxes = nint (nctr / real(nobs_box,wp))
        nboxes = nint (nboxes / real(nproc_use,wp)) * nproc_use
        nboxes = max  (nboxes, nproc_use)
      end if
    endif
    nctr2 = nctr; nctr2 = nctr2 * nctr2
    if (dace% lpio) then
      select case (iscan)
      case (1)
        print *,'set_boxes : estimated total cost   = ',cost_tot
!       print *,'set_boxes : estimated total cost   = ',nctr
        print *,'set_boxes : number of spots used   = ',nspota
        print *,'set_boxes : number of boxes        = ',nboxes
        print *,'set_boxes : nproc_use              = ',nproc_use
        print *
      case (2)
        print *,'set_boxes : estimated total cost   = ',cost_tot
        print *,'set_boxes : number of deg. freedom = ',nctr
        print *,'set_boxes : number of spots used   = ',nspota
        print *,'set_boxes : number of boxes        = ',nboxes
        print *,'set_boxes : nproc_use              = ',nproc_use
        print *,'set_boxes : size of osas matrix    = ',nctr2
        print *
      end select
    endif

    ndgftot = nctr
    ndgf1   = 0
    cost    = 0
    nullify (obx)

    allocate (pes      (        nboxes))
    allocate (box_x    (3,      nboxes))
    allocate (cost_bx  (        nboxes))
    allocate (nobs_bx  (n_ot, 0:nboxes))
    nobs_bx = 0
    cost_bx = 0
    cost_pe = 0
    pes = -1
    !-----------------------
    ! loop over boxes to use
    !-----------------------
!if (iscan==1) call stop_time ("set_boxes: sort")
    mboxes = 0
    do nbox=1,nboxes
      ndgf2 = real (ndgftot,wp) / nboxes * nbox
      if (nbox==nboxes) ndgf2 = ndgftot
!     npe   = mod (nboxes-1, nprocs) + 1
      npe = (nbox-0.5_wp) / nboxes * nprocs
      if (npe >= pio) npe = npe + 1
      pes (nbox) = npe
      !-----------------------
      ! derive next box center
      !-----------------------
      col% c% x = nextloc (ibx, nibx)
      box_x (:,nbox) = col% c% x(:)
      do i = 1, nibx
        ibx(i)% d2 = rearth * sum((ibx(i)%x(:) - col% c% x)**2)
      end do
      !-----------------------------------------------------
      ! prepare list of neighbours in vicinity of box center
      !-----------------------------------------------------
      cost_box = nint (real (cost_tot-cost,wp) / (nboxes+1-nbox)) ! Cost/box
      cost2    = nint (real (cost_tot,wp)      /  nboxes * nbox)  ! Cost goal
      cost0    = cost
      call p_sort (obx, nobx, ibx, nibx, ndgf2-ndgf1, max(cost2-cost0, cost_box/2))

      !----------------------------------------------------------------------
      ! Fill box.  Try to keep observations at same distance in the same box.
      !----------------------------------------------------------------------
!     d2_prev  = 0._wp
      crit     = .false.
      do i=1,nobx
        bx    => obx(i)
        lastd = .true.
        if (i<nobx) then
          if (obx(i+1)%d2 == bx% d2) lastd = .false.
        endif
        !----------------
        ! update counters
        !----------------
!       d2_prev = bx% d2
        ndgf1 = ndgf1 + bx% count
        cost  = cost  + bx% cost
        nobs_bx(bx% otyp,nbox) = nobs_bx(bx% otyp,nbox) + bx% count
        crit  = .false.
        !--------------------------
        ! fill box on respective PE
        !--------------------------
        if (bx% box% pe == dace% pe) then
          ib = bx% box% box
          is = bx% box% spot
          j  = bx%      isrc
          ibx(j)% valid = .false.
          obs(ib)% spot(is)% pre1% box = nbox
          obs(ib)% spot(is)% pre1% pe  = npe
        endif
        !-----------------------
        ! check if box is 'full'
        !-----------------------
        if (ndgf1 >= ndgf2) then
           !----------------------------
           ! Degrees of freedom exceeded
           !----------------------------
           crit(1) = .true.
           if (lastd) exit
        else if (cost > cost2) then
           !-----------------------------------------------------------------
           ! Cost goal exceeded, make sure there is at least 1 obs in the box
           !-----------------------------------------------------------------
!          if (cost > cost0) then
!             if (bx% d2 > d2_prev) exit
!          end if
           !------------------------------------------------------
           ! Cost goal exceeded, try to avoid excessive imbalances
           !------------------------------------------------------
           if (cost > cost0 + cost_box/2) then
             crit(2) = .true.
             if (lastd) exit
           end if
        end if
        crit = .false.
      end do
      cost_bx(nbox) = cost - cost0
      cost_pe(npe)  = cost_pe(npe) + cost_bx(nbox)
      !----------------------------
      ! rearrange local report list
      !----------------------------
      call  release_bx (ibx, nibx)
      !---------
      ! printout
      !---------
      nibxs = p_sum(nibx)
      if (dace% lpio) then
        if (nbox==1) then
          write(6,*) '     box       pe         cost          observations   reports'
          write(6,*) '                     proc.     goal    proc.     goal     left'
        end if
        write(6,'(8i9,2l2)') nbox,    npe,  cost ,  cost2,  ndgf1,  ndgf2,  nibxs, nobx, crit
      endif
      !-----------------------------
      ! Determine last non-empty box
      !-----------------------------
      if (nibxs > 0 .or. nobx > 0) mboxes = nbox
    end do
    if (ndgf1 < ndgftot .or. cost < cost_tot) then
      if (dace% lpio) then
        write(0,'(A)') repeat('-',79)
        write(0,'(A)') "set_boxes: WARNING: observations may have been lost!"
        write(0,'(2(A,i10))') "Degrees of freedom:",ndgf1," out of",ndgftot
        write(0,'(2(A,i10))') "Total cost (est.) :",cost, " out of",cost_tot
        write(0,'(A)') repeat('-',79)
      end if
    end if
    !------------------------------------
    ! Show distribution of load estimates
    !------------------------------------
    call print_dist ()

    !---------------------------------------------------------------------
    ! Reorder boxes to achieve approximate balance of cost over processors
    !---------------------------------------------------------------------
    if (lsort .and. nboxes > nprocs .and. .not. freeio) then
      if (dace% lpio) then
         write (6,'()')
         write (6,'(A)') " Shuffling boxes"
         if (netcdf_verb > 0) then
            write(6,'()')
            write(6,'(A)') "  newbox  oldbox    cost  new_pe"
         end if
      end if
      allocate (ip     (mboxes))
      allocate (newbox (mboxes))
      bpp = ceiling (real (mboxes,wp) / nprocs)     ! Max. # of boxes/proc.
      !------------------------
      ! Sort by decreasing cost
      !------------------------
      call sort (cost_bx(1:mboxes), ip, -1, ier)
      pes(1:mboxes) = -1
      newbox  = -1
      cost_pe = 0
      idx_pe  = (/ (i, i=0,dace% npe-1) /)   ! Initial processor assignment
      ib      = 0
      nboxr   = mboxes
      do i = 1, bpp
         ib1 = min (nboxr, nprocs)
         do np1 = 1, ib1
            ib = ib + 1
            j              = idx_pe(np1)    ! Processor #
            newbox(ip(ib)) = ib
            pes   (ip(ib)) = j
            cost_pe(j)     = cost_pe(j) + cost_bx(ip(ib))
            if (dace% lpio .and. netcdf_verb > 0) then
               write(6,'(4i8)') ib,ip(ib),cost_bx(ip(ib)),j
            end if
         end do
         nboxr = nboxr - ib1
         if (nboxr > 0) then
            !-----------------------------------
            ! Sort processors by increasing load
            !-----------------------------------
            call sort (cost_pe, idx_pe, 1, ier)
            idx_pe = idx_pe-1
         end if
      end do
      if (any (newbox < 1) .or. any (newbox > mboxes)) then
         if (dace% lpio) write (0,*) "set_boxes: newbox =", newbox
         call finish ("set_boxes","invalid box remapping")
      end if
      if (any (pes < 0) .or. any (pes > nprocs-1)) then
         if (dace% lpio) write (0,*) "set_boxes: pes =", pes
         call finish ("set_boxes","invalid PE assignment after box reordering")
      end if
      if (sum (cost_pe) /= cost_tot) then
         if (dace% lpio) write (0,*) "cost before/after:",cost_tot,sum(cost_pe)
         call finish ("set_boxes","bad total cost after box reordering")
      end if
      do ib = 1, size(obs)
        if (associated(obs(ib)% spot)) then
          do is = 1, obs(ib)% n_spot
            if (obs(ib)% spot(is)% pre1% box > 0) then
              obs(ib)%spot(is)%pre1% pe  = pes   (obs(ib)%spot(is)%pre1% box)
              obs(ib)%spot(is)%pre1% box = newbox(obs(ib)%spot(is)%pre1% box)
            endif
          end do
        endif
      end do
      ip          = pes(1:mboxes)
      pes(newbox) = ip
      !----------------------------------------
      ! Adjust observation count to permutation
      !----------------------------------------
      do i = 1, n_ot
         ip                = nobs_bx(i,1:mboxes)
         nobs_bx(i,newbox) = ip
      end do
      if (dace% lpio) then
         write (*,'()')
         write (*,'(a)') " After box reordering:"
      end if
      call print_dist (full=.false.)
      deallocate (ip, newbox)
    end if


!   if (lsort .and. iscan == 2 .and. nboxes > nprocs) then
!     allocate (d2    (mboxes,mboxes))
!     allocate (d     (mboxes))
!     allocate (ip    (mboxes))
!     allocate (newbox(mboxes))
!     pes(1:mboxes) = -1
!     newbox = -1
!     nboxr  = mboxes
!     nprocr = nprocs
!     ib     = 0
!     do np1 = 1,nprocs
!       bpp = nboxr / nprocr
!       do ib1=1,mboxes
!         do ib2=1,mboxes
!           d2 (ib1,ib2) = sum((box_x(:,ib1) - box_x(:,ib2))**2)
!           if (d2 (ib1,ib2)==0._wp) d2(ib1,ib2)= 100._wp
!         end do
!       end do
!       i2 = minloc(d2(:,:))
!       do nbox=1,mboxes
!           d (nbox) = sum((box_x(:,nbox) - box_x(:,i2(1)))**2)
!       end do
!       call sort (d, ip, 1, ier)
!       do i = 1,bpp
!         ib = ib + 1
!         newbox  (ip(i)) = ib
!         pes     (ip(i)) = np1-1
!         box_x (:,ip(i)) = 100._wp
!       end do
!       nboxr  = nboxr - bpp
!       nprocr = nprocr - 1
!     enddo
!     if (any (newbox < 1) .or. any (newbox > mboxes)) then
!        write (0,*) dace% pe, "set_boxes: newbox =", newbox
!        call finish ("set_boxes","invalid box remapping")
!     end if
!     if (any (pes < 0) .or. any (pes > nprocs)) then
!        write (0,*) dace% pe, "set_boxes: pes =", pes
!        call finish ("set_boxes","invalid PE assignment after box reordering")
!     end if
!     do ib = 1, size(obs)
!       if (associated(obs(ib)% spot)) then
!         do is = 1, obs(ib)% n_spot
!           if (obs(ib)% spot(is)% pre1% box > 0) then
!             obs(ib)%spot(is)%pre1% pe  = pes   (obs(ib)%spot(is)%pre1% box)
!             obs(ib)%spot(is)%pre1% box = newbox(obs(ib)%spot(is)%pre1% box)
!           endif
!         end do
!       endif
!     end do
!     ip          = pes(1:mboxes)
!     pes(newbox) = ip
!     deallocate (ip, newbox, d, d2)
!   endif

    deallocate (box_x)
    deallocate (ibx,obx)

    !---------------------------------------------
    ! FG scan: set obs(:)% spot(:)% fgbx, not pre1
    !---------------------------------------------
    if (iscan==1) then
      do ib = 1, size(obs)
        obs(ib)% spot% fgbx = obs(ib)% spot% pre1
        obs(ib)% spot% pre1% box = -1
        obs(ib)% spot% pre1% pe  = -1
      end do
    endif

!if (iscan==1) call stop_time ("set_boxes: write_pending")
    !----------------------
    ! write pending reports
    !----------------------
    call write_pending

  contains

    subroutine print_dist (full)
      logical, intent(in), optional :: full

      integer :: j, nbox
      logical :: ful
      logical :: mask (n_ot)       ! Obstype present?

      ful = .true.; if (present (full)) ful = full
      if (dace% lpio) then
        nobs_bx(:,0) = sum (nobs_bx(:,1:mboxes),dim=2)
        mask         = nobs_bx(:,0) > 0
        if (ful) then
          write(6,'()')
          write(6,*) "Distribution of observation count vs. type over boxes"
          write(6,'()')
!#if defined (__SX__)
!         ! Workaround for bug in sxf90 revs. 400,430,440
!         if (count (mask) > 0) then
!            write(6,'(a,20a8)') "   Box    PE", adjustr (pack (obstyp% name, mask))
!         else
!            write(6,'(a,20a8)') "   Box    PE"
!         end if
!#else
          write(6,'(a,20a8)') "   Box    PE", adjustr (pack (obstyp% name, mask))
!#endif
          do nbox=1, mboxes
             write(6,'(2i6,20i8)') nbox, pes(nbox), pack (nobs_bx(:,nbox),mask)
          end do
          write(6,'(a,20i8)') " Total      ", pack (nobs_bx(:,0)   ,mask)
        end if
        write(6,'()')
        write(6,*) "Distribution of observation count vs. type over processors"
        write(6,'()')
!#if defined (__SX__)
!       ! Workaround for bug in sxf90 revs. 400,430,440
!       if (count (mask) > 0) then
!          write(6,'(a,20a8)') "    pe", adjustr (pack (obstyp% name, mask))
!       else
!          write(6,'(a,20a8)') "    pe"
!       end if
!#else
        write(6,'(a,20a8)') "    pe", adjustr (pack (obstyp% name, mask))
!#endif
        do j = 0, nprocs-1
          nobs_bx(:,0) = sum (nobs_bx(:,1:),dim=2,&
                              mask=spread (pes(:)==j,dim=1,ncopies=n_ot))
          write(6,'(i6,20i8)') j, pack (nobs_bx(:,0),mask)
        end do
        if (netcdf_verb > 0) then
          write (*,'()')
          write (*,'(A)') " Cost distribution"
          write (*,'()')
          write (*,'(a)') "    pe      cost"
          do j = 0, nprocs-1
             write (*,'(i6,i10)') j, cost_pe(j)
          end do
          write (*,'(a,i10)') " Total", sum (cost_pe)
          write (*,'()')
        end if
      end if

    end subroutine print_dist

  end subroutine set_boxes
!------------------------------------------------------------------------------
  subroutine p_sort (obx, nobx, ibx, nibx, maxcnt, maxcst)
  type(t_bx) ,pointer     :: obx(:) ! list of reports in box (gathered)
  integer    ,intent(out) :: nobx   ! number of reports in box
  type(t_bx) ,intent(in)  :: ibx(:) ! list of reports considered (distributed)
  integer    ,intent(in)  :: nibx   ! number of reports considered
  integer    ,intent(in)  :: maxcnt ! maximum count of reports wanted
  integer    ,intent(in)  :: maxcst ! maximum cost  of reports wanted

    integer    ,parameter   :: nbx = 50       ! number of ranges for pre-sort
    integer                 :: ip1(nibx)      ! index array
    integer                 :: cbx(nbx,3)     ! count/cost in pre-sort boxes
    real(wp)                :: edstep
    integer    ,allocatable :: ipo(:)
    integer                 :: i, j
    integer                 :: count
    integer                 :: ier
    type(t_bx) ,allocatable :: sndbx  (:)
    type(t_bx) ,allocatable :: rcvbx  (:)
    integer                 :: nsend1, nsend2, nrecv1, nrecv2
    real(wp)                :: d2_max
    integer                 :: k, ii            ! auxiliary indices for sort

    !---------------------------------------------
    ! pre-sort local reports into nbx classes
    ! with increasing distance from starting point
    !---------------------------------------------
    d2_max = maxval (ibx(1:nibx)% d2)
    d2_max = p_max (d2_max)
    edstep = nbx / (d2_max + 1._wp)
    cbx    = 0
    do i = 1, nibx
      ip1(i) = max (1, ceiling (ibx(i)% d2 * edstep))
    end do
    do i = 1, nibx
      j        = ip1(i)
      cbx(j,1) = cbx(j,1) + 1
      cbx(j,2) = cbx(j,2) + ibx(i)% count
      cbx(j,3) = cbx(j,3) + ibx(i)% cost
    end do
    cbx(:,2:3) = p_sum(cbx(:,2:3))
    do i = 2, nbx
      cbx(i,:) = cbx(i,:) + cbx(i-1,:)
    end do
    do i = 1, nbx
      j = i
      if (cbx(i,2) >= maxcnt) exit
      if (cbx(i,3) >  maxcst) exit
    end do
    !------------------------------------------
    ! option to sort only the last 'class' used
    ! currently not used: !!!
    !------------------------------------------
    nsend2 = cbx(j,1)
    nsend1 = 0
!!! if (j>1) nsend1 = cbx(j-1,1)

    nrecv1 = p_sum (nsend1)
    nrecv2 = p_sum (nsend2)
    allocate (sndbx(nsend2))
    allocate (rcvbx(nrecv2))
    allocate (ipo  (nrecv2))
    !-------------------------
    ! fill send buffer, gather
    !-------------------------
    k = 0
    do i = 1, nibx
      if (ip1(i) >= j) cycle
      k = k + 1
      sndbx(k) = ibx(i)
    end do
    do i = 1, nibx
      if (ip1(i) /= j) cycle
      k = k + 1
      sndbx(k) = ibx(i)
    end do
!!! call p_gather_bx (sndbx(       1:nsend1), rcvbx(       1:nrecv1), root=dace% pio)
    call p_gather_bx (sndbx(nsend1+1:nsend2), rcvbx(nrecv1+1:nrecv2), root=dace% pio)
    !------------
    ! global sort
    !------------
    if (dace% lpio) then
!!!    ipo(1:nrecv1) = [(i,i=1,nrecv1)]
       call sort (rcvbx(       1:nrecv2)% d2, ipo(       1:nrecv2), 1, ier) ! sort spots according to penalty
!!!    call sort (rcvbx(nrecv1+1:nrecv2)% d2, ipo(nrecv1+1:nrecv2), 1, ier) ! sort spots according to penalty
!!!    ipo (nrecv1+1:nrecv2) = ipo (nrecv1+1:nrecv2) + nrecv1
    end if
    !------------------------------------
    ! fill output buffer on I/O processor
    !------------------------------------
    if (dace% lpio) then
       count = 0
       nobx = nrecv2
       do i = 1, nrecv2
          j = ipo(i)
          count = count + rcvbx(j)% count
          if (count > maxcnt) then
             nobx = i
             exit
          endif
       end do
       !------------------------------------------------------------------
       ! For deterministic distribution of observations over boxes, append
       ! all observations collocated with the last one selected above
       !------------------------------------------------------------------
       if (nobx < nrecv2) then
          j = ipo(nobx)
          d2_max = rcvbx(j)% d2
          do i = nobx+1, nrecv2
             j = ipo(i)
             if (rcvbx(j)% d2 == d2_max) then
                count = count + rcvbx(j)% count
                nobx = i
             else
                exit
             end if
          end do
       end if
    end if
    call p_bcast (nobx, dace% pio)
    if (.not. associated (obx)) allocate (obx(nobx))
    if (size(obx) < nobx) then
       deallocate (obx)
       allocate (obx(nobx))
    endif
    if (dace% lpio) then
       do i = 1, nobx
          j = ipo(i)
          !------------------------------------------------------
          ! Fill output buffer, performing a final insertion sort
          ! on the selected observations:
          ! 1) increasing distance
          ! 2) decreasing no. of reports
          !------------------------------------------------------
          k = i
          if (i > 1) then
             do ii = i-1, 1, -1
                if (rcvbx(j)% d2    >  obx(ii)% d2  .or. &
                    rcvbx(j)% count <= obx(ii)% count    ) exit
                obx(ii+1) = obx(ii)
                k = ii
             end do
          end if
          obx(k) = rcvbx(j)
       end do
    end if
    !------------------
    ! Broadcast results
    !------------------
    call p_bcast_bx (obx, dace% pio)

    deallocate (rcvbx, ipo)
    deallocate (sndbx)

  end subroutine p_sort
!------------------------------------------------------------------------------
  subroutine fill_bx (ibx, nibx, obs, ibox, mask)
  type(t_bx)  ,pointer              :: ibx(:)   ! list of reports to extend
  integer     ,intent(inout)        :: nibx     ! number of reports in list
  type(t_obs) ,intent(in)           :: obs      ! observation variable
  integer     ,intent(in)           :: ibox     ! actual observation 'box'
  logical     ,intent(in) ,optional :: mask(:)  ! mask for reports to use

    integer              :: i, j, n
    type (t_bx) ,pointer :: tmp (:)

    if (obs% pe /= dace% pe) return

    n = obs% n_spot; if (present(mask)) n=count (mask)
    if (.not.associated(ibx)) then
      allocate (ibx(n))
    endif

    if (size(ibx) < nibx + n) then
      tmp => ibx
      allocate (ibx (nibx+n))
      ibx (1:nibx) = tmp (1:nibx)
      deallocate (tmp)
    endif

    j = 0
    do i=1,obs% n_spot
      if (present(mask)) then
        if(.not.mask(i)) cycle
      endif
      j = j + 1
      ibx(nibx+j)% box% pe   = dace% pe
      ibx(nibx+j)% box% box  = ibox
      ibx(nibx+j)% box% spot = i
      ibx(nibx+j)% isrc      = nibx+j
      ibx(nibx+j)% valid     = .true.
    end do
    nibx = nibx + n
  end subroutine fill_bx
!------------------------------------------------------------------------------
  subroutine release_bx (ibx, nbx)
  type(t_bx) ,intent(inout) :: ibx(:)
  integer    ,intent(inout) :: nbx
    integer :: i, j
    j = 0
    do i=1,nbx
      if (ibx(i)% valid) then
        j = j + 1
        if (j < i) then
          ibx(j) = ibx(i)
        endif
      endif
    end do
    do i=1,j
       ibx(i)% isrc = i
    end do
    ibx(j+1:)% valid = .false.
    nbx = j
  end subroutine release_bx
!------------------------------------------------------------------------------
  function nextloc (ibx, nbx) result (x)
  type(t_bx)  ,intent(in) :: ibx(:) ! list of reports
  integer     ,intent(in) :: nbx    ! length of list to consider
  real(wp)                :: x (3)  !
  !-------------------------------------------------------------
  ! search the location of the closest (to the initial location)
  ! report not used so far to the initial location
  !-------------------------------------------------------------
    integer    :: i(1)
    type(t_bx) :: bx
    type(t_bx) :: bxs (dace% npe)
    bx% d1 = huge(bx%d1)
    bx% x  = 0._wp
    if (any(ibx(1:nbx)% valid)) then
      i  = minloc (ibx(1:nbx)% d1, mask=ibx(1:nbx)% valid)
      bx = ibx(i(1))
    endif
    call p_gather_bx ((/bx/), bxs, root=dace% pio)
    if (dace% lpio) then
      i = minloc (bxs% d1)
      x = bxs(i(1))% x
    end if
    call p_bcast (x, dace% pio)
  end function nextloc
!==============================================================================
! preliminary allgather routine
!------------------------------
  subroutine p_allgather (sendbuf, recvbuf)
  type (t_bx)       ,INTENT(in)  :: sendbuf (:)
  type (t_bx)       ,INTENT(out) :: recvbuf (:)
    integer    :: sendcounts (dace% npe)
    integer    :: i
    integer    :: n
    type(t_bx) :: tmp (size(sendbuf)* dace% npe)
    n          = size(sendbuf)
    sendcounts = n
    do i = 0, dace% npe-1
      tmp(i*n+1:i*n+n) = sendbuf
    end do
    call p_alltoall_bx (tmp, recvbuf, sendcounts=sendcounts)
  end subroutine p_allgather
!==============================================================================
#define DERIVED type(t_bx)
#undef  MPI_TYPE
#define p_alltoall_DERIVED p_alltoall_bx
#include "p_alltoall_derived.incf"
#undef  DERIVED
!==============================================================================
#define DERIVED type(t_bx)
#undef  MPI_TYPE
#define p_gather_DERIVED p_gather_bx
#include "p_gather_derived.incf"
#undef  DERIVED
!==============================================================================
!------------------------
! Broadcast array of t_bx
!------------------------
#define DERIVED type(t_bx),dimension(:)
#define VECTOR
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_bx
#include "p_bcast.incf"
#undef  DERIVED
!==============================================================================
end module mo_boxes
