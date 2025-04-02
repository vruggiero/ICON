!
!+ Read observation data from BUFR file
!
MODULE mo_obs_bufr_dwd
!
! Description:
!   Read observation data from BUFR file.
!   Subroutine read_obs_bufr first has to be called in a scanning mode to
!   get an inventory of the BUFR files. Based on the inventory it is
!   decided which records will be read (on which processor element) in the
!   second call of the routine. The routine calls the observation type
!   specific routines in the respective modules.
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
!  thinning: implement new parameters 'instr', 'pass
! V1_20        2012-06-18 Harald Anlauf
!  Framework for arbitrary number of input feedback files
! V1_28        2014/02/26 Andreas Rhodin
!  comments added
! V1_29        2014/04/02 Andreas Rhodin
!  consistently use n_ot instead of n_obstype
! V1_44        2015-09-30 Harald Anlauf
!  move type t_bufr_inv to mo_t_obs to simplify dependencies
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Author: Andreas Rhodin MPIfM 2000      Original code (use EMOS library)
!         Andreas Rhodin DWD   2003-2008 use DWD BUFR3 routines
!=========================================================================
  !=============
  ! Modules used
  !=============
  use mo_exception, only: finish           ! abort routine
  use mo_time,      only: t_time,         &! time data type
                          init_time,      &! initialise time data type
                          operator(<),    &! compare times
                          operator(<=)     ! compare times
  use mo_run_params,only: p_readbufr,     &! PE to read BUFRs on
                          path_file        ! add path to a file name
  !--------------------------------
  ! acess to observation data types
  !--------------------------------
  use mo_t_obs,     only: t_obs,          &! observation data type
                          t_spot,         &! component of t_obs
                          t_head,         &! component of t_spot
                          bufr_pause,     &! wait (for CR) after record read
                          bufr_verb,      &! verbosity flag
                          derive_dbkz,    &! derive DBKZ if not present
                          source,         &! list   of Report source files
                          n_source,       &! number of Report source files
!                         m_source,       &! max no of Report source files
                          t_bufr_inv,     &! observation inventory data type
                          bufr_inv,       &! observation inventory data
                          release_mem      ! release unused memory
  use mo_wmo_tables,only: wmo_mnem,       &! get mnemonics for table A
                          WMO0_ECMWF       ! generating center
  use mo_dwd_tables,                      &! Datenbank-Kennziffer:
                  only: dbkz_mnem,        &! Datenbank-Kennziffer mnem.
                        DK_SATOB_2,       &! Satellitenbeob, SATOB, Section 2
                        DK_SATOB_3,       &! Satellitenbeob, SATOB, Section 3
                        DK_AMV_EUMETSAT,  &!
                        DK_AMV_GOES,      &!
                        DK_AMV_MODIS,     &!
!                       DK_AMV_MODIS_OLD, &!
!                       DK_AMV_MODIS_A,   &!
                        DK_AMV_FY_X,      &!
                        DK_AMV_MTV,       &!
                        DK_AMV_NOAA        !
  use mo_obstypes,only: t_obsid,          &! observation id table entry
                        obstype_dbkz,     &! derive obsids from dbkz
                        obstype_bufr       ! derive obsids from bufr type
  use mo_obs_tables,                      &!
                  only: idb_dbk,          &! index in table rept_stat
                        rept_char,        &! observation type characteristics
                        check_report_0,   &! simple checks on report validity
                        dism_report_0,    &! dismiss report
!                       print_rept_stat,  &! print observation type statistics
                        write_pending,    &! write pending dismissed reports
                        obstyp             ! table
  use mo_temp,    only: read_temp_bufr     ! read TEMP observation from BUFR
  use mo_synop,   only: read_synop_bufr    ! read SYNOP message from BUFR
  use mo_amv,     only: read_amv_bufr      ! read AMV   message from BUFR
  use mo_airep,   only: read_airep_bufr    ! read AIREP message from BUFR
  use mo_occ,     only: read_gpsro_bufr    ! read GPSRO message from BUFR
  use mo_t_use,   only: t_use,            &! status variable data type
                        STAT_DISMISS,     &! status: dismissed
                        STAT_PAS_REJ,     &! status: passive
                        CHK_NOIMPL         ! reason for dissmission
  use mo_test_obs,only: test_obs           ! test consistency of obs.datatype
  use mo_fdbk_tables,                     &!
                  only: n_ot,             &! number of observation types
                        OT_SYNOP,         &! report type identifyer
                        OT_SATOB,         &! report type identifyer
                        OT_TEMP,          &! report type identifyer
                        OT_PILOT,         &! report type identifyer
                        OT_AIREP,         &! report type identifyer
                        OT_DRIBU,         &! report type identifyer
                        OT_PAOB,          &! report type identifyer
                        OT_SCATT,         &! report type identifyer
                        OT_GPSRO           ! report type identifyer
  !-------------------
  ! DWD BUFR3 routines
  !-------------------
  use mo_bufr_dwd, only: t_bufr,            &!
                         bufr_open_file,    &!
                         bufr_close_file,   &!
                         bufr_read_bufr,    &!
                         bufr_destroy,      &!
                         bufr_get_sections, &!
                         bufr_get_data       !
  !----------------------------------------------------
  ! required for generic communication routines (bcast)
  !----------------------------------------------------
  use mo_mpi_dace, only: dace                ! MPI group info

  implicit none
  !----------------
  ! Public entities
  !----------------
  private
  public :: read_obs_bufr ! module procedure: set observation from BUFR
! public :: p_bcast       ! generic bcast routines
  !-----------------
  ! Private entities
  !-----------------
  integer ,parameter :: n_temp =   1000 ! elements to allocate in a chunk
  integer ,parameter :: n_obs  = 100000
  !----------------------------------------
  ! Parameters to pass to DWD BUFR routines
  !----------------------------------------
  integer :: ieof          ! output parameter: 0:OK, 1:EOF

!=============================================================================
contains
!=============================================================================
  !-----------------------------------------------
  ! Routine to set  observation data type variable
  !            from BUFR data
  ! new records are appended to present data
  !-----------------------------------------------
  subroutine read_obs_bufr (obs, cc)
  type (t_obs)     ,intent(inout),optional :: obs     ! observations data type
  integer          ,intent(in)   ,optional :: cc      ! part of year ccyy

    !================
    ! local variables
    !================
    character(len=128)   :: file         !
    type (t_spot)        :: spt          !
    type (t_spot) ,save  :: empty        !
    type (t_head)        :: head         !
    type (t_bufr) ,save  :: bufr         ! BUFR record
    type (t_time)        :: time         ! observation time
    type (t_time)        :: db_time      ! data bank time
    type (t_use)         :: use          ! status varuable
    integer              :: bufr_type    ! BUFR message type    read
    integer              :: bufr_subtype ! BUFR message subtype read
    integer              :: centre       ! generating centre
    integer              :: ccyy         ! year
!   integer              :: degf         ! degrees of freedom
    logical              :: lkeep        ! true to keep current observation
    integer              :: nkeep        ! number of observations to keep
    integer              :: qcf          ! quality control flag
    integer              :: ifile        ! file counter
    integer              :: obstype      ! observaton type
    integer              :: report_subt  ! observaton subtype (Datenbankkennz.)
    integer              :: report_subti ! observaton subtype index
    type (t_obsid)       :: obsid        ! observation id table entry
    integer              :: i_source     ! position in source file (record)
    integer              :: entry1,entry ! position in source file (subset)
    integer              :: pass         ! 1:read inventory, 2:read data
    type(t_bufr_inv),pointer :: binv
    !-------------------------------------
    ! statistics on bufr records processed
    !-------------------------------------
    integer           :: irecs           ! number of records accepted
    integer           :: iobst   (n_ot)  ! number of reports accepted
    integer           :: jobst   (n_ot)  ! number of reports read
    integer           :: jother          ! number of other reports read
    integer           :: isubset (n_ot)  ! number of subsets read
    integer           :: nsubset         ! number of subsets in record
    integer           :: tsubset
    integer           :: types (0:256,0:256)
    integer           :: dbkz  (-1:10000)
    character(len=16) :: mnemonic
    character(len=80) :: comment
    integer           :: i,j
    !------------------------------------
    ! set pass, only run on dedicated PEs
    !------------------------------------
    pass   = 1; if (present(obs)) pass = 2
    if (pass == 1 .and. dace% pe /= p_readbufr) return
    !-------------
    ! set counters
    !-------------
    irecs  = 0
    iobst  = 0
    jobst  = 0
    jother = 0
    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write (6,*) repeat('-',79)
      write (6,*)
      write (6,*) 'mo_obs_bufr_dwd: reading BUFR files'
    endif
    do ifile = 1, n_source
      !--------------------------
      ! scan list of source files
      !--------------------------
      if (     source(ifile)% filetype /= 1    ) cycle  ! is BUFR file ?
      if (     source(ifile)% obstype  /= 0    ) cycle  ! all types    ?
      if (.not.source(ifile)% used .and.pass==2) goto 97! read file from thisPE
      !---------------------
      ! update record offset
      !---------------------
      if (pass==1) then
        if(ifile>1) bufr_inv% subseto(ifile) = bufr_inv% subseto(ifile-1)
        entry    = 0
      else
        isubset = 0
        if(ifile>1) isubset = bufr_inv% subseto(ifile-1)
        entry   = sum (source(1:ifile-1)% entries)
      endif
      i_source = 0
      !---------------
      ! open bufr file
      !---------------
      file =            source(ifile)% file
      file = path_file (source(ifile)% path,  file)
      call bufr_open_file (bufr, file, output=0)
      if (dace% lpio) then
        write (6,*)
        write (6,*) 'opened : ',trim(file)
      endif
      !----------------
      ! zero statistics
      !----------------
      types = 0
      dbkz  = 0
      !-----------------------
      ! loop over BUFR records
      !-----------------------
      do
        !------------------
        ! preset some flags
        !------------------
!       degf  = 0      ! number of degrees of freedom of observation
        qcf   = 0      ! quality control flag
        lkeep = .true. ! keep this observation for assimilation
        !---------------
        ! read bufr file
        !---------------
        call bufr_read_bufr (bufr, ieof=ieof)
        select case (ieof)
        case (0) ! OK
        case (1) ! EOF
          exit
        case default
          write(0,*) 'read_obs_bufr: read_bufr returns ieof=',ieof
          call finish('read_obs_bufr','read_bufr returns ieof/=0')
        end select
        !----------------------
        ! decode sections 0,1,2
        !----------------------
        call bufr_get_sections (bufr)
        nsubset  = max(1,bufr% sec3% num_subsets)
        i_source = i_source + 1
        irecs    = irecs    + 1
        entry1   = entry    + 1
        entry    = entry    + nsubset
        if (bufr_verb > 2) then
          write (6,*)
          write (6,'(i8,a)') irecs,'. BUFR record'
        endif
        !---------------------------------------
        ! handle empty section 2 (no DWD bufr ?)
        !---------------------------------------
        if (bufr% sec2% len == 0) then
            bufr% sec2% bank = -1
        endif
        !--------------------------------
        ! update statistics on BUFRs read
        !--------------------------------
        types   (bufr% sec1% message_type, bufr% sec1% message_subtype) = &
          types (bufr% sec1% message_type, bufr% sec1% message_subtype) + 1
        dbkz (bufr% sec2% bank) = dbkz (bufr% sec2% bank) + 1
        !=======================================
        ! derive observation type specifications
        !=======================================
        report_subt  = bufr% sec2% bank
        bufr_type    = bufr% sec1% message_type
        bufr_subtype = bufr% sec1% message_subtype
        centre       = bufr% sec1% center
        if (report_subt >= 0) then
          !---------------------------------
          ! derive information from DWD dbkz
          !---------------------------------
          obsid      = obstype_dbkz (report_subt)
        else
          !---------------------------------------------------------------
          ! or from bufr_type, bufr_subtype specified by generating center
          !---------------------------------------------------------------
          obsid      = obstype_bufr (bufr_type, bufr_subtype, centre)
          !------------------------
          ! optionally set DWD dbkz
          !------------------------
          if (derive_dbkz) report_subt   = obsid% dbkz
        endif
        !-------------------------------------------------------
        ! set CMA obstype, BUFR type, BUFR subtype, if not given
        !-------------------------------------------------------
        obstype                          = obsid% obstype
        if (bufr_type   <0) bufr_type    = obsid% bufrtype
        if (bufr_subtype<1 .and. obsid% centre == WMO0_ECMWF) &
                            bufr_subtype = obsid% subtype
        if (obstype < 0) goto 99
        report_subti  = idb_dbk (report_subt, obstype)
        !----------------------------------------------
        ! still some observation type specific counting
        !----------------------------------------------
        select case (obstype)
        case (OT_GPSRO, OT_SYNOP, OT_DRIBU, OT_PAOB, OT_SCATT)
          jobst(obstype) = jobst(obstype) + 1
        case (OT_SATOB)
          select case (bufr% sec2% bank)
          case (DK_SATOB_2, DK_SATOB_3, &
                DK_AMV_EUMETSAT, DK_AMV_GOES, DK_AMV_MODIS, &
                DK_AMV_FY_X, DK_AMV_MTV,                    &
                DK_AMV_NOAA)
!>200705        DK_AMV_MODIS_OLD == DK_AMV_GOES
!         case (DK_SATOB_2, DK_SATOB_3, &
!               DK_AMV_EUMETSAT, DK_AMV_GOES, DK_AMV_MODIS, &
!               DK_AMV_FY_X)
            jobst(obstype) = jobst(obstype) + nsubset
          case default
            jother = jother + 1
          end select
        case (OT_AIREP)
          jobst(obstype) = jobst(obstype) + nsubset
        case (OT_TEMP, OT_PILOT)
          select case (bufr% sec2% bank)
          case (508, 509, 510, 511, &! PILOT      A,B,C,D geopotentielle Hoehe
                512, 513, 514, 515, &! PILOT      A,B,C,D
                516, 517, 518, 519, &! TEMP       A,B,C,D Mobil
                520, 521, 522, 523, &! TEMP       A,B,C,D
                524, 525, 526, 527, &! TEMP       A,B,C,D Mobil
                764, 765, 766, 767, &! PILOT SHIP A,B,C,D geopotentielle Hoehe
                768, 769, 770, 771, &! PILOT SHIP A,B,C,D
                776, 777, 778, 779, &! TEMP  SHIP A,B,C,D
                780, 781, 782, 783)  ! TEMP  DROP A,B,C,D
!         case (512, 513, 514, 515, &! PILOT      A,B,C,D
!               520, 521, 522, 523, &! TEMP       A,B,C,D
!               768, 769, 770, 771, &! PILOT SHIP A,B,C,D
!               776, 777, 778, 779)  ! TEMP  SHIP A,B,C,D
            jobst(obstype) = jobst(obstype) + 1
          case default
            jother = jother + 1
          end select
        end select
        !----------------------------
        ! derive time, data bank time
        !----------------------------
        ccyy = bufr% sec1% year
        if(ccyy< 100 .and. present(cc)) ccyy = ccyy + cc * 100
        if(ccyy==100) ccyy = ccyy + 1900
        if(ccyy< 100) ccyy = ccyy + 2000
        call init_time (time, yyyy = ccyy,              &
                                mo = bufr% sec1% month, &
                                dd = bufr% sec1% day,   &
                                hh = bufr% sec1% hour,  &
                                mi = bufr% sec1% minute )
        ccyy = bufr% sec2% year_dec
        if(ccyy< 100 .and. present(cc)) ccyy = ccyy + cc * 100
        if(ccyy==100) ccyy = ccyy + 1900
        if(ccyy< 100) ccyy = ccyy + 2000
        call init_time (db_time, yyyy = ccyy,                  &
                                   mo = bufr% sec2% month_dec, &
                                   dd = bufr% sec2% day_dec,   &
                                   hh = bufr% sec2% hour_dec,  &
                                   mi = bufr% sec2% min_dec    )
        head% obstype     = obstype
        head% dbkz        = report_subt
        head% modtype     = rept_char(obstype)% mod
        head% buf_type    = bufr_type
        head% buf_subtype = bufr_subtype
        head% codetype    = obsid% codetype
        head% time        = time
        head% db_time     = db_time
        head% idbk        = report_subti
        head% source      = ifile
        head% record      = i_source
        head% id          = entry1
        head% center      = bufr% sec1% center
        head% subcenter   = bufr% sec1% sub_center
        !-----------------
        ! update inventory
        !-----------------
        binv => bufr_inv(obstype)
        if (pass==1) then
          binv% file(head% source) = .true.
          binv% nrec               = binv% nrec           + 1
          binv% nsubset            = binv% nsubset        + nsubset
          binv% subseto(ifile)     = binv% subseto(ifile) + nsubset
          source(ifile)% entries   = entry
          goto 98
        endif
        !----------------------------------------
        ! filter subsets in case of parallel read
        !----------------------------------------
        tsubset          = isubset(obstype) + (nsubset+1)/2
        isubset(obstype) = isubset(obstype) +  nsubset
        if (tsubset <= binv% subsets) then
          goto 98
        endif
        if (tsubset >  binv% subsetl) then
          if (all(isubset >= bufr_inv% subsetl)) then
            call bufr_destroy (bufr)
            exit
          endif
          goto 98
        endif
        !--------------------------------------------
        ! perform simple generic check on report type
        !--------------------------------------------
        call check_report_0 (use, head, nsubset)
        if (use% state <= STAT_DISMISS) then
          lkeep = .false.
          goto 99
        endif
        !------------
        ! decode data
        !------------
        call bufr_get_data (bufr)
        !------------------
        ! create new report
        !------------------
        spt                  = empty
        spt% use             = use
        spt% hd              = head
        !--------------------------------------------
        ! call observation type specific read-routine
        !--------------------------------------------
        select case (obstype)
        case (OT_GPSRO)
          call read_gpsro_bufr (bufr, spt, obs, lkeep)
          if (lkeep) iobst(obstype) = iobst(obstype) + 1
        case (OT_SYNOP, OT_DRIBU, OT_PAOB, OT_SCATT)
          !------------------------------
          ! SYNOP surface data - land,sea
          !------------------------------
          call read_synop_bufr (bufr, spt, obs, lkeep, cc)
          if (lkeep) iobst(obstype) = iobst(obstype) + 1
        case (OT_SATOB)
          !----------------------------
          ! satellite wind observations
          !----------------------------
          call read_amv_bufr (bufr, spt, obs, lkeep, nkeep)
          if (lkeep) iobst(obstype) = iobst(obstype) + nkeep
        case (OT_TEMP, OT_PILOT)
          !------------------------------------------
          ! vertical soundings (other than satellite)
          !------------------------------------------
          call read_temp_bufr (bufr, spt, obs, lkeep, cc)
          if (lkeep) iobst(obstype) = iobst(obstype) + 1
        case (OT_AIREP)
          !-----------------
          ! Aircraft reports
          !-----------------
          call read_airep_bufr (bufr, spt, obs, lkeep, nkeep)
          if (lkeep) iobst(obstype) = iobst(obstype) + nkeep
        case default
          call dism_report_0 (use, head, CHK_NOIMPL) !Buchfuehrung stimmt nicht
          jother = jother + 1                        !++++ lost < 0 ++++
        end select
        !------------
        ! end of loop
        !------------
99      continue
        !---------
        ! printout
        !---------
        if (bufr_verb > 0) then
          write(6,*)irecs,qcf,lkeep
          if (bufr_pause) read(5,*)
        endif
98      continue
        !------------------------------
        ! release memory of BUFR record
        !------------------------------
        call bufr_destroy (bufr)
      end do
      !----------------
      ! close bufr file
      !----------------
      call bufr_close_file (bufr)
      if (dace% lpio) then
        print *
        print *,'read_obs_bufr: ',trim(file)
      endif
      !-----------------
      ! print statistics
      !-----------------
      if (dace% lpio) then
        write (6,*)
        write (6,*) 'BUFR file : ',trim(file)
        write (6,*)
        write (6,'(a)') ' type,subt:records  comment'
        do i=0,256
          do j=0,256
            if (types (i,j) > 0) then
              call wmo_mnem ('A', i, mnemonic, comment)
              write (6,'(1x,i4,a,i4,a,i6," ( ",a," )")') &
              i,',',j,':',types (i,j), trim(mnemonic)
            endif
          end do
        end do
        write (6,*)
        write (6,'(a)') ' DB-Kennz.:records  comment'
        do i=-1,10000
          if (dbkz(i) > 0) then
            call dbkz_mnem (i, mnemonic, comment)
            write (6,'(1x,i9,a,i6," - ",a)') &
              i,':',dbkz(i), trim (comment)
          endif
        end do
        write (6,*)
      endif
97    continue
      !--------------------------------
      ! write pending dismissed reports
      !--------------------------------
      if (pass==2) call write_pending
      !-------------------------------------------
      ! thinning, release unused memory
      ! (temporarily unless input is parallelized)
      !-------------------------------------------
      if (present(obs)) then
        if (obs% n_obs > 0) then
!         call thin (obs)
          call test_obs (obs, 0, 'read_obs_bufr before release_mem',0)
          call release_mem (obs, keep = obs% spot% use% state >= STAT_PAS_REJ)
          call test_obs (obs, 0, 'read_obs_bufr after release_mem',0)
        end if
      end if
      !--------------------------------
      ! write pending dismissed reports
      !--------------------------------
      if (pass==2) call write_pending
    end do
    if (present(obs)) call test_obs (obs, 0, 'read_obs_bufr',0)
    !--------------
    ! write summary
    !--------------
    if (dace% lpio) then
      write(6,'(a,i8,a)') 'read_obs_bufr:',irecs, ' records    read'
      do i = 1, n_ot
        if (jobst(i) > 0) write(6,'(a,i8,1x,a,a)')                &
            'read_obs_bufr:',jobst(i), obstyp(i)%name, ' data read'
        if (iobst(i) > 0) write(6,'(a,i8,1x,a,a)')                    &
            'read_obs_bufr:',iobst(i), obstyp(i)%name, ' data accepted'
      end do
      write(6,'(a,i8,a)') 'read_obs_bufr:',jother,' other data read'
      write(6,'()')
!     call print_rept_stat &
!       (comment='Observation REPORT Statistics after BUFR read')
    endif
  end subroutine read_obs_bufr
!==============================================================================
end module mo_obs_bufr_dwd
