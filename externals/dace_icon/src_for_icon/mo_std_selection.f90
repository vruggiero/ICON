!
!+ organise selection of active reports for ZTD and STD operator
!
Module mo_std_selection
!
! Description:
!   organise selection of active reports for ZTD and STD operator
!
! Current Code Owner: DWD, Michael Bender
!    phone: +49 69 8062 2710
!    fax:   +49 69 8062 3721
!    email: michael.bender@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_47        2016-06-06 Andreas Rhodin
!  new module to organise selection of active reports for ZTD and STD operator
! V1_50        2017-01-09 Andreas Rhodin
!  namelist /gpsgb_select/ for selecting GPSGB observations
! V1_51        2017-02-24 Michael Bender
!  set table of GNSS processing centers and their products
!
! Code Description:
! Language: Fortran 2008.
! Software Standards:
!
!==============================================================================

  !=============
  ! modules used
  !=============
  use kind_parameters, only: wp, i1, i2       ! working precision kind parameter
  use environment,     only: model_abort      ! abort in case of error
#ifndef __COSMO__
  use mo_namelist,     only: position_nml,   &! routine to position nml group
                             nnml,           &! namelist fortran unit number
                             POSITIONED       ! position_nml: OK    return flag
  use mo_mpi_dace,     only: dace,           &! MPI group info
                             p_gather,       &! generic MPI gather
                             p_bcast          ! broadcast routine
#endif
  use mo_t_table,      only: value_name,     &! derive value from name
                             INVALID_VALUE    ! return value from value_name
  use mo_fdbk_tables,  only: ST_ACTIVE,      &! 'active'  status flag value
                             ST_PASSIVE,     &! 'passive' status flag value
               tab_status => status,         &! status variable table
               tab_varno  => varno,          &! variable number table
                             init_fdbk_tables ! initialise the tables
#ifndef __COSMO__
  use mo_dace_string,  only: toupper          ! convert to UPPERCASE
  use mo_algorithms,   only: index            ! sort
#else
  use utilities,       only: to_upper         ! convert to UPPERCASE
#endif

  implicit none

  !================
  ! public entities
  !================
  private
  public :: read_nml_gpsgb_select ! read namelist /GPSGB_SELECT/
  public :: std_status            ! set std report status after input
  public :: name2spec             ! derive processing specifications
  public :: print_tables          ! printout of station and centers used
  public :: t_gnss_center         ! type defining processing center and product
  public :: cntr_pre              ! list of processing centers and products
  public :: list                  ! centers/products selected in namelist
  public :: name2prod             ! convert EGVAP name to center/product ID
  public :: prod2name             ! convert center/product ID to  EGVAP name

  !------------------------------------------------------
  ! list element of namelist group /GPSGB_SELECT/ entries
  !------------------------------------------------------
  integer ,parameter :: ml = 3000  ! max. number of list entries
  integer            :: nl =    0  ! actual number of list entries
  integer ,parameter :: mc =   30  ! max. number of criteria
  integer ,parameter :: ms =  100  ! max. number of stations

  type t_std_select
    !-------------------
    ! selection criteria
    !-------------------
    integer          :: center     (mc)  ! data provider
    integer          :: sub_center (mc)  ! data provide
    integer          :: varno            ! varno: ZTD or STD
    integer          :: product    (mc)  ! provider specific
    character(len=4) :: cnt_prd    (mc)  ! center/product code
    character(len=8) :: station    (ms)  ! station (w/o product code)
    integer          :: ns               ! number of station entries
    !-------
    ! action
    !-------
    integer :: state      ! dismiss, passive, active
    integer :: preference ! for later use
  end type t_std_select

  !----------------------------------------------
  ! list of namelist group /GPSGB_SELECT/ entries
  !----------------------------------------------
  type(t_std_select) ,save :: list (ml)

  !---------------------------------------------------------
  ! table to relate e-gvap name to station/center/processing
  !---------------------------------------------------------
  type t_gnss_station
    integer           :: i
    character(len=10) :: name_center
    character(len=10) :: station
    character(len=4)  :: center4
    integer           :: ilat     ! * 100000
    integer           :: ilon     ! * 100000
    integer           :: n
  end type t_gnss_station

  type t_gnss_center
    integer           :: i
    character(len=4)  :: center4
    integer           :: wmo_center
    integer           :: wmo_subcenter
    character(len=16) :: center
    integer           :: icenter
    integer           :: processing
    integer           :: n
  end type t_gnss_center

  integer :: n_stat = 0
  integer :: n_cntr = 0
  logical :: first  = .true.

  type (t_gnss_station) :: stat (15000)
  type (t_gnss_center)  :: cntr ( 1000)

  !--------------------------
  ! pre-defined table entries (examples, currently not used)
  !--------------------------
  type (t_gnss_station) ,parameter :: stat_pre (1) =     &
       [t_gnss_station (0, 'zzzz-zzzz ','zzzz    ','zzzz',0,0,0)]

  !----------------------------------------------------------------
  ! list of GNSS processing centers and their products, 2015 - 2022
  !----------------------------------------------------------------
  ! GFZ Potsdam, list of products:
  ! GFZ_  - EPOS 6 ZTDs (E-GVAP)
  ! GFS0  - EPOS 6 STDs, ASCII per FTP
  ! GFS1  - EPOS 8 STDs, ASCII per FTP
  ! GF1G  - EPOS 8 ZTDs (E-GVAP), global product, 1.5 h delay
  ! GF1R  - EPOS 8 ZTDs (E-GVAP), rapid product, 1h 15m delay

  ! Experimental, DWD
  ! GF1U  - Ultra rapid => Slants in SkyT, ZTDs per E-GVAP ???

  ! GS1G  - EPOS 8 STDs, global product, DWD AFD => database
  ! GS1R  - EPOS 8 STDs, rapid product, DWD AFD => database
  ! GS1U  - EPOS 8 STDs, ultra rapid product, DWD AFD => database
  !
  ! ZTDs and STDs are never the same product!
  !----------------------------------------------------------------
  type (t_gnss_center) ,parameter :: cntr_pre (75) = [               &
        t_gnss_center (0,'ASI_', 74, 21,'Telespazio      ', 21, 1, 0)&
       ,t_gnss_center (0,'ASIC', 74, 21,'Telespazio      ', 21, 2, 0)&
       ,t_gnss_center (0,'AUT1', 74, 24,'Aristo., Greece ',102, 1, 0)&
       ,t_gnss_center (0,'BKG_', 74, 30,'Amt Kartographie', 30, 1, 0)&
       ,t_gnss_center (0,'BKGA', 74, 30,'Amt Kartographie', 30, 2, 0)&
       ,t_gnss_center (0,'BMEG', 74,  0,'Budapest,Hungary',101, 1, 0)&
       ,t_gnss_center (0,'CONH', 74, 43,'USA             ', 43, 1, 0)&
       ! GFZ Potsdam -------------------------------------------------
       ,t_gnss_center (0,'GFZ_', 74, 23,'GFZ, Potsdam    ', 23, 1, 0)&
       ,t_gnss_center (0,'GFS0', 74, 23,'GFZ, Potsdam    ', 23, 2, 0)&
       ,t_gnss_center (0,'GFS1', 74, 23,'GFZ, Potsdam    ', 23, 3, 0)&
       ,t_gnss_center (0,'GF1R', 74, 23,'GFZ, Potsdam    ', 23, 4, 0)&
       ,t_gnss_center (0,'GF1G', 74, 23,'GFZ, Potsdam    ', 23, 5, 0)&
       ,t_gnss_center (0,'GF1U', 78,  0,'GFZ, Potsdam    ', 23, 6, 0)&
       ,t_gnss_center (0,'GS1R', 78,  0,'GFZ, Potsdam    ', 23, 7, 0)&
       ,t_gnss_center (0,'GS1G', 78,  0,'GFZ, Potsdam    ', 23, 8, 0)&
       ,t_gnss_center (0,'GS1U', 78,  0,'GFZ, Potsdam    ', 23, 9, 0)&
       ! GFZ reprocessing, test
       ! ZTD products - GF*
       ,t_gnss_center (0,'GF01', 78,  0,'GFZ, Potsdam    ', 23,99, 0)&
       ,t_gnss_center (0,'GF02', 78,  0,'GFZ, Potsdam    ', 23,98, 0)&
       ,t_gnss_center (0,'GF03', 78,  0,'GFZ, Potsdam    ', 23,97, 0)&
       ,t_gnss_center (0,'GF04', 78,  0,'GFZ, Potsdam    ', 23,96, 0)&
       ,t_gnss_center (0,'GF05', 78,  0,'GFZ, Potsdam    ', 23,95, 0)&
       ,t_gnss_center (0,'GF06', 78,  0,'GFZ, Potsdam    ', 23,94, 0)&
       ,t_gnss_center (0,'GF07', 78,  0,'GFZ, Potsdam    ', 23,93, 0)&
       ,t_gnss_center (0,'GF08', 78,  0,'GFZ, Potsdam    ', 23,92, 0)&
       ,t_gnss_center (0,'GF09', 78,  0,'GFZ, Potsdam    ', 23,91, 0)&
       ! STD products - GS*
       ,t_gnss_center (0,'GS01', 78,  0,'GFZ, Potsdam    ', 23,90, 0)&
       ,t_gnss_center (0,'GS02', 78,  0,'GFZ, Potsdam    ', 23,89, 0)&
       ,t_gnss_center (0,'GS03', 78,  0,'GFZ, Potsdam    ', 23,88, 0)&
       ,t_gnss_center (0,'GS04', 78,  0,'GFZ, Potsdam    ', 23,87, 0)&
       ,t_gnss_center (0,'GS05', 78,  0,'GFZ, Potsdam    ', 23,86, 0)&
       ,t_gnss_center (0,'GS06', 78,  0,'GFZ, Potsdam    ', 23,85, 0)&
       ,t_gnss_center (0,'GS07', 78,  0,'GFZ, Potsdam    ', 23,84, 0)&
       ,t_gnss_center (0,'GS08', 78,  0,'GFZ, Potsdam    ', 23,83, 0)&
       ,t_gnss_center (0,'GS09', 78,  0,'GFZ, Potsdam    ', 23,82, 0)&
       ! GFZ Potsdam -------------------------------------------------
       ,t_gnss_center (0,'GOPG', 74, 24,'Geodetic Observ.', 24, 1, 0)&
       ,t_gnss_center (0,'GOP1', 74, 24,'Geodetic Observ.', 24, 2, 0)&
       ,t_gnss_center (0,'IGE_', 74, 35,'I. Geografico   ', 35, 1, 0)&
       ,t_gnss_center (0,'IGE2', 74, 35,'I. Geografico   ', 35, 2, 0)&
       ,t_gnss_center (0,'IGER', 74, 35,'I. Geografico   ', 35, 3, 0)&
       ,t_gnss_center (0,'JMA',  34,  0,'JMA, Japan      ',134, 1, 0)&
       ,t_gnss_center (0,'KNM3', 74, 33,'Royal Met. I.   ', 33, 3, 0)&
       ,t_gnss_center (0,'KNM4', 74, 33,'Royal Met. I.   ', 33, 4, 0)&
       ,t_gnss_center (0,'LPT_', 74, 26,'SwissTopo       ', 26, 1, 0)&
       ,t_gnss_center (0,'LPTR', 74, 26,'SwissTopo       ', 26, 2, 0)&
       ,t_gnss_center (0,'LPTX', 74, 26,'SwissTopo       ', 26, 3, 0)&
       ,t_gnss_center (0,'METO', 74,  0,'UK Metoffice    ',  0, 1, 0)&
       ,t_gnss_center (0,'METG', 74,  0,'UK Metoffice    ',  0, 2, 0)&
       ,t_gnss_center (0,'METG', 74,256,'UK Metoffice    ',  0, 2, 0)&
       ,t_gnss_center (0,'MTGH', 74,  0,'UK Metoffice    ',  0, 3, 0)&
       ,t_gnss_center (0,'MTRH', 74,  0,'UK Metoffice    ',  0, 4, 0)&
       ,t_gnss_center (0,'MTRS', 74,  0,'UK Metoffice    ',  0, 5, 0)&
       ,t_gnss_center (0,'NGA1', 74, 34,'Lantmateriet    ', 34, 1, 0)&
       ,t_gnss_center (0,'NGA2', 74, 34,'Lantmateriet    ', 34, 2, 0)&
       ,t_gnss_center (0,'NGII', 40,245,'NGII, Korea     ',140, 1, 0)&
       ,t_gnss_center (0,'NMPT', 40,245,'NGII, Korea     ',140, 2, 0)&
       ,t_gnss_center (0,'NMSC', 40,245,'NGII, Korea     ',140, 3, 0)&
       ,t_gnss_center (0,'NOAA', 74, 40,'NOAA, USA       ', 40, 1, 0)&
       ,t_gnss_center (0,'NOAA', 59,  0,'NOAA, USA       ', 40, 1, 0)&
       ,t_gnss_center (0,'NOAA', 59,256,'NOAA, USA       ', 40, 1, 0)&
       ,t_gnss_center (0,'ROBH', 74, 37,'Royal Observ.   ', 37, 2, 0)&
       ,t_gnss_center (0,'ROBG', 74, 37,'Royal Observ.   ', 37, 1, 0)&
       ,t_gnss_center (0,'ROBQ', 74, 37,'Royal Observ.   ', 37, 3, 0)&
       ,t_gnss_center (0,'ROBT', 74, 37,'Royal Observ.   ', 37,99, 0)&
       ,t_gnss_center (0,'SGN_', 74, 29,'I. Geographique ', 29, 1, 0)&
       ,t_gnss_center (0,'SGN1', 74, 29,'I. Geographique ', 29, 2, 0)&
       ,t_gnss_center (0,'SGNC', 74, 29,'I. Geographique ', 29, 3, 0)&
       ,t_gnss_center (0,'SGN2', 74, 29,'I. Geographique ', 29, 4, 0)&
       ,t_gnss_center (0,'SGN3', 74, 29,'I. Geographique ', 29, 5, 0)&
       ,t_gnss_center (0,'SGN4', 74, 29,'I. Geographique ', 29, 6, 0)&
       ,t_gnss_center (0,'SGNR', 74, 29,'I. Geographique ', 29, 7, 0)&
       ,t_gnss_center (0,'SGNN', 74, 29,'I. Geographique ', 29, 8, 0)&
       ,t_gnss_center (0,'SGNP', 74, 29,'I. Geographique ', 29, 9, 0)&
       ,t_gnss_center (0,'WLIT', 74, 41,'Wroclaw, Poland ', 41, 1, 0)&
       ,t_gnss_center (0,'WUEL', 74, 41,'Wroclaw, Poland ', 41, 2, 0)&
       ,t_gnss_center (0,'WUHN', 74,  0,'Wuhan, China    ',201, 1, 0)]

!=====================================================================
contains
!=====================================================================

  subroutine read_nml_gpsgb_select (unit)
  !-----------------------------
  ! read namelist /GPSGB_SELECT/
  !-----------------------------
    integer, optional, intent(in) :: unit  ! Fortran unit number for namelist
                                           ! for COSMO

    !----------------
    ! local variables
    !----------------
    integer           :: ierr, ipos
    logical ,save     :: first      = .true.
    integer ,save     :: preference = 0
    integer           :: state_
    integer           :: varno_
#ifdef __COSMO__
    integer           :: nnml   ! Fortran unit number for namelist
#endif

    !========================
    ! namelist /GPSGB_SELECT/
    !========================
    !-------------------
    ! selection criteria
    !-------------------
    integer           :: center     (mc) ! data provider
    integer           :: sub_center (mc) ! data provider
    character(len=8)  :: varno           ! varno: ZTD or STD
    integer           :: product    (mc) ! provider specific
    character(len=4)  :: cnt_prd    (mc) ! center/product code
    character(len=8)  :: station    (ms) ! station (w/o product code)
    !-------
    ! action
    !-------
    character(len=8) :: state      ! dismiss, passive, active

    namelist /GPSGB_SELECT/ center, sub_center, varno, product, &
                            state, cnt_prd, station

    !-------------------------------
    ! COSMO: Set Fortran unit number
    !-------------------------------
#ifdef __COSMO__
    if (present(unit)) then
       nnml = unit
    end if
#endif

    !------------------------
    ! read namelist only once
    !------------------------
    if (.not. first) return
    call init_fdbk_tables ()

    do
      !-------------
      ! set defaults
      !-------------
      center     = -1
      sub_center = -1
      varno      = ''
      product    = -1
      state      = ''
      cnt_prd    = ''
      station    = ''
      !--------------
      ! read namelist
      !--------------
#ifndef __COSMO__
      if (dace% lpio) then
        call position_nml ('GPSGB_SELECT' ,lrewind=first ,status=ipos)
        select case (ipos)
        case (POSITIONED)
#if defined(__ibm__)
          !-------------------
          ! catch error on IBM
          !-------------------
          read (nnml ,nml=GPSGB_SELECT ,iostat=ierr)
          if (ierr/=0) call model_abort (-1,-1,                             &
                                 'GPSGB: ERROR in namelist /GPSGB_SELECT/', &
                                 'read_nml_gpsgb_select'           )
#else
          read (nnml ,nml=GPSGB_SELECT)
#endif
        case default
          if (first) then
            write (6,'(a)') repeat('-',79)
            write (6,'()')
            write (6,'(a)')      '  no namelist /GPSGB_SELECT/ present!'
            write (6,'()')
          endif
        end select
      endif
      first = .false.
      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ipos, dace% pio)
      if (ipos /= POSITIONED) exit
      if (dace% lpio) then
#else
         ! COSMO: Assume that namelist file has just been opened
         !        and was positioned at the beginning.

         read (nnml ,nml=GPSGB_SELECT ,iostat=ierr)
         if (ierr < 0) then
            ! end of file
            if (first) then
               write (6,'(a)') repeat('-',79)
               write (6,'()')
               write (6,'(a)')      '  no namelist /GPSGB_SELECT/ present!'
               write (6,'()')
            endif
            exit
         else if (ierr > 0) then
            ! read error
            call model_abort (-1,-1,                                     &
                              'GPSGB: ERROR in namelist /GPSGB_SELECT/', &
                              'read_nml_gpsgb_select'             )
         end if
         first = .false.
#endif

        !----------------------------------------
        ! convert names to integer representation
        !----------------------------------------
        state_ = -1
        state  = toupper (state)
        if (state /= '') then
          state_ = value_name (tab_status, state)
          if (state_ == INVALID_VALUE)                                   &
            call model_abort (-1,-1,                                     &
                              'GPSGB: Invalid value for state: '//state, &
                              'read_nml_gpsgb_select'               )
        endif
        if (state_ == ST_ACTIVE) preference = preference + 1
        varno_ = -1
        varno  = toupper (varno)
        if (varno /= '') then
          varno_ = value_name (tab_varno, varno)
          if (varno_ == INVALID_VALUE)                                   &
            call model_abort (-1,-1,                                     &
                              'GPSGB: Invalid value for varno: '//varno, &
                              'read_nml_gpsgb_select'               )
        endif
        !---------
        ! printout
        !---------
        write (6,'(a)') repeat('-',79)
        write (6,'()')
        write (6,'(a)')      '  namelist /GPSGB_SELECT/ read:'
        write (6,'()')
        write (6,'(a,100i6)')    '    center     = ',pack(center,     center     > -1)
        write (6,'(a,100i6)')    '    sub_center = ',pack(sub_center, sub_center > -1)
        write (6,'(a,100(2x,a))')'    cnt_prd    = ',pack(cnt_prd,    cnt_prd   /= '')
        write (6,'(a,100(2x,a))')'    station    = ',pack(station,    station   /= '')
        if (varno_ > -1) then
          write (6,'(a, a6)')    '    varno      = ',     varno
          write (6,'(a, i6)')    '    varno_     = ',     varno_
        endif
        write (6,'(a,100i6)')    '    product    = ',pack(product,    product    > -1)
        write (6,'(a,   a8)')    '    state      = ',     state
        write (6,'(a,   i6)')    '    preference = ',     preference
        write (6,'()')
#ifndef __COSMO__
      endif
#endif

      !----------
      ! broadcast
      !----------
#ifndef __COSMO__
      call p_bcast (center     ,dace% pio)
      call p_bcast (sub_center ,dace% pio)
      call p_bcast (varno_     ,dace% pio)
      call p_bcast (product    ,dace% pio)
      call p_bcast (state_     ,dace% pio)
      call p_bcast (preference ,dace% pio)
      call p_bcast (cnt_prd    ,dace% pio)
      call p_bcast (station    ,dace% pio)
#endif

      !----------------
      ! include in list
      !----------------
      nl = nl + 1
      if (nl > ml) call model_abort (-1,-1,                  &
                                     'GPSGB: nl > ml',       &
                                     'read_nml_gpsgb_select' )
      list (nl)% center     = center
      list (nl)% sub_center = sub_center
      list (nl)% varno      = varno_
      list (nl)% product    = product
      list (nl)% state      = state_
      list (nl)% preference = preference
      list (nl)% cnt_prd    = cnt_prd
      list (nl)% station    = station
      list (nl)% ns         = count (station /= '')
    end do

  end subroutine read_nml_gpsgb_select

!=====================================================================

  subroutine std_status (center, sub_center, varno, product, cnt_prd,&
                         station, status, prefer)
  !----------------------------------------------
  ! set std report status after observation input
  !----------------------------------------------
  integer         ,intent(in)  :: center     ! data provider
  integer         ,intent(in)  :: sub_center ! data provider
  integer         ,intent(in)  :: varno      ! varno: ZTD or STD
  integer         ,intent(in)  :: product    ! provider specific
  character(len=*),intent(in)  :: cnt_prd    ! center/product code
  character(len=*),intent(in)  :: station    ! station (w/o product code)
  integer         ,intent(out) :: status     ! dismiss, passive, active
  integer(i1)     ,intent(out) :: prefer     ! preference for thinning

    integer :: i, ns

    status = -1

    do i = 1, nl
      if (       list(i)% center    (1) > -1 .and.    &
            all (list(i)% center        /= center    )) cycle
      if (       list(i)% sub_center(1) > -1 .and.    &
            all (list(i)% sub_center    /= sub_center )) cycle
      if (       list(i)% varno         > -1 .and.    &
                 list(i)% varno         /= varno      ) cycle
      if (       list(i)% product   (1) > -1 .and.    &
            all (list(i)% product       /= product   )) cycle
      if (       list(i)% cnt_prd   (1) /= ''.and.    &
            all (list(i)% cnt_prd       /= cnt_prd   )) cycle
      if (       list(i)% station   (1) /= ''         ) then
        ns =     list(i)% ns
        if (all (list(i)% station (:ns) /= station   )) cycle
      endif
      if (       list(i)% state         > -1)         &
        status = list(i)% state
        prefer = list(i)% preference
    end do

  end subroutine std_status

!=====================================================================
  !--------------------
  !subroutine name2spec
  !--------------------
  ! Input from E-GVAP ZTDs, BUFR sequence 307022
  ! name_center: NNNN-PPPP - NNNN 4 character station name
  !                          PPPP 4 charcter center/product identifier
  !              The string provides all information about the
  !              station name, the processing center and the
  !              product of the processing center.
  ! cenpro     : Not present
  !              All information is given in "name_center", no
  !              specfic center/product information available in BUFR.
  ! ...
  !
  !
  ! Input from STD BUFR, BUFR sequence 307024
  ! name_center: Station name, up to 20 characters in BUFR
  !              The string is just the station name and does not provide
  !              any further information, e.g. about center/product.
  ! cenpro     : Center/product identifier from BUFR, needs to be present
  !              as this is the only source of product information.
  !              The product provides the same information as the PPPP part
  !              of the E-GVAP "name_center" and can be used in the same way.
  !              6 characters in BUFR.
  ! ...
  subroutine name2spec (name_center, wmo_center, wmo_subcenter,        &
                        cenpro, lat, lon, station, processing, center)
  !----------------------------------------------------
  ! For GNSS observation
  ! derive unique specifications to be used within DACE
  ! from information provided in the BUFR input
  !----------------------------------------------------
  character(len=*) ,intent(in)  :: name_center   ! name/product (EGVAP)
  integer(i2)      ,intent(in)  :: wmo_center    ! center       (WMO)
  integer(i2)      ,intent(in)  :: wmo_subcenter ! subcenter    (WMO)
  character(len=*) ,optional ,intent(in)  :: &
                                   cenpro        ! center/product identifier
  real(wp)         ,intent(in)  :: lat           ! latitude
  real(wp)         ,intent(in)  :: lon           ! lingitude
  character(len=*) ,intent(out) :: station       ! station name  (DACE)
  integer          ,intent(out) :: processing    ! product number(DACE)
  integer          ,intent(out) :: center        ! center  number(DACE)

  integer :: i, n, ilat, ilon, iprod, icntr

  ! Warning
  ! Center/product in STD BUFR is 6 characters
  ! => Update if more than 4 characters are used !!!
  character(len=4) :: prd   ! Center/product identifier


  if (present(cenpro)) then
     ! Use center/product identifier from BUFR
     prd = cenpro
  else
     ! Assume E-GVAP name NNNN-PPPP
     prd = name_center(6:9)
  end if

    !-----------------------------
    ! set up initial table entries
    !-----------------------------
    if (first) then
       first      = .false.
       n          = size (cntr_pre)
       cntr (1:n) = cntr_pre
       n_cntr     = n
       stat(:)    = stat_pre(1)
    endif

    !---------------------------------
    ! set defaults for missing entries
    !---------------------------------
    station    = ''
    processing = -1
    center     = -1

    !-------------------------------------------
    ! check for combined identifier: name+center
    !-------------------------------------------
    if (.not. present(cenpro)) then
       ! Assume E-GVAP name NNNN-PPPP
       ! The E-GVAP product can be derived only if there is a '-'
       ! at the right position. Without the '-' and without 'cenpro'
       ! there is no product information: Return with defaults.
       if (name_center(5:5) /= '-') return
    end if

    !-----------------------------
    ! try: table lookup of station
    !-----------------------------
    ilat = nint (lat * 100000)
    ilon = nint (lon * 100000)
    if (present(cenpro)) then
       ! Use center/product identifier from BUFR
       ! => compare station name and center/product identifier
       do i = 1, n_stat
          if (stat(i)% name_center == name_center .and.  &
              stat(i)% center4 == prd                  ) then
             station    = stat(i)% station
             stat(i)% n = stat(i)% n + 1
             exit
          endif
       end do
    else
       ! Assume E-GVAP name NNNN-PPPP
       ! => "name_center" provides all information about center/product/station
       do i = 1, n_stat
          if (stat(i)% name_center == name_center ) then
             station    = stat(i)% station
             stat(i)% n = stat(i)% n + 1
             exit
          endif
       end do
    end if

    !---------------------------------
    ! failed: generate new table entry
    !---------------------------------
    if (station == '' .and. n_stat < size(stat)) then
      n_stat = n_stat + 1
      stat(n_stat)% name_center = name_center
      if (present(cenpro)) then
         ! name_center is station name, nothing else
         stat(n_stat)% station     = name_center
         station                   = name_center
         stat(n_stat)% center4     = prd
      else
         ! Assume E-GVAP name NNNN-PPPP
         stat(n_stat)% station     = name_center (1:4)
         station                   = name_center (1:4)
         stat(n_stat)% center4     = name_center (6:9)
      end if
      stat(n_stat)% ilat        = ilat
      stat(n_stat)% ilon        = ilon
      stat(n_stat)% n           = 1
    else if (station == '') then
       call model_abort (-1,-1,                                             &
                         "GPSGB, Overflow: Array stat in mo_std_selection", &
                         'name2spec'                                )
    endif

    !---------------------------------------
    ! try: table lookup of center/processing
    !---------------------------------------
    !iprod = 501  ! next free product ID
    !icntr = 501  ! next free center  ID
    iprod = maxval(cntr(1:n_cntr)% icenter)
    if (iprod < 500) then
       iprod = 501        ! next free product ID
    else
       iprod = iprod + 1  ! next free product ID
    end if
    icntr = iprod         ! next free center  ID

    do i = 1, n_cntr
      icntr = max (icntr, cntr(i)% icenter + 1)

      ! Product handling modified at the end of 2022:
      ! Prior to this date 3 criteria have been used to identify
      ! the center and the product: wmo_center, wmo_subcenter and center4.
      ! This check was redundant and has been simplified. For some E-GVAP
      ! products the wmo_subcenter was changed in BUFR and the products
      ! could not be found. None of theses products have been used
      ! operationally and no conflict with old operational settings is
      ! to be expected.
      ! if   (cntr(i)% wmo_center    == wmo_center   &
      ! .and. cntr(i)% wmo_subcenter == wmo_subcenter) then
      if (cntr(i)% center4 == prd) then
         center = cntr(i)% icenter
         iprod  = max (iprod, cntr(i)% processing + 1)
         processing = cntr(i)% processing
         cntr(i)% n = cntr(i)% n + 1
         exit
      endif
    end do
    !---------------------------------
    ! failed: generate new table entry
    !---------------------------------
    if (processing == -1 .and. n_cntr < size(cntr)) then
      processing                  = iprod
      if (center == -1) center    = icntr
      n_cntr = n_cntr + 1
      !cntr(n_cntr)% center4       = name_center(6:9)
      cntr(n_cntr)% center4       = prd
      cntr(n_cntr)% wmo_center    = wmo_center
      cntr(n_cntr)% wmo_subcenter = wmo_subcenter
      cntr(n_cntr)% icenter       = center
      cntr(n_cntr)% n             = 1
      cntr(n_cntr)% processing    = processing
      cntr(n_cntr)% center        = ''
    else if (processing == -1) then
       call model_abort (-1,-1,                                             &
                         "GPSGB, Overflow: Array cntr in mo_std_selection", &
                                'name2spec'                                )
    endif

  end subroutine name2spec

!---------------------------------------------------------------------

  subroutine name2prod (center4, icenter, processing)
  !----------------------------------------------------
  ! For GNSS observation
  ! convert product string identifier (t_gnss_center - center4) to
  ! to processing center ID (t_gnss_center - icenter) and to
  ! product ID (t_gnss_center - processing).
  !----------------------------------------------------
  character(len=*) ,intent(in)  :: center4    ! product (EGVAP)
  integer          ,intent(out) :: icenter    ! center  (DWD)
  integer          ,intent(out) :: processing ! product (DWD)

  logical :: first  = .true.
  character (len=5*size(cntr_pre)), save :: centers
  integer :: i, idx

  !-------------------------------------------------
  ! Create string with all product names in cntr_pre
  !-------------------------------------------------
  if (first) then
     first = .false.
     do i=1, size(cntr_pre)
        idx = 1 + (i-1)*5
        centers(idx:) = cntr_pre(i)%center4
     end do
  end if

  !-------------
  ! Find product
  !-------------
  idx = index(centers, center4)
  if (idx < 1) then
     ! unknown product
     icenter    = -1
     processing = -1
  else
     i = 1 + (idx-1)/5
     icenter    = cntr_pre(i)%icenter
     processing = cntr_pre(i)%processing
  end if

  end subroutine name2prod

!---------------------------------------------------------------------

  subroutine prod2name (icenter, processing, center4)
  !-----------------------------------------------------------
  ! For GNSS observation
  ! convert processing center ID (t_gnss_center - icenter) and
  ! product ID (t_gnss_center - processing) to
  ! product string identifier (t_gnss_center - center4).
  !-----------------------------------------------------------
  integer          ,intent(in)  :: icenter    ! center  (DWD)
  integer          ,intent(in)  :: processing ! product (DWD)
  character(len=*) ,intent(out) :: center4    ! product (EGVAP)

  logical :: first  = .true.
  integer, dimension(size(cntr_pre)), save :: IDlist
  integer :: i, idx, id

  !---------------------------------------------------------
  ! Create integer array with combined center and product ID
  !---------------------------------------------------------
  if (first) then
     first = .false.
     do i=1, size(cntr_pre)
        IDlist(i) = 1000*cntr_pre(i)%icenter + cntr_pre(i)%processing
     end do
  end if

  !------------------------
  ! Find center and product
  !------------------------
  id = 1000*icenter + processing
  idx = -1
  do i=1, size(cntr_pre)
     if (id == IDlist(i)) then
        idx = i
        exit
     end if
  end do
  if (idx > 0) then
     ! center and product found
     center4 = cntr_pre(idx)%center4
  else
     ! unknown center and product
     center4 = ''
  end if

  end subroutine prod2name

!---------------------------------------------------------------------

  subroutine print_tables
    integer :: i, j

    integer :: verbosity

    verbosity = 1

!#if !defined (__COSMO__)
!    integer :: read_pe              ! PE where tables were initialized
!    integer :: stat_size            ! temporary sizes
!    integer :: cntr_size
!    logical :: need_tbl(dace% npe)
!    !--------------------------------------------------
!    ! Ensure that tables are broadcast.  If not, do so.
!    !--------------------------------------------------
!    call p_gather (first, need_tbl, dace% pio)
!    read_pe = -1
!    do i = 1, size (need_tbl)
!       if (.not. need_tbl(i)) then
!          read_pe = i - 1
!          exit
!       end if
!    end do
!    call p_bcast (read_pe, dace% pio)
!    if (read_pe < 0) then
!!      call model_abort (-1,-1, 'tables not initialized',       &
!!                               'mo_std_selection::print_tables')
!       return
!    end if
!    call p_bcast (n_stat, read_pe)
!    call p_bcast (n_cntr, read_pe)
!    stat_size = size (transfer (stat, ["*"]))
!    cntr_size = size (transfer (cntr, ["*"]))
!    call p_bcast_derivedtype   (stat, stat_size, read_pe, (dace% comm))
!    call p_bcast_derivedtype   (cntr, cntr_size, read_pe, (dace% comm))
!#endif

#ifndef __COSMO__
!   if (dace% lpio) then
    if (.not. first) then

      ! no sorting in COSMO, might be added later if this routine
      ! is really used in COSMO ...
      if (n_stat > 0) then 
         stat(:n_stat)% i = index (stat(:n_stat)% name_center)
         stat(:n_stat)    = stat  (stat(:n_stat)% i)
         stat(:n_stat)% i = index (stat(:n_stat)% station)
         stat(:n_stat)    = stat  (stat(:n_stat)% i)
      end if
#endif
      write(6,'()')
      write(6,'(a)') ' GNSS stations:'
      write(6,'()')

      if (verbosity == 1) then

         write(6,'(a,i5)') 'Number of GNSS stations            : ', n_stat
         write(6,'(a,i5)') 'Number of predefined GNSS products : ', n_cntr

      else if (verbosity > 1) then

         write(6,'(a)') '          station    product       lat       lon   n'
         write(6,'()')
         j = 1
         do i = 1, n_stat
            if (i>1) then
               if (stat(i)% station /= stat(i-1)% station) then
                  write(6,'()')
                  j = j + 1
               endif
            endif
            write(6,'(2i6,1x,a,1x,a,2i10,i4)')                     &
                     j,i,stat(i)% station,stat(i)% center4         ,&
                     stat(i)% ilat,   stat(i)% ilon,  stat(i)% n
         end do

      end if  ! if (verbosity == 1) then

      write(6,'()')
      write(6,'(a)') ' GNSS processing centers:'
      write(6,'()')
      write(6,'(a)') ' WMO center/subcenter center-ID product product-ID entries center full-name'
      write(6,'()')

#ifndef __COSMO__
      if (n_cntr > 0) then
         cntr(:n_cntr)% i = index (cntr(:n_cntr)% processing)
         cntr(:n_cntr)    = cntr  (cntr(:n_cntr)% i)
         cntr(:n_cntr)% i = index (cntr(:n_cntr)% wmo_subcenter)
         cntr(:n_cntr)    = cntr  (cntr(:n_cntr)% i)
         cntr(:n_cntr)% i = index (cntr(:n_cntr)% icenter)
         cntr(:n_cntr)    = cntr  (cntr(:n_cntr)% i)
         cntr(:n_cntr)% i = index (cntr(:n_cntr)% wmo_center)
         cntr(:n_cntr)    = cntr  (cntr(:n_cntr)% i)
      end if
#endif
      j = 0
      do i = 1, n_cntr
        if (j /= cntr(i)% icenter) then
          j = cntr(i)% icenter
          write(6,*)
        endif
        write(6,'(i11,i10,i10,1x,a4,i11,i8,4x,a)')     &
           cntr(i)% wmo_center, cntr(i)% wmo_subcenter &
          ,cntr(i)% icenter,    cntr(i)% center4       &
          ,cntr(i)% processing, cntr(i)% n             &
          ,cntr(i)% center
      end do
      write(6,*)

#ifndef __COSMO__
    endif
#endif

  end subroutine print_tables

#ifdef __COSMO__
  ! Wrapper for COSMO subroutine to_upper
  FUNCTION toupper (lower)
    CHARACTER(LEN=*)              ,INTENT(in) :: lower
    CHARACTER(LEN=LEN_TRIM(lower))            :: toupper

    toupper = lower
    CALL to_upper(toupper)

  END FUNCTION toupper
#endif
!=====================================================================
end module mo_std_selection
!=====================================================================
