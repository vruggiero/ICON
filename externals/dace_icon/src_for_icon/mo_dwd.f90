!
!+ Optimum interpolation correlation error parameters
!
MODULE mo_dwd
!
! Description:
!   Provides optimum (OI) interpolation correlation error parameters
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
! V1_8         2009/12/09 Andreas Rhodin
!  changed error handling for subroutine get_inventory
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  adapt to new module mo_grib_invt
! V1_15        2011/12/06 Harald Anlauf
!  formatted printing of some diagnostic output
! V1_22        2013-02-13 Harald Anlauf
!  bugfix to properly read land-sea-mask
! V1_26        2013/06/27 Harald Anlauf
!  Enable reading of analysis error from GRIB2
! V1_43        2015-08-19 Andreas Rhodin
!  remove dependency from external coefficient file 'rszcoef_f'
! V1_50        2017-01-09 Harald Anlauf
!  mo_dwd: remove dependency on coefficient file btfcoef_f by inlining the data
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! A. Rhodin MPIfM/DWD  2000-2005  original source
!
! References: Research Manual 5 (hereafter RM5)
!             ECMWF Data Assimilation
!             Program Documentation
!             Edt: J.Haseler
!             7/88 Revision 1
!
!             Research Manual 1 (hereafter RM1)
!             ECMWF Data Assimilation
!             Scientific documentation
!             3/92 3rd Edition
!==============================================================================
  !-------------
  ! Modules used
  !-------------
  use mo_kind,          only: wp, dp            ! kind parameters
  use mo_mpi_dace,      only: dace              ! MPI group info
  use mo_exception,     only: finish            ! error abort routine
  use mo_fortran_units, only: get_unit_number, &! reserve a FORTRAN unit number
                              return_unit_number! release a FORTRAN unit number
  use mo_emos_grib1,    only: t_grib1           !
  use mo_grib_invt,     only: ver_time          ! verification time from GRIB
  use mo_grib_handling, only: t_inventory,     &!
                              open_gribfile,   &!
                              close_gribfile,  &!
                              get_inventory,   &!
                              print_inventory   !
  use mo_grib12,        only: read_gribex       ! read GRIB message and decode
  use mo_atm_grid,      only: t_grid, print, destruct
  use mo_grib,          only: read, read_multi, BCAST
  use mo_atm_state,     only: t_atm, construct, destruct, allocate, print
  use mo_time,          only: t_time, p_bcast
  use mo_run_params,    only: oldanerr_file,   &! name of analysis error file
                              data,            &! path to constant data
                              path_file         ! concatenate path/filename
  implicit none
  !----------------
  ! Public entities
  !----------------
  private
  public :: t_ccf              ! data type: climatological error coefficients
  public :: t_vsfc             ! data type: climatological 6h forecasterror cf.
  public :: nclerr             ! fitted climatological error coefficients
  public :: nfcoun             ! fitted climatological 6h forecast error cf.
! public :: bt                 ! barotropic height error coefficients
  public :: rs                 ! radiosonde height observation error coeffs.
  public :: ba_err             ! background analysis error
  public :: lba_err            ! background analysis error is present
  public :: emodsm             ! background error scaling factor
  public :: read_OI_coef       ! initialisation routine for this module
  public :: get_analysis_error ! read get_analysis_error
  public :: clean_OI_coef      ! clean up this module
!==============================================================================
  !=================
  ! Type Definitions
  !=================
!------------------------------------------------------------------------------
  !=============================================
  ! Content of correlation model parameter files
  !=============================================
  !--------------------------------------------
  ! Commpressed Climate File, cf. RM8 page 7.44
  !--------------------------------------------
  type t_ccf
    !------------------------
    ! Data Description Record
    !------------------------
    integer          :: nfddln    ! length of ddr                 (240       )
    integer          :: nftims    ! number of data times          ( 12 months)
    integer          :: nfrows    ! number of latitude rows       ( 30)
    integer          :: iflon (30)! no. of longitude points/row   (  1)
    integer          :: nfvar     ! no. of variables              ( 12)
    integer          :: nfpars(12)! number of parameters/variable (var dep)
    integer          :: icdate    ! creation date                 (yymmdd)
    integer          :: ictime    ! creation time                 (hhmmss)
    real(dp)         :: fnlat     ! northern latitude  boundary   ( 87.)
    real(dp)         :: fslat     ! southern latitude  boundary   (-87.)
    real(dp)         :: fwlon     ! western  longitude boundary   (  0.)
    real(dp)         :: felon     ! eastern  longitude boundary   (  0.)
    character(len=80):: ztitle    ! file title                    (file dep)
    character(len=80):: zch   (12)! description of each variable  (var dep)
    !------------
    ! data record
    !------------
    integer     :: im       (12) ! month
    real(dp)    :: cz  (7,30,12) ! height error,6th order polynominal in ln(p)
    real(dp)    :: cx  (7,30,12) ! height vertical transformation coeffs.
    real(dp)    :: czc   (30,12) ! height vertical correlation coefficients
    real(dp)    :: cu  (7,30,12) ! wind error coefficients
    real(dp)    :: cy  (7,30,12) ! wind vertical transformation coeffcients
    real(dp)    :: cuc   (30,12) ! wind vertical correlation coefficients
    real(dp)    :: mu    (30,12) ! geostrophic factor
    real(dp)    :: b     (30,12) ! horizontal correlation length
    real(dp)    :: cr  (7,30,12) ! relative humidity error coefficients
    real(dp)    :: ca  (7,30,12) ! relative humidity transformation coeffs.
    real(dp)    :: crc   (30,12) ! relative humidity vertical correlation coef.
    real(dp)    :: brh   (30,12) ! horizontal relative humidity corr. length
  end type t_ccf
  !------------------------------------
  ! vertical structure file coefficient
  !------------------------------------
  type t_vsfc
    !------------------------
    ! Data Description Record
    !------------------------
    integer     :: nfddln    ! length of ddr                 (240)
    integer     :: nfrows    ! number of latitude rows       (  1)
    integer     :: iflon (1) ! no. of longitude points/row   (  1)
    integer     :: nfvar     ! no. of variables              (  3)
    integer     :: nfpars(3) ! number of parameters/variable (var dep)
    !------------
    ! data record
    !------------
    real(dp)    :: cz    (7) ! height error coefficients
    real(dp)    :: cx    (7) ! height error corr transformation coefficients
    real(dp)    :: czcor     ! height error correlation coefficients
  end type t_vsfc
!==============================================================================
  !=================
  ! Module variables
  !=================
  !--------------------------------------------
  ! Data type for DWD fitted error coefficients
  ! read from files
  !--------------------------------------------
  type (t_ccf)  ,save,target   :: nclerr ! climatological error coefficients
  type (t_ccf)  ,save,target   :: nfcoun ! climatological 6h forecast error cf.
  type (t_atm)  ,save          :: ba_err ! background error
  type (t_grid) ,save ,pointer :: an_grd  => NULL() ! background error grid
  logical                      :: lba_err = .false. ! analysis error is present
  real(wp)  ,save ,allocatable :: emodsm(:,:)       ! analysis error factor
  !-------------------------------------------------
  ! radiosonde height observation error coefficients
  ! formerly read from 'rszcoef_f'
  ! (still used for geopotential height monitoring)
  !-------------------------------------------------
  type (t_vsfc) ,parameter     :: rs                              &
    = t_vsfc ( 240, 1, [1], 3, [7,7,1],                           &
              [-.414225937309887E+05_dp, 0.272553998151483E+05_dp,&
               -.734365484515589E+04_dp, 0.104056180838218E+04_dp,&
               -.819303229905822E+02_dp, 0.340287631576209E+01_dp,&
               -.582990963220769E-01_dp],                         &
              [-.453618838167086E+04_dp, 0.305997529611301E+04_dp,&
               -.852460504900315E+03_dp, 0.125786055220227E+03_dp,&
               -.103751614837180E+02_dp, 0.453801859023436E+00_dp,&
               -.822891534897197E-02_dp],                         &
               0.963621748665073E+00_dp                           )
  !-------------------------------------
  ! barotropic height error coefficients
  ! formerly read from 'btfcoef_f'
  !-------------------------------------
! type (t_vsfc) ,save          :: bt     ! barotropic height error coefficients
  type (t_vsfc) ,parameter     :: bt                              &
    = t_vsfc ( 240, 1, [1], 3, [7,7,1],                           &
              [-.414225937309887E+05_dp, 0.272553998151483E+05_dp,&
               -.734365484515589E+04_dp, 0.104056180838218E+04_dp,&
               -.819303229905822E+02_dp, 0.340287631576209E+01_dp,&
               -.582990963220769E-01_dp],                         &
              [0.635648597637835E+03_dp, -.287857832111344E+03_dp,&
               0.407110858364174E+02_dp, -.116218995521578E-01_dp,&
               -.505829087614377E+00_dp, 0.446542001395749E-01_dp,&
               -.122181844316789E-02_dp],                         &
               0.789331165913676E+00_dp                           )
!==============================================================================
contains
!==============================================================================
  subroutine get_analysis_error

    type (t_inventory) ,pointer :: ana_inv(:)
    type (t_time)               :: time
    !---------------------------
    ! parameters to get_pds_info
    !---------------------------
    type (t_grib1)   :: grib

    !=============
    ! print header
    !=============
    if(dace% lpio) then
      write(6,'()')
      write(6,'(a)') repeat('=',79)
      write(6,'()')
      write(6,'(a)') '  reading analysis error:'
      write(6,'()')
      write(6,'(a)') '  content of grib file '//trim(oldanerr_file)
      write(6,'()')
    endif

    !--------------------------------------
    ! check for analysis error file to read
    !--------------------------------------
    if (oldanerr_file == '/none/') then
      if(dace% lpio) then
        write(6,'(a)') '  reading NO analysis error'
        write(6,'()')
      endif
      return
    endif

    !=====================
    ! read analysis errors
    !=====================
    lba_err = .true.
    !--------------
    ! get inventory
    !--------------
    nullify (ana_inv)
    call get_inventory (ana_inv, trim(oldanerr_file))
    if (size(ana_inv) == 0)  call finish ('get_analysis_error',&
            'file empty or not present: '//trim(oldanerr_file) )
    call print_inventory (ana_inv)
    !----------
    ! read grid
    !----------
    if(dace% lpio) then
      write(6,'()')
      write(6,'(a)') '  reading grid information:'
      write(6,'()')
    endif
    if(.not.associated(an_grd)) allocate (an_grd)
    call read  (an_grd, oldanerr_file, ana_inv, lsm=.false., geosp=.false.)
!               nproc1=1, nproc2=1, comm=dace% comm)
!   call print (an_grd, verbose=.true.)
    !----------------------------------------
    ! set verification time from first record
    !----------------------------------------
    call construct (ba_err, an_grd)
    if (dace% lpio) then
      call open_gribfile (grib, oldanerr_file, 'r')
      call read_gribex   (grib,'J')
      time = ver_time    (grib)
!     time = ana_inv(1)% ti% ver_time
    endif
    call p_bcast (time, dace% pio)
    ba_err% time = time
    !------------
    ! read fields
    !------------
    if(dace% lpio) then
      write(6,'()')
      write(6,'(a)') '  reading analysis errors:'
      write(6,'()')
    endif
    call allocate   (ba_err,'u')
    call read_m     ('u',ba_err% u)
    call allocate   (ba_err,'v')
    call read_m     ('v',ba_err% v)
    call allocate   (ba_err,'geof')
    call read_m     ('geof',ba_err% geof)
    call print (ba_err, verbose=.true.)
    !-----------
    ! close file
    !-----------
    if (dace% lpio) call close_gribfile (grib)
    deallocate (ana_inv)
    !---------------------------------------------------------
    ! allocate analysis error correction factor
    ! preset with values of 1.
    ! the final factor is calculated in subroutine init_fg_cov
    !                                  (module mo_fg_cov)
    !---------------------------------------------------------
    if (.not.allocated(emodsm)) allocate (emodsm(an_grd% nx, an_grd% ny))
    emodsm (:,:) = 1._wp

  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine read_m (name, x)
    character(len=*) ,intent(in)            :: name
    real(wp)         ,intent(out)           :: x (:,:,:,:)

      integer          :: code, table
      integer          :: i

      code = 255; table = 255
      call read_multi (x, ana_inv, ba_err% grid, grib,     &
                        name=name, code=code, table=table, &
                        pio=dace% pio, pmode=BCAST         )

!print *,'read_m: ',name, code,table

      do i=1,size(ba_err% m)
        if (ba_err% m(i)% i% name == name) then
            ba_err% m(i)% i% code  = code
            ba_err% m(i)% i% table = table
!print *,'read_m: set code:',name,code,table
            exit
        endif
      end do

    end subroutine read_m
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine get_analysis_error
!==============================================================================
  subroutine clean_OI_coef
    if (allocated(emodsm)) &
      deallocate (emodsm)
    if (associated(an_grd)) then
      call destruct (an_grd)
      deallocate    (an_grd)
    endif
    call destruct (ba_err)
  end subroutine clean_OI_coef
!==============================================================================
  subroutine read_OI_coef
  !----------------------------------------
  ! Read DWD OI Covariance model parameters
  !----------------------------------------
    !----------------
    ! local variables
    !----------------
    logical :: first = .true.
    integer :: ccf_size !, vsfc_size
    !------------------------------------------------
    ! interfaces for bcast routines
    ! (external routines, defined in mo_mpi_dace.f90)
    !------------------------------------------------
!    interface p_bcast
!      subroutine p_bcast_derivedtype (buffer, count, source, comm)
!      use mo_dwd, only: t_ccf
!      type(t_ccf)  ,intent(inout)        :: buffer    ! variable to bcast
!      integer      ,intent(in)           :: count     ! len(byte) of variable
!      integer      ,intent(in)           :: source    ! source processor index
!      integer      ,intent(in) ,optional :: comm      ! communicator
!      end subroutine p_bcast_derivedtype
!    end interface p_bcast

!    interface p_bcast
!      subroutine p_bcast_derivedtype2 (buffer, count, source, comm)
!      use mo_dwd, only: t_vsfc
!      type(t_vsfc) ,intent(inout)        :: buffer    ! variable to bcast
!      integer      ,intent(in)           :: count     ! len(byte) of variable
!      integer      ,intent(in)           :: source    ! source processor index
!      integer      ,intent(in) ,optional :: comm      ! communicator
!      end subroutine p_bcast_derivedtype2
!    end interface p_bcast
    !--------------------
    ! read data only once
    !--------------------
    if (first) then
      first = .false.
      !----------------------------
      ! read files on I/O processor
      !----------------------------
      if (dace% lpio) then
        print *, 'reading DWD OI correlation model on PE',dace% pe
        call read_ccf   (nclerr, path_file (data,'clerr_f'))
        call print_ccf  (nclerr)
        call read_ccf   (nfcoun, path_file (data,'fgerr_f'))
        call print_ccf  (nfcoun)
!       call read_vsfc  (bt,     path_file (data,'btfcoef_f'))
!       call print_vsfc (bt,     'barotropic height error coefficients')
!       call print_vsfc (rs,'radiosonde height observation error coefficients')
      endif
      !---------------------------------
      ! broadcast data to all processors
      !---------------------------------
      ccf_size  = size (transfer (nclerr ,(/' '/)))
!     vsfc_size = size (transfer (bt     ,(/' '/)))

      call p_bcast_derivedtype  (nclerr, ccf_size, (dace% pio), (dace% comm))
      call p_bcast_derivedtype  (nfcoun, ccf_size, (dace% pio), (dace% comm))
!     call p_bcast_derivedtype2 (bt,    vsfc_size, (dace% pio), (dace% comm))
!     call p_bcast_derivedtype2 (rs,    vsfc_size, (dace% pio), (dace% comm))

    endif
  end subroutine read_OI_coef
!==============================================================================
  subroutine read_ccf (ccf, file)
  !---------------------------------------------------------------
  ! Read Compressed Climate Files (formatted, machine independent)
  !---------------------------------------------------------------
  type(t_ccf)      ,intent(out) :: ccf
  character(len=*) ,intent(in)  :: file
    !----------------
    ! local variables
    !----------------
    integer :: iu ! unit number
    integer :: i  ! loop index (months)
    integer :: l  ! loop index (latitudes)
    !-------------------
    ! explizites format
    !  fuer 64 bit worte
    !-------------------
9010 format (5i20)    ! integer
9020 format (5e21.15) ! real
9030 format (a80)     ! character
    !----------
    ! open file
    !----------
    if (dace% lpio) print *,'read_ccf: reading ',trim (file)
    iu = get_unit_number()
    open (iu ,file=file ,form='formatted' ,status='old' ,action='read')
    !-----------------------------
    ! Read Data Description Record
    !-----------------------------
    read (iu ,9010) ccf% nfddln ,ccf% nftims ,ccf% nfrows ,ccf% iflon &
                   ,ccf% nfvar  ,ccf% nfpars ,ccf% icdate ,ccf% ictime
    read (iu ,9020) ccf% fnlat  ,ccf% fslat  ,ccf% fwlon  ,ccf% felon
    read (iu ,9030) ccf% ztitle ,ccf% zch
    !------------------------------------------------------
    ! Read Data Record (*,30,12) (parameter,latitude,month)
    ! store South to North
    !------------------------------------------------------
    do i=1,12
      READ (iu,9010) ccf% im    (i)
      READ (iu,9020)(ccf% cz(:,l,i) ,ccf% cx(:,l,i) ,ccf% czc (l,i) &
                    ,ccf% cu(:,l,i) ,ccf% cy(:,l,i) ,ccf% cuc (l,i) &
                    ,ccf% mu  (l,i) ,ccf% b   (l,i) ,ccf% cr(:,l,i) &
                    ,ccf% ca(:,l,i) ,ccf% crc (l,i) ,ccf% brh (l,i),l=30,1,-1)
    end do
    !-----------
    ! close file
    !-----------
    close (iu)
    call return_unit_number(iu)
  end subroutine read_ccf
!------------------------------------------------------------------------------
  subroutine print_ccf (ccf)
  !-----------------------------------------
  ! Print Content of Compressed Climate File
  !-----------------------------------------
  type(t_ccf)       ,intent(in) :: ccf   ! dataset to print
   integer :: i, j
   print *,'-----------------------------------------------------------------'
   write(6,'(  a      )')trim(ccf% ztitle)
   write(6,'(a,i10  ,a)')' nfddln=',ccf% nfddln,' length of ddr          (240)'
   write(6,'(a,i10  ,a)')' nftims=',ccf% nftims,' number of data times(12mon.)'
   write(6,'(a,i10  ,a)')' nfrows=',ccf% nfrows,' number of latitude rows (30)'
   write(6,'(a,10x,a)')  ' iflon =',' no. of longitude points/row'
   write(6,'(a,(30i2))') '        ',ccf% iflon
   write(6,'(a,i10  ,a)')' nfvar =',ccf% nfvar ,' no. of variables        (12)'
!  do i=1,size(ccf% nfpars)
!  write(6,'(a,i10,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
!  end do
   write(6,'(a,i10  ,a)')' icdate=',ccf% icdate,' creation date       (yymmdd)'
   write(6,'(a,i10  ,a)')' ictime=',ccf% ictime,' creation time       (hhmmss)'
   write(6,'(a,f10.3,a)')' fnlat =',ccf% fnlat ,' nort. latit. boundary ( 87.)'
   write(6,'(a,f10.3,a)')' fslat =',ccf% fslat ,' sout. latit. boundary (-87.)'
   write(6,'(a,f10.3,a)')' fwlon =',ccf% fwlon ,' west. long.  boundary (  0.)'
   write(6,'(a,f10.3,a)')' felon =',ccf% felon ,' east. long.  boundary (  0.)'
!  do i=1,size(ccf% zch)
!  write(6,'(a,i10,1x,a)')'        ',ccf% im(i),trim(ccf% zch(i))
!  end do

   i=1
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' cz    =',minval(ccf% cz), maxval(ccf% cz)
   i=2
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' cx    =',minval(ccf% cx), maxval(ccf% cx)
   i=3
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' czc   =',minval(ccf% czc), maxval(ccf% czc)
   i=4
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' cu    =',minval(ccf% cu), maxval(ccf% cu)
   i=5
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' cy    =',minval(ccf% cy), maxval(ccf% cy)
   i=6
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' cuc   =',minval(ccf% cuc), maxval(ccf% cuc)
   i=7
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' mu    =',minval(ccf% mu), maxval(ccf% mu)
   i=8
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' b     =',minval(ccf% b), maxval(ccf% b)
   i=9
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' cr    =',minval(ccf% cr), maxval(ccf% cr)
   i=10
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' ca    =',minval(ccf% ca), maxval(ccf% ca)
   i=11
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' crc   =',minval(ccf% crc), maxval(ccf% crc)
   i=12
   write(6,'(i3,1x   ,a)') ccf% im(i),trim(ccf% zch(i))
   write(6,'(a,i10   ,a)')' nfpars=',ccf%nfpars(i),' number of parametrs/variable'
   write(6,'(a,2f15.3,a)')' brh   =',minval(ccf% brh), maxval(ccf% brh)
   print *
   i=1
   do j=1, ccf% nfrows
     write(6,'(i3,3f6.2,f10.0,f6.2,f10.0)') &
       j, ccf% czc(j,1),ccf% cuc(j,1),ccf% mu(j,1),ccf% b(j,1),ccf% crc(j,1),ccf% brh(j,1)
   end do
   print *,'-----------------------------------------------------------------'
  end subroutine print_ccf
!==============================================================================
  subroutine read_vsfc (vsfc, file)
  !-----------------------------------------
  ! read vertical structure coefficient file
  !-----------------------------------------
  type (t_vsfc)    ,intent(out) :: vsfc
  character(len=*) ,intent(in)  :: file
    !----------------
    ! local variables
    !----------------
    integer :: iu   ! unit number
    integer :: i    ! loop index
    integer :: irow ! dummy row index
    !-------------------
    ! explizites format
    !  fuer 64 bit worte
    !-------------------
9010 format (5i20)    ! integer
9020 format (5e21.15) ! real
    !----------
    ! open file
    !----------
    if (dace% lpio) print *,'read_vsfc: reading ',trim (file)
    iu = get_unit_number()
    open (iu ,file=file ,form='formatted' ,status='old' ,action='read')
    !-----------------------------
    ! Read Data Description Record
    !-----------------------------
    read (iu ,9010) vsfc% nfddln ,vsfc% nfrows ,vsfc% iflon &
                   ,vsfc% nfvar  ,vsfc% nfpars
    do i=1,46
      read (iu ,*)
    end do
    !-----------------
    ! Read Data Record
    !-----------------
    READ (iu,9010) irow
    READ (iu,9020) vsfc% cz, vsfc% cx, vsfc% czcor
    !-----------
    ! close file
    !-----------
    close (iu)
    call return_unit_number(iu)
  end subroutine read_vsfc
!------------------------------------------------------------------------------
  subroutine print_vsfc (vsfc, title)
  !-----------------------------------------
  ! read vertical structure coefficient file
  !-----------------------------------------
  type(t_vsfc)     ,intent(in) :: vsfc  ! dataset to print
  character(len=*) ,intent(in) :: title ! title string
   integer :: i,k
   integer ,parameter :: stp(15) = &      ! standard pressure levels to plot
     (/10,20,30,50,70,100,150,200,250,300,400,500,700,850,1000/)
   real(wp) :: p, lnp, e, x
   print *,'-----------------------------------------------------------------'
   print *,title
   write(6,'(a,i10 ,a)')' nfddln=',vsfc% nfddln,' length of ddr          (240)'
   write(6,'(a,i10 ,a)')' nfrows=',vsfc% nfrows,' number of latitude rows ( 1)'
   write(6,'(a,i10 ,a)')' iflon =',vsfc%iflon(1),' no. of longitude points/row'
   write(6,'(a,i10 ,a)')' nfvar =',vsfc% nfvar ,' no. of variables        ( 3)'
   do i=1,size(vsfc% nfpars)
   write(6,'(a,i10,a)')' nfpars=',vsfc%nfpars(i),' number of parametrs/var.'
   end do
   write(6,'(a,(7f10.3))')' cz    =',vsfc% cz
   write(6,'(a,(7f10.3))')' cx    =',vsfc% cx
   write(6,'(a,(7f10.3))')' czcor =',vsfc% czcor
   write(6,*)
   write(6,'(a)')'  k      p      ln(p)     error  x(ln(p))'
   do k=size(stp),1,-1
     p = stp(k) * 100
     lnp = log(p)
     e = 0._wp
     x = 0._wp
     do i = vsfc% nfpars(1), 1, -1
       e = e * lnp + vsfc% cz (i)
     end do
     do i = vsfc% nfpars(2), 1, -1
      x = x * lnp + vsfc% cx (i)
     end do
     write(6,'(i3,f8.0,3f10.3)') k, p, lnp, e, x
   end do
   print *,'-----------------------------------------------------------------'
  end subroutine print_vsfc
!==============================================================================
end module mo_dwd
