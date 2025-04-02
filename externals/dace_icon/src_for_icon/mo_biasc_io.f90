!
!+ Manage IO of bias correction statistics files
!
module mo_biasc_io
!
! Description:
!
!   Manage IO of bias correction statistics files used by RADIANCE and AIRCRAFT
!   bias correction routines.
!   Namelist /BIASCOR/ provides the names of input bias correction statistics
!   files retrieved from the data base.
!   Routines are provided to derive suitable file names for output, and to
!   handle (open and close NetCDF files for reading or writing). The latter
!   are currently used by the AIRCRAFT bias correction only. RADIANCE bias
!   correction still uses a specific ASCII format.
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_13        2011/11/01 Andreas Rhodin
!  Manage IO of bias correction statistics files
! V1_15        2011/12/06 Harald Anlauf
!  Cleanup
! V1_20        2012-06-18 Harald Anlauf
!  mo_biasc_io, mo_conv_bc: verbosity level for bias correction
! V1_22        2013-02-13 Andreas Rhodin
!  ignore unknown observation type in bias correction file name (fake NUMEX)
! V1_27        2013-11-08 Andreas Rhodin
!  add comment lines
! V1_29        2014/04/02 Andreas Rhodin
!  consistently use n_ot instead of n_obstype
! V1_48        2016-10-06 Robin Faulwetter
!  Allow different t_decay values for different channels
! V1_49        2016-10-25 Harald Anlauf
!  don't write zero entries of t_decay
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================

!-------------
! Modules used
!-------------
use mo_kind,         only: wp                ! working precision kind parameter
use mo_mpi_dace,     only: dace,            &! MPI group info
                           p_bcast           ! broadcast routine
use mo_exception,    only: finish,          &! abort routine
                           message           ! write warning
use mo_namelist,     only: position_nml,    &! position namelist
                           nnml,            &! namelist Fortran unit number
                           POSITIONED        ! ok    code from position_nml
use mo_run_params,   only: ana_time,        &! analysis time
                           input,           &! default input path
                           output,          &! default output path
                           nex,             &! experiment number
                           path_file,       &! concatenate: path / file
                           run_time,        &! run-time of 3dvar
                           ana_time          ! analysis time
use mo_time,         only: cyyyymmddhhmm     ! derive string from time variable
use mo_dace_string,  only: char4             ! Conversion: INTEGER->CHAR(LEN=4)
use mo_fdbk_tables,  only: n_ot              ! number of report types defined
use mo_obs_tables,   only: obstyp            ! table of observation types
use mo_rad,          only: m_chan            ! max. number of channels
use mo_fortran_units,only: get_unit_number, &!  reserve a unit number
                           return_unit_number!  release a unit number
use netcdf,          only:                  &! NetCDF f90 interface
                           nf90_open,       &! open   NetCDF file
                           nf90_create,     &! create NetCDF file
                           nf90_close,      &! close NetCDF file
                           nf90_strerror,   &! derive error character string
                           nf90_put_att,    &! write attribute
                           nf90_get_att,    &! read  attribute
                           NF90_CLOBBER,    &! over-write flag for nf90_create
                           NF90_NOWRITE,    &! read       flag for nf90_open
                           NF90_NOERR,      &! NetCDF return code for no error
                           NF90_GLOBAL       ! global attribute flag
implicit none

!================
! Public entities
!================
private
!--------------
! derived types
!--------------
public :: t_bcor_head    ! file header data
!---------
! routines
!---------
public :: bc_filename    ! function: derive file name from observation type etc
public :: read_nml_biasc ! subroutine  :  read namelist /BIASCOR/
public :: new_bc_head    ! subroutine  :  construct new file header
public :: open_bc_write  ! subroutine  :  open  file for writing
public :: open_bc_read   ! subroutine  :  open  file for reading
public :: close_bc       ! subroutine  :  close file
!------------------------------------------
! namelist variables and derived quantities
!------------------------------------------
public :: bc_path        ! path (optional) to shortcut NUMEX/database
public :: bc_files       ! bias correction files to read
public :: bc_paths       ! full pathnames of files to read
public :: bc_obstyp      ! observation type related to each file
public :: nbcf           ! number of files in input list
public :: fallback       ! flag to use template if input file not present
public :: verbose        ! Verbosity level of bias correction

!========================
! derived type definition
!========================

  !---------------------------------------------
  ! len of global attributes (character strings)
  !---------------------------------------------
  integer, parameter ::  TLEN = 28 ! len of title       attribute
  !nteger, parameter ::  HLEN = 80 ! len of history     attribute
  integer, parameter ::  ILEN = 24 ! len of institution attribute
  integer, parameter ::  SLEN = 16 ! len of source      attribute

  !-----------------------------------
  ! Common bias correction file header
  !-----------------------------------
  type t_bcor_head
    character(len=128)  :: path            = ''      ! path/file-name
    character(len=TLEN) :: title           = ''      ! (CF) title
    character(len=ILEN) :: institution     = ''      ! (CF) institution
    character(len=SLEN) :: source          = ''      ! (CF) source
!   character(len=HLEN) &                        ! (CF) history
!              ,pointer :: history(:)      =>NULL()  !
    character(len=5)    :: version         = '00.00' ! (CF) file version number
    integer             :: exp             = -1      ! experiment number
    character(len=12)   :: created         = ''      ! creation date 'yyyymmddhhmm'
    character(len=12)   :: modified        = ''      ! modification date
    character(len=12)   :: first_date      = ''      ! first training date
    character(len=12)   :: last_date       = ''      ! last  training date
    real(wp)            :: t_decay(m_chan) = 0._wp   ! background info decay time
    integer             :: ncid            = -1      ! NetCDF file Id
  end type t_bcor_head

!=================
! Module variables
!=================

  !--------------------
  ! Namelist /BIASCOR/
  !--------------------
  integer ,parameter :: mbcf            = 64  ! max. number of bias corr. files
  character(len= 32) :: bc_files (mbcf) = ''  ! bias correction files to read
  character(len=128) :: bc_path         = ''  ! path to shortcut data base
  logical            :: fallback    = .false. ! use template if no input file
  integer            :: verbose         = 2   ! Verbosity level of bias corr.

  namelist /BIASCOR/ bc_files, bc_path, fallback, verbose

  !-------------------------------------------
  ! quantities derived from namelist variables
  !-------------------------------------------
  integer            :: bc_obstyp(mbcf) = -1  ! observation type related to file
  character(len=256) :: bc_paths (mbcf) = ''  ! full pathname
  integer            :: nbcf            = 0   ! number of files in input list

  !-------------
  ! internal use
  !-------------
  character(len=32)  :: rname                 ! name of active subroutine

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

  function bc_filename (obstype, satid, instrid)
  !---------------------------------------------------------------------------
  ! Derive a suitable file name for writing a bias correction statistics file.
  ! In order to be processed by the data base the files have to follow
  ! the given file name convention derived from observation type and
  ! analysis time (and optionally satellite and instrument id):
  !
  !   bias_TTTT[_SSSSIIII].yyyymmddhhmm
  !
  ! with:
  !       TTTT     = obbservation type       : S2 database specifier
  !       SSSSIIII = satellite/instrument id : S3 specifier (optional)
  !   yyyymmddhhmm = analysis time
  !---------------------------------------------------------------------------
  character(len=32)             :: bc_filename  ! file name to derive
  integer ,intent(in)           :: obstype      ! 3dVar observation type
  integer ,intent(in) ,optional :: satid        ! WMO satellite id
  integer ,intent(in) ,optional :: instrid      ! WMO instrument id

  if   (present(satid)  ) then
    if (present(instrid)) then
      !--------------------------
      ! bias_TTT[TT].yyyymmddhhmm
      !--------------------------
      bc_filename = 'bias_'//trim(obstyp  (obstype)% name) &
                      //'_'//char4        (satid)          &
                           //char4        (instrid)        &
                      //'.'//cyyyymmddhhmm(ana_time)
    else
      !-------------------------------
      ! bias_TTT[TT]_SSSS.yyyymmddhhmm
      !-------------------------------
      bc_filename = 'bias_'//trim(obstyp  (obstype)% name) &
                      //'_'//char4        (satid)          &
                      //'.'//cyyyymmddhhmm(ana_time)
    endif
  else
    !-------------------------------------
    !   bias_TTT[TT]_SSSSIIII.yyyymmddhhmm
    !-------------------------------------
    bc_filename   = 'bias_'//trim(obstyp  (obstype)% name) &
                      //'.'//cyyyymmddhhmm(ana_time)

  endif
  end function bc_filename

!------------------------------------------------------------------------------

  subroutine read_nml_biasc
  !----------------------------
  ! read the namelist /BIASCOR/
  !----------------------------
    integer :: ierr  ! error return parameter
    integer :: i     ! file index
    integer :: j     ! character position index
    integer :: k     ! index
    integer :: iu    ! Fortran unit number
    integer :: ios   ! IOstat return code
    integer :: o     ! position of last '/' in name
    !--------------------------------
    ! set namelist variables defaults
    !--------------------------------
    bc_files  = ''      ! bias correction files to read
    bc_path   = ''      ! path to shortcut data base
    bc_paths  = ''      ! full pathname
    bc_obstyp = -1      ! observation type related to file
    nbcf      =  0      ! number of files in list
    fallback  = .false. ! use template if no input file
    verbose   =  2      ! Verbosity level of bias corr.
    !-------------------------
    ! read namelist /BIASCOR/
    !-------------------------
    if (dace% lpio) then
      call position_nml ('BIASCOR', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=BIASCOR ,iostat=ierr)
        if (ierr/=0) call finish ('read_nml_biasc',&
                                  'ERROR in namelist /BIASCOR/')
#else
        read (nnml ,nml=BIASCOR)
#endif
        !---------------------------------------
        ! derive observation type from file name
        !---------------------------------------
        do i = 1, mbcf
          if (bc_files (i)       == '') cycle
          o = index (bc_files (i), '/', back=.true.)
          if (bc_files (i) (o+1:o+5) /= 'bias_') &
            call finish ('read_nml_biasc',                      &
                         'invalid filename: '//trim(bc_files(i)))
          j          = index (bc_files (i) (o+6:) ,'.')
          k          = index (bc_files (i) (o+6:) ,'_')
          if (k>0) j = min(j,k)
          if (j>1) then
            j = j + 4
            do k = 1, n_ot
              if (bc_files (i) (o+6:o+j) == obstyp(k)% name) bc_obstyp(i) = k
            end do
          endif
          if (bc_obstyp(i) < 1) then
!           call finish  ('read_nml_biasc',                                 &
!                         'invalid obstype in filename: '//trim(bc_files(i)))
            call message ('read_nml_biasc',                                 &
                          'invalid obstype in filename: '//trim(bc_files(i)))
            cycle
          endif
          nbcf        = i
        end do
      end select
      !------------------
      ! printout namelist
      !------------------
      write(6,'(a)') repeat('-',79)
      write(6,'( )')
      write(6,'(a)')  '  Namelist /BIASCOR/ read:'
      write(6,'( )')
      write(6,'(a,l1)')       '    fallback    = ',fallback
      write(6,'(2a)')         '    bc_path     = ',trim (bc_path)
      write(6,'(a,i3)')       '    # of files  = ',nbcf
      write(6,'(a,i3)')       '    verbose     = ',verbose
      iu = get_unit_number ()
      do i = 1, nbcf
        write(6,'(4x,i3,1x,a,"  obstype=",i3)') i, bc_files(i), bc_obstyp(i)
        !------------------------------
        ! check presence of input files
        !------------------------------
        bc_paths(i) = path_file (input, path_file (bc_path, bc_files(i)))
        open  (iu, file=bc_paths(i), action='read', status='old', iostat=ios)
        if (ios/=0) call finish ('read_nml_biasc',                        &
                                 'cannot access file: '//trim(bc_paths(i)))
        close (iu)
      end do
      call return_unit_number (iu)
      write(6,'( )')
    endif
    !-------------------
    ! broadcast namelist
    !-------------------
    call p_bcast (bc_path,   dace% pio)
    call p_bcast (nbcf,      dace% pio)
    call p_bcast (bc_files,  dace% pio)
    call p_bcast (bc_paths,  dace% pio)
    call p_bcast (bc_obstyp, dace% pio)
    call p_bcast (fallback,  dace% pio)
    call p_bcast (verbose,   dace% pio)
  end subroutine read_nml_biasc

!------------------------------------------------------------------------------

  subroutine new_bc_head (head, obstype, satid, instrid, t_decay)
  !-------------------------------------------------
  ! preset a new bias correction file header
  ! currently only used by AIRCRAFT bias correction.
  !-------------------------------------------------
  type(t_bcor_head) ,intent(out)          :: head       ! file header data
  integer           ,intent(in)           :: obstype    ! 3dVar observation type
  integer           ,intent(in) ,optional :: satid      ! WMO satellite id
  integer           ,intent(in) ,optional :: instrid    ! WMO instrument id
  real(wp)          ,intent(in) ,optional :: t_decay(:) ! info decay time

  integer :: n

    !---------------------
    ! set file header data
    !---------------------
    head% path        =  path_file (output,                   &
                         path_file (bc_path,                  &
                         bc_filename(obstype, satid, instrid)))
    head% title       = 'Bias Correction Cycling File'
    head% institution = 'German Weather Service'
    head% source      = '3DVAR'
    head% created     =  cyyyymmddhhmm (run_time)
    head% modified    =  head% created
    head% first_date  =  cyyyymmddhhmm (ana_time)
    head% last_date   =  head% first_date

    if (present (t_decay)) then
      n = min (size(t_decay), size(head% t_decay))
      head% t_decay(1:n)   = t_decay(1:n)
      if (size(head% t_decay) > n) head% t_decay(n+1:) = head% t_decay(n)
    end if

  end subroutine new_bc_head

!------------------------------------------------------------------------------

  subroutine open_bc_write (h, obstype, satid, instrid)
  !-------------------------------------------------
  ! open bias correction file for writing
  ! and write the header to the NetCDF file.
  ! currently only used by AIRCRAFT bias correction.
  !-------------------------------------------------
  type(t_bcor_head) ,intent(inout)        :: h        ! file header data
  integer           ,intent(in)           :: obstype  ! 3dVar observation type
  integer           ,intent(in) ,optional :: satid    ! WMO satellite id
  integer           ,intent(in) ,optional :: instrid  ! WMO instrument id

    !--------------------------
    ! set some file header data
    !--------------------------
    h% exp         =  nex
    h% path        =  path_file (output,                   &
                         path_file (bc_path,                  &
                         bc_filename(obstype, satid, instrid)))
    h% title       = 'Bias Correction Cycling File'
    h% institution = 'German Weather Service'
    h% source      = '3DVAR'
    h% modified    =  cyyyymmddhhmm (run_time)
    h% last_date   =  cyyyymmddhhmm (ana_time)

    !---------------------
    ! open the NetCDF file
    !---------------------
    write(6,'(a)') repeat('-',79)
    write(6,'( )')
    write(6,'(a)')  '  Creating: '//trim(h% path)
    write(6,'( )')
    rname = 'open_bc_write'
    call chk (nf90_create (h% path, NF90_CLOBBER, h% ncid))

    !------------------------
    ! write global attributes
    !------------------------
    call chk(nf90_put_att (h% ncid, NF90_GLOBAL, "title",      h% title))
    call chk(nf90_put_att (h% ncid, NF90_GLOBAL, "institution",h% institution))
    call chk(nf90_put_att (h% ncid, NF90_GLOBAL, "source",     h% source))
    call chk(nf90_put_att (h% ncid, NF90_GLOBAL, "version",    h% version))
    call chk(nf90_put_att (h% ncid, NF90_GLOBAL, "created",    h% created))
    call chk(nf90_put_att (h% ncid, NF90_GLOBAL, "modified",   h% modified))
    call chk(nf90_put_att (h% ncid, NF90_GLOBAL, "first_date", h% first_date))
    call chk(nf90_put_att (h% ncid, NF90_GLOBAL, "last_date",  h% last_date))
    call chk(nf90_put_att (h% ncid, NF90_GLOBAL, "experiment", h% exp))
    if (all (h% t_decay(2:) == h% t_decay(1)) .or. &
        all (h% t_decay(2:) == 0._wp        )      ) then
      call chk(nf90_put_att (h% ncid, NF90_GLOBAL, "t_decay",  h% t_decay(1)))
    else
      call chk(nf90_put_att (h% ncid, NF90_GLOBAL, "t_decay",  h% t_decay))
    end if

  end subroutine open_bc_write

!------------------------------------------------------------------------------

  subroutine close_bc (h)
  !-------------------------------------------------
  ! close bias correction NetCDF file.
  ! currently only used by AIRCRAFT bias correction.
  !-------------------------------------------------
  type(t_bcor_head) ,intent(inout) :: h        ! file header data

    rname = 'close_bc'
    call chk (nf90_close (h% ncid))
    h% ncid = -1

  end subroutine close_bc

!------------------------------------------------------------------------------

  subroutine open_bc_read (h)
  !-------------------------------------------------
  ! open NetCDF bias correction file for reading
  ! and read the file header.
  ! currently only used by AIRCRAFT bias correction.
  !-------------------------------------------------
  type(t_bcor_head) ,intent(inout) :: h

    !---------------------
    ! open the NetCDF file
    !---------------------
    write(6,'(a)') repeat('-',79)
    write(6,'( )')
    write(6,'(a)')  '  Opening: '//trim(h% path)
    write(6,'( )')
    rname = 'open_bc_read'
    call chk (nf90_open (h% path, NF90_NOWRITE, h% ncid))

    !-----------------------
    ! read global attributes
    !-----------------------
    call chk(nf90_get_att (h% ncid, NF90_GLOBAL, "title",      h% title))
    call chk(nf90_get_att (h% ncid, NF90_GLOBAL, "institution",h% institution))
    call chk(nf90_get_att (h% ncid, NF90_GLOBAL, "source",     h% source))
    call chk(nf90_get_att (h% ncid, NF90_GLOBAL, "version",    h% version))
    call chk(nf90_get_att (h% ncid, NF90_GLOBAL, "created",    h% created))
    call chk(nf90_get_att (h% ncid, NF90_GLOBAL, "modified",   h% modified))
    call chk(nf90_get_att (h% ncid, NF90_GLOBAL, "first_date", h% first_date))
    call chk(nf90_get_att (h% ncid, NF90_GLOBAL, "last_date",  h% last_date))
    call chk(nf90_get_att (h% ncid, NF90_GLOBAL, "experiment", h% exp))
    call chk(nf90_get_att (h% ncid, NF90_GLOBAL, "t_decay",    h% t_decay))

  end subroutine open_bc_read

!------------------------------------------------------------------------------

  subroutine chk (status)
    integer          ,intent(in)           :: status
    !------------------------------------------------------------------------
    ! checks error status after each call to NetCDF routines,
    ! prints out text message each time an error code is returned and aborts.
    !------------------------------------------------------------------------
      if(status /= NF90_NOERR) then
        call finish (trim(rname), trim(nf90_strerror(status)))
      end if
  end subroutine chk

!==============================================================================
end module mo_biasc_io
