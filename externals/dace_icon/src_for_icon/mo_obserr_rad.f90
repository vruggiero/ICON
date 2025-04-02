!
!+ correlated observation errors for radiances
!
MODULE mo_obserr_rad
!
! Description:
!   correlated observation errors for radiances
!
! Current Code Owner: DWD, Olaf Stiller
!    phone: +49 69 8062 2910
!    fax:   +49 69 8062 3721
!    email: olaf.stiller@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================

  !=============
  ! modules used
  !=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_exception, only: finish           ! abort routine
  use mo_kind,      only: wp, sp           ! working, single precision
  use mo_namelist,  only: position_nml,   &! position namelist
                          nnml,           &! namelist Fortran unit number
                          POSITIONED       ! ok    code from position_nml
  use mo_mpi_dace,  only: dace,           &! MPI group info
                          p_bcast          ! broadcast routine
  use mo_run_params,only: data,           &! path name of data  files
                          path_file        ! concatenate: path/file.sufx
  use netcdf,       only: nf90_open,      &! open   NetCDF file
                          nf90_close,     &! close  NetCDF file
                          NF90_NOWRITE     ! NetCDF read flag
  use mo_t_netcdf,  only: ncid,           &! NetCDF file id
                          stanc,          &! NetCDF start  parameter
                          counc,          &! NetCDF count  parameter
                          strnc,          &! NetCDF stride parameter
                          chk,            &! checks error status
                          get_dim,        &! read NetCDF dimension
                          get_var          ! read NetCDF variable
  use mo_matrix,    only: check_rs,       &! check for eigenvalues
                          scale            ! scale matrix for diagonal -> 1
  use mo_instrid,   only: instr_rttov      ! WMO instrument id from RTTOV number
  implicit none

!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: read_obserr_rad_nml  ! read namelist /obserr_rad/
  public :: obs_corr_rad         ! full R matrix for different instruments
  public :: n_obserr_rad         ! # of entries in obs_corr_rad

!------------------------------------------------------------------------------
  !-----------------
  ! module variables
  !-----------------
  type t_obserr_rad
    integer               :: instrid                ! WMO instrument ID
    character(len=128)    :: file
    real(wp)              :: w_min_r                ! minimal eigenvalue for R matrix
    real(wp)              :: sgm_vq                 ! variational quality control bound
    integer               :: frm_vq                 ! var. quality control formulation
    integer               :: n
    real(wp) ,allocatable :: corr (:,:)             ! correlation
    real(wp) ,allocatable :: chan (:)               ! channel
    logical  ,allocatable :: use_chan(:)            ! whether the channel shall be considered as
                                                    ! correlated with others
  end type t_obserr_rad

  integer ,parameter         :: m_obserr_rad = 10           ! max. number of correlation matrices
  integer                    :: n_obserr_rad =  0           ! actual number of corr. matrices
  type(t_obserr_rad), target :: obs_corr_rad (m_obserr_rad) ! correlation matrices
  integer ,parameter         :: m_chan_range = 10           ! max. number of channel ranges

!==============================================================================
contains
!==============================================================================
  subroutine read_obserr_rad_nml
  !===========================
  ! read namelist /OBSERR_RAD/
  !===========================

    !----------------------
    ! namelist /OBSERR_RAD/
    !----------------------
    integer               :: instr     ! RTTOV instrument Id
    integer               :: instrid   ! WMO instrument Id
    character(len=128)    :: file      ! file with estimated covariances
    real(wp)              :: w_min_r   ! minimum eigenvalue to keep
    real(wp)              :: sgm_vq    ! variational quality control bound
    integer               :: frm_vq    ! var. quality control formulation
    integer               :: start_chan(m_chan_range)
    integer               :: end_chan  (m_chan_range)
    logical               :: use_vari ! whether the channel shall be considered as

    namelist  /OBSERR_RAD/ instr, instrid, file, w_min_r, sgm_vq, frm_vq, start_chan, end_chan, &
                           use_vari

    !----------------
    ! local variables
    !----------------
    type(t_obserr_rad), pointer     :: oep => null()
    integer                         :: ierr
    integer                         :: ii
    logical                         :: first
    character(len=256)              :: path = ''
    logical                         :: exists
    integer                         :: n_lev
    real(wp),           allocatable :: corr_sym (:,:) ! temporary
    real(wp),           allocatable :: vari (:)       ! temporary
    integer                         :: hits           ! mumber of eigenvalues modified
    real(wp)                        :: evmin          ! smallest eigenvalue in matrix
    real(wp)                        :: evmax          ! largest  eigenvalue in matrix

    first = .true.
    if (dace% lpio) then
      write (6,*)
      write (6,*) '  reading namelist /OBSERR_RAD/'
      write (6,*)
    endif
    do
      !-------------
      ! set defaults
      !-------------
      instr      = -1      ! RTTOV instrument Id
      instrid    = -1      ! WMO instrument Id
      file       = ''      ! file with estimated covariances
      w_min_r    = 0.1_wp  ! minimum eigenvalue to keep
      sgm_vq     = -1._wp  ! variational quality control bound
      frm_vq     = -1      ! var. quality control formulation
      start_chan = -1      ! start channel
      end_chan   = -1      ! end channel
      use_vari   = .false.

      !---------------------------------
      ! read namelist, consistency check
      !---------------------------------
      if (dace% lpio) then
        call position_nml ('OBSERR_RAD' ,lrewind=first ,status=ierr)
        first=.false.
        select case (ierr)
        case (POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=OBSERR_RAD, iostat=ios)
          if(ios/=0)call finish('read_obserr_rad_nml','ERROR in namelist /OBSERR_RAD/')
#else
          read (nnml ,nml=OBSERR_RAD)
#endif
        end select
      endif
      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ierr, dace% pio)
      if (ierr /= POSITIONED) exit
      !------------------------
      ! broadcast to other PE's
      !------------------------
      call p_bcast (instr,      dace% pio)
      call p_bcast (instrid,    dace% pio)
      call p_bcast (file,       dace% pio)
      call p_bcast (w_min_r,    dace% pio)
      call p_bcast (sgm_vq,     dace% pio)
      call p_bcast (frm_vq,     dace% pio)
      call p_bcast (start_chan, dace% pio)
      call p_bcast (end_chan,   dace% pio)
      call p_bcast (use_vari,   dace% pio)
      !------------
      ! store entry
      !------------
      n_obserr_rad = n_obserr_rad + 1
      if (n_obserr_rad > m_obserr_rad)                                 &
        call finish('read_obserr_rad_nml','n_obserr_rad > m_obserr_rad')
      oep => obs_corr_rad (n_obserr_rad)

      if (instrid < 0 .and. instr >= 0) instrid = instr_rttov(instr, -1)
      oep% instrid = instrid
      oep% file    = file
      oep% w_min_r = w_min_r
      oep% sgm_vq  = sgm_vq
      oep% frm_vq  = frm_vq
      if (dace% lpio) then
        write (6,*) '    instrid   = ' ,instrid
        write (6,*) '    file      = ' ,trim(file)
        write (6,*) '    w_min_r   = ' ,w_min_r
        write (6,*) '    sgm_vq    = ' ,sgm_vq
        write (6,*) '    frm_vq    = ' ,frm_vq
        write (6,*) '    use_vari  = ' ,use_vari
        !-----------------------
        ! read covariance matrix
        !-----------------------
        path = path_file (data, file)
        inquire (file=trim (path), exist=exists)
        if (.not. exists) &
          call finish('read_obserr_rad_nml','can not open '//trim(path))

        call chk (nf90_open (path, NF90_NOWRITE, ncid))
        call get_dim (n_lev ,'lev')

        oep% n = n_lev
        allocate (oep% chan (n_lev))
        allocate (oep% corr (n_lev,n_lev))

        stanc = 1
        strnc = 1
        counc = n_lev

        call get_var (oep% chan,'lev')
        call get_var (oep% corr,'cov_oa_of')

        call chk(nf90_close (ncid))
      endif
      call p_bcast (oep% n    ,dace% pio)
      if (.not.dace% lpio) then
        allocate (oep% chan (oep% n))
        allocate (oep% corr (oep% n,oep% n))
      endif
      call p_bcast (oep% chan ,dace% pio)
      call p_bcast (oep% corr ,dace% pio)

      !--------------------------------------------------------
      ! check, which channels shall be considered as correlated
      !--------------------------------------------------------
      allocate (oep% use_chan (oep% n))
      if (any(start_chan > 0) .or. any(end_chan > 0)) then
         oep% use_chan(:) = .false.
         do ii = 1, size(start_chan)
            if (start_chan(ii) > 0 .or. end_chan(ii) > 0) then
               if (end_chan(ii) <= 0) end_chan(ii) = oep% chan(oep%n)
               where(oep% chan(1:oep%n) >= start_chan(ii) .and. oep% chan(1:oep%n) <= end_chan(ii)) &
                    oep% use_chan(1:oep%n) = .true.
            end if
         end do
      else
         oep% use_chan(:) = .true.
      end if
      do ii=1, oep% n
         if (.not.oep% use_chan(ii)) then
            oep% corr(:ii-1,ii) = 0._wp
            oep% corr(ii+1:,ii) = 0._wp
            oep% corr(ii,:ii-1) = 0._wp
            oep% corr(ii,ii+1:) = 0._wp
         end if
      enddo
      !----------------------------------------------------------------
      ! symetrize, scale for correlation matrix, set minimal eigenvalue
      !----------------------------------------------------------------
      allocate (corr_sym (oep% n, oep% n) )
      allocate (vari (oep% n))
      do ii=1, oep% n
        corr_sym(:,ii)=(oep% corr(:,ii)+oep% corr(ii,:))/2
        vari(ii) = sqrt(max(oep% corr(ii,ii),0._wp))
      enddo

      call scale    (corr_sym)  ! normalize for better interpretation of w_min_r
      call check_rs (corr_sym , min_ev= w_min_r, y= oep% corr, &
                                hits= hits, evmin= evmin, evmax= evmax                 )
      if(use_vari) then
        call scale(oep% corr, vari,vari)
      else
        call scale    (oep% corr)  ! normalize for correct correlation matrix
      endif
      if (dace% lpio) then
        write (6,*) "============ BEGIN: variance of covariance matrix ================"
        write (6,*)  "index chan.  old variance  new variance"
        do ii=1, oep% n
          write (6,'(2(1x,I5),2(1x,F13.9),2x,A)') ii, int(oep% chan(ii)),  vari(ii)**2, oep% corr(ii,ii)
        enddo
        write (6,*) "============ END: variance of covariance matrix ================"
      endif

      deallocate    (corr_sym)
      deallocate    (vari)

      !---------------------------------
      ! some printout on modified matrix
      !---------------------------------
      if (dace% lpio) then
        write (6,*) '    min. eigenvalue           = ' ,evmin
        write (6,*) '    max. eigenvalue           = ' ,evmax
        write (6,*) '    # of eigenvalues modified = ' ,hits
        write (6,*) '    # correlated channels     = ' ,count(oep% use_chan(1:oep%n))
        write (6,*) '    # uncorrelated channels   = ' ,count(.not.oep% use_chan(1:oep%n))
        write (6,*)
      endif
    end do

  end subroutine read_obserr_rad_nml
!==============================================================================
end module mo_obserr_rad
