!
!+ Verification of forecasts vs. observations: low-level routines
!
MODULE mo_t_veri
!
! Description:
!   Verification of forecasts vs. observations: low-level routines
!
! Current Code Owner: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! DACE-2.01    2019-05-22 Harald Anlauf
!  Split-off from mo_veri_obs
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!-----------------------------------------------------------------------

!=============
! Modules used
!=============
!------------------------
! general purpose modules
!------------------------
use mo_kind,         only: wp, sp            ! working precision kind parameter
use mo_exception,    only: finish            ! abort in case of error
use mo_mpi_dace,     only: dace              ! MPI group info
use mo_dace_string,  only: char4,           &! convert int -> char(len=4)
                           char5             ! convert int -> char(len=5)
use mo_time,         only: t_time,          &! date+time derived type
                           operator(-),     &! subtract time
!                          operator(+),     &! add      time
!                          operator(==),    &! compare times
                           ihhhmm,          &! derive integer representation
                           cyyyymmddhhmm,   &! derive string from time
                           cyyyymmddhhmmss   ! derive string from time
use mo_run_params,   only: output,          &! output files path
                           path_file,       &! concatenate path+file name
                           model,           &! 'ICON', 'COSMO', ...
                           method,          &! 'ENKF', 'LETKF', ...
!                          ana_time,        &! analysis reference time
!                          fc_ref_time,     &! start of verification period
                           run_time          ! program run time
!-----------------------
! observation processing
!-----------------------
use mo_t_obs,        only: t_obs,           &! observation data type
!                          t_spot,          &! component of t_obs
                           source,          &! list of source-files
                           n_source,        &! number of source files
                           m_source,        &! max. number of source files
                           construct,       &! constructor routine for t_obs
                           destruct,        &! destructor  routine for t_obs
                           FT_FEEDBACK,     &! flag for feedback file
!                          FT_MISSING,      &! flag for missing  file
                           TSK_SET_CHR,     &!
                           TSK_SETUP_OP,    &!
                           TSK_SETUP_COLS    !
!                          TSK_SETUP_FULL,  &!
!                          TSK_INIT,        &!
!                          TSK_SHRINK,      &!
!                          TSK_READ,        &!
!                          TSK_R,           &!
!                          TSK_Y             ! flag to process_obs
use mo_obs_set,      only: t_obs_set         ! observation data type set
use mo_obs,          only: process_obs       ! general purpose routine
use mo_test_obs,     only: test_obs          ! test integrity of pbs.data
use mo_set_matrix,   only: PB_NONE,         &! no calculation of Pb flag
                           set_flags         ! set flags in t_obs
use mo_dec_matrix,   only: t_vector,        &! decomposed vector type
                           t_ivector,       &! integer    vector type
                           t_dec_info,      &! decomposition meta data
                           construct,       &! t_vector constructor routine
                           destruct,        &! t_vector  destructor routine
                           gather            ! gather vector components on PE
!                          scatter           ! scatter vector components
!---------------------
! model state handling
!---------------------
use mo_atm_grid,     only: t_grid,          &! grid meta data derived type
                           cmodel            ! model character representation
use mo_atm_state,    only: t_atm,           &! atmospheric state type
                           print             ! print atmospheric state
use mo_wmo_tables,   only: WMO6_LATLON,     &! latitude/Longitude         grid
                           WMO6_GAUSSIAN,   &! Gaussian                   grid
                           WMO6_ROTLL,      &! rotated latitude/longitude grid
                           DWD6_ICOSAHEDRON,&! icosahedral triangular     grid
                           DWD6_ICON         ! ICON (triangular)          grid
!-----------------------
! feedback file handling
!-----------------------
use mo_fdbk,         only: t_fdbk,          &! feedback file data type
                           t_fdbk_meta,     &! verification meta data type
                           open_fdbk_write, &! open feedback file for writing
!                          open_fdbk_read,  &! open feedback file for read access
                           read_meta,       &! read  attributes + meta data
                           add_history,     &! add history entry
                           add_verification,&! add verification data
                           write_history,   &! write history entry to file
                           get_veri_index,  &! get index of verification run
                           get_fillvalue,   &! get NetCDF fillvalue
                           close_fdbk,      &! close a feedback file
                           cleanup_fdbk      ! deallocate t_fdbk components
use mo_t_table,      only: text_value,      &! find name of table entry
!                          value_name,      &! find value of table entry
                           DLEN              ! lengt of text
use mo_fdbk_tables,  only: VT_FORECAST,     &! forecast             flag
                           VT_INIT_ANA,     &! initialized analysis flag
                           VT_ANALYSIS,     &!             analysis flag
                           VT_FIRSTGUESS,   &! first guess          flag
                           RC_ASS,          &! analysis cycle       flag
                           VE_DETERM,       &! deterministic run    flag
!                          VE_MEMBER,       &! ensemble             flag
                           VE_VQC_WEIGHT,   &! var. quality control weight
                           OF_BT_CLEAR_SKY, &! operator flag for clear-sky Tbs
!                          OT_RAD,          &! obstype RAD
                           ensmem,          &! member special meaning table
                           obstype           ! obstype table
use mo_t_netcdf_file,only: OPEN_W            ! file status flag
!---------------------
! netCDF f90 interface
!---------------------
use netcdf,          only: nf90_inq_varid,  &! get variable id by name
                           nf90_get_var,    &!
                           nf90_put_var,    &!
                           nf90_strerror,   &! returns NetCDF error message
                           NF90_NOERR        ! status return value: no error

implicit none

!================
! Public entities
!================
private
!------------
! subroutines
!------------
public :: setup_veri_obs    ! set up observation interpolation operators
public :: add_veri          ! add verification entry to feedback file
public :: open_veri_files   ! open  all feedback files (for appending)
public :: close_veri_files  ! close all feedback files
!---------------------
! (namelist) variables
!---------------------
public :: prefix_in         ! statistic/feedback input  file prefix
public :: prefix_out        ! statistic/feedback output file prefix
public :: rm_old            ! remove old verification data?

!===========
! interfaces
!===========
interface add_veri
  module procedure add_veri     ! scalar version
  module procedure add_veri_ens ! ensemble version
end interface add_veri

!=================
! Module variables
!=================
!--------------------
! namelist /veri_obs/
!--------------------
character(len=3)   :: prefix_out   ='ver'  ! statistic output file prefix
character(len=3)   :: prefix_in    ='cof'  ! statistic input  file prefix
integer            :: rm_old       = 0     ! remove old verification data:
                                           ! 0: do not remove
                                           ! 1: overwrite same runtype, vv, exp
                                           ! 2: overwrite same runtype, vv
                                           ! 3: remove all old data
                                           ! 4: overwrite data from diff. model,
                                           !    keeping orig. VQC weight entries

!=========================
! Derived type definitions
!=========================
!-----------------------------
! Container for feedback files
!-----------------------------
type t_fdbk_cntr
   character(256) :: fname  = ""        ! Feedback file full path
   type(t_fdbk)   :: fb                 ! Feedback file handle
   logical        :: opened = .false.   ! status
   logical        :: flush  = .true.    ! flush history at each update
!  integer        :: pe     = -1        ! handling processor
end type t_fdbk_cntr

type(t_fdbk_cntr), pointer :: fb_files(:) => NULL()

!==============================================================================
contains
!==============================================================================
  subroutine setup_veri_obs (obs, atm)
  type(t_obs_set) ,intent(inout) :: obs    ! observational data
  type(t_atm)     ,intent(in)    :: atm    ! forecast model
  !--------------------------------------------
  ! prepare observation interpolation operators
  !--------------------------------------------
    call process_obs (TSK_SETUP_COLS,           obs, atm, local=.true.)
    call process_obs (TSK_SET_CHR+TSK_SETUP_OP, obs, atm, local=.true.)
    call test_obs    (obs% o,'TSK_SETUP_CHR',0)
    call set_flags   (obs% o, atm% time, PB_NONE)
  end subroutine setup_veri_obs
!------------------------------------------------------------------------------
  subroutine add_veri (H_fc, obs, fc, meta, ensm, vv, op_flag)
  type (t_vector)   ,intent(in)           :: H_fc     ! model equivalents
  type(t_obs_set)   ,intent(in)           :: obs      ! observation data
  type(t_atm)       ,intent(in) ,optional :: fc       ! forecast model state
  type(t_fdbk_meta) ,intent(in) ,optional :: meta     ! verification meta data
  integer           ,intent(in) ,optional :: ensm     ! special meaning
  integer           ,intent(in) ,optional :: vv       ! nominal forecast time
  integer           ,intent(in) ,optional :: op_flag(:)!operator flag
  !---------------------------------------
  ! add verification data to feedback file
  ! (wrapper for deterministic state)
  !---------------------------------------
    type (t_vector) :: H_ens (1)

    H_ens (1) = H_fc              ! this is a flat copy, assignment(=) not used
    call add_veri_ens (H_ens, obs, fc, meta, ensm, vv, op_flag)

  end subroutine add_veri
!------------------------------------------------------------------------------
  subroutine add_veri_ens (H_ens, obs, fc, meta, ensm, vv, op_flag)
  type(t_vector)    ,intent(in)           :: H_ens(:) ! model equivalents
  type(t_obs_set)   ,intent(in)           :: obs      ! observation data
  type(t_atm)       ,intent(in) ,optional :: fc       ! forecast model state
  type(t_fdbk_meta) ,intent(in) ,optional :: meta     ! verification meta data
  integer           ,intent(in) ,optional :: ensm     ! special meaning
  integer           ,intent(in) ,optional :: vv       ! nominal forecast time
  integer           ,intent(in) ,optional :: op_flag(:)!operator flag
  !-------------------------------------------------------
  ! add verification data to feedback file
  ! (generic routine for deterministic and ensemble state)
  !-------------------------------------------------------
    integer                   :: i, j, l
    integer                   :: k
    type(t_fdbk)     ,pointer :: fb
    type(t_fdbk_cntr),pointer :: fb_cntr    ! Feedback file container entry
    integer                   :: varid
    integer                   :: io, i0
    integer                   :: ih, ib, is, isrc, n
    integer                   :: iv(1)
    integer                   :: nidx
    integer                   :: op_flag_
    real(sp) ,allocatable     :: buffer (:)
    integer  ,allocatable     :: ix(:), nix(:)
    integer                   :: status
    character(len=256)        :: outfile
    type (t_ivector)          :: on
    type (t_ivector)          :: oi
    type (t_ivector)          :: of
    type (t_ivector)          :: or
    type (t_ivector)          :: op
    type(t_dec_info) ,pointer :: dec_info   ! decomposition meta data
    type(t_dec_info) ,pointer :: dec_infb   ! decomposition meta data
    integer          ,save    :: n_veri (m_source) = 0
    logical                   :: preopened
    logical                   :: flush
    !-----------------------------------------------
    ! meta data derived from the model grid or state
    !-----------------------------------------------
    integer                   :: ke         ! member
    character(len=DLEN)       :: cmember    ! description: ' member IENS'
    type(t_fdbk_meta)         :: lmeta      ! verification meta data

    !------------------------------------------------------------
    ! optionally derive feedback file verification meta data from
    ! atmospheric grid and state
    !------------------------------------------------------------
    if (present (meta)) then
      lmeta = meta
    else if (present (fc)) then
      lmeta = fdbk_meta_atm (fc, ensm=ensm, vv=vv)
    else
      call finish ('add_veri_ens','neither "fc" nor "meta" is present')
    endif
    !------------------------------------
    ! special meaning for ensemble member
    !------------------------------------
    ke = VE_DETERM
    cmember = ' '
    if (present(ensm)) then
      ke = ensm
      cmember(2:) = text_value (ensmem, ensm)
    endif

    !---------------------------------------------------------------
    ! Check whether feedback files were preopened by open_veri_files
    !---------------------------------------------------------------
    preopened = associated (fb_files)
    if (preopened) then
       ! Sanity check
       if (size (fb_files) /= n_source) &
            call finish ("add_veri_ens","size (fb_files) /= n_source")
    else
       allocate (fb)
    end if

    !----------------------
    ! allocate index arrays
    !----------------------
    allocate       (dec_info)
    allocate       (dec_infb)
    call construct (dec_info, obs% oi% n,     obs% oi% n_b, &
                              obs% o% n_spot, obs% oi% b% pe)
    call construct (dec_infb, obs% oi% n,     obs% oi% n_b, &
                              obs% o% n_obs , obs% oi% b% pe)
    call construct (on, dec_info)
    call construct (oi, dec_info)
    call construct (of, dec_info)
    call construct (or, dec_info)
    call construct (op, dec_infb)

    !----------------------
    ! gather data on I/O PE
    !----------------------
    do ib = 1, H_ens(1)% n_s
      if (on% s(ib)% pe == dace% pe) then
        if (obs% o(ib)% n_obs > 0) then
          on% s(ib)%x = obs% o(ib)% spot% o% n
          oi% s(ib)%x = obs% o(ib)% spot% o% i
          of% s(ib)%x = obs% o(ib)% spot% hd% mon_file
          or% s(ib)%x = obs% o(ib)% spot% hd% mon_rec
          op% s(ib)%x = obs% o(ib)% body%     mon_pos
        endif
      endif
    end do
    call gather (H_ens, dest=dace% pio)
    call gather (on,    dest=dace% pio)
    call gather (oi,    dest=dace% pio)
    call gather (of,    dest=dace% pio)
    call gather (or,    dest=dace% pio)
    call gather (op,    dest=dace% pio)

    !-------------------------
    ! loop over feedback files
    !-------------------------
    do i = 1, n_source
      if (source(i)% filetype /= FT_FEEDBACK) cycle
      if (present(op_flag)) then
        if (op_flag(i) /= OF_BT_CLEAR_SKY) cycle
        op_flag_ = op_flag(i)
      else
        op_flag_ = 0
      end if
      if (preopened) then
        !-------------------
        ! use preopened file
        !-------------------
        fb_cntr => fb_files(i)
        fb      => fb_cntr% fb
        if (dace% lpio .and. .not. fb_cntr% opened)         &
          call finish ("add_veri_ens","file not preopened!?")
      else if (dace% lpio) then
        !----------
        ! open file
        !----------
        outfile  = path_file (output, source(i)% file)
        write(6,'(/a,a,i3)') &
          '  writing verification data to: ',trim(source(i)% file),i
        call open_fdbk_write (fb, outfile)
        if (fb% nc% status /= OPEN_W) &
          call finish ('add_veri_ens','cannot open: '//trim(outfile))
        !---------------
        ! read meta data
        !---------------
        call read_meta (fb)
      end if

      if (dace% lpio) then
        allocate (buffer (fb% n_body))
        allocate (ix     (fb% n_hdr ))
        allocate (nix    (fb% n_hdr ))
        status = nf90_inq_varid (fb% nc% ncid, 'i_body', varid)
        status = nf90_get_var   (fb% nc% ncid,                &
                                         varid,               &
                                         ix,                  &
                                         [1],                 &
                                         [fb% n_hdr]          )
        status = nf90_inq_varid (fb% nc% ncid, 'l_body', varid)
        status = nf90_get_var   (fb% nc% ncid,                &
                                         varid,               &
                                         nix,                 &
                                         [1],                 &
                                         [fb% n_hdr]          )
        !----------------------------
        ! write history and meta data
        !----------------------------
        if (preopened) then
           flush = fb_cntr% flush
        else
           flush = .true.
        end if
        call add_history                          &
          (fb,                                    &
           trim(lmeta% model)//' '//trim(method), &
           lmeta% initial_date,                   &
           lmeta% description,                    &
           cyyyymmddhhmm(run_time), write = flush )

        !----------------------------
        ! bookkeeping for rm_old == 3
        !----------------------------
        if (rm_old == 3) then
          fb% n_veri = 0
          n_veri(i)  = 0
        endif
        if (rm_old == -3) fb% n_veri = n_veri(i)
        !--------------------------------
        ! mark candidates for replacement
        ! keep   : VQC weight entries
        ! replace: different model
        !--------------------------------
        if (rm_old == 4) then
           do j = 1, min (fb% n_veri, size (fb% veri))
              if (any (fb% veri(j)% ens_member == [ VE_VQC_WEIGHT ])) cycle
              if (     fb% veri(j)% model      == model             ) cycle
              fb% veri(j)% ens_member = -999
           end do
        end if

        !---------------------------
        ! loop over ensemble members
        !---------------------------
        do k = 1, size (H_ens)
          if (size (H_ens) > 1) then
            ke      = k
            cmember = ' member '//char4(ke)
          endif

          !-------------------------
          ! replace or append data ?
          !-------------------------
          iv = 0
          select case (rm_old)
          case default
          case (1)
            if (lmeta% run_type==VT_ANALYSIS) then
              call get_veri_index (iv, nidx, fb,                 &
                                   run_type   = lmeta% run_type, &
                                   ens_member = ke,              &
                                   exp_id     = lmeta% exp_id    )
            else
              call get_veri_index (iv, nidx, fb,                         &
                                   run_type      = lmeta% run_type,      &
                                   forecast_time = lmeta% forecast_time, &
                                   ens_member    = ke,                   &
                                   exp_id        = lmeta% exp_id         )
            endif
          case (2)
            if (lmeta% run_type==VT_ANALYSIS) then
              call get_veri_index (iv, nidx, fb,                   &
                                   run_type      =lmeta% run_type, &
                                   ens_member    =ke )
            else
              call get_veri_index (iv, nidx, fb,                        &
                                   run_type      =lmeta% run_type,      &
                                   forecast_time =lmeta% forecast_time, &
                                   ens_member    =ke )
            endif
          case (4)
            if (lmeta% run_type==VT_ANALYSIS) then
              call get_veri_index (iv, nidx, fb,                   &
                                   run_type      =lmeta% run_type, &
                                   ens_member    =ke )
            else
              call get_veri_index (iv, nidx, fb,                        &
                                   run_type      =lmeta% run_type,      &
                                   forecast_time =lmeta% forecast_time, &
                                   ens_member    =ke )
            endif
            !--------------------------------------
            ! check for data marked for replacement
            !--------------------------------------
            if (iv(1) == 0) then
               do j = 1, min (fb% n_veri, size (fb% veri))
                  if (fb% veri(j)% ens_member == -999) then
                     fb% veri(j)% ens_member = ke
                     iv = j
                     exit
                  end if
               end do
            end if
          end select

          !----------------
          ! write meta data
          !----------------
          call add_verification                 &
            (fb,                                &! feedback file meta data
             lmeta% model,                      &! model (GME or ICON)
             lmeta% run_type,                   &! runtype
             lmeta% run_class,                  &! haupt=0,(vor=1),ass=2,test=3
             lmeta% initial_date,               &! initial_date
             lmeta% forecast_time,              &! forecast time (hhhmm)
             lmeta% resolution,                 &! resolution
             lmeta% domain_size,                &! domain
             trim(lmeta% description)//cmember, &! description
             ke,                                &! ensemble member
             lmeta% exp_id,                     &! experiment id
             varid,                             &! var.id. of 'veri_data'
             replace = iv(1),                   &! verification data index
             operator_flag=op_flag_)           ! operator flag

          if (fb%nc% error /= NF90_NOERR) &
            call finish('add_verification','NetCDF error')

          !---------
          ! add data
          !---------
          buffer = get_fillvalue (fb, 'veri_data')
          do ib = 1, size(obs% o)
            do is = 1, obs% o(ib)% n_spot
              ih   = or% s(ib)% x(is)  ! hd% record
              n    = on% s(ib)% x(is)  ! o% n
              i0   = oi% s(ib)% x(is)  ! o% i
              isrc = of% s(ib)% x(is)  ! hd% source
              if (isrc /= i) cycle
              io   = ix(ih)
              do j = 1, n
                l = op% s(ib)% x(i0+j)
                if (l > nix(ih)) call finish ('add_veri','index mismatch')
                if (l > 0)  buffer (io+l-1) = H_ens(k)% s(ib)% x(i0+j)
              end do
            end do
          end do

          status = nf90_put_var (fb% nc% ncid,     &
                                 varid,            &
                                 buffer,           &
                                 (/1,iv(1)/),      &
                                 (/fb% n_body, 1/) )
          if (status/=NF90_NOERR) call finish ('add_veri:gme_stat:',&
                                               nf90_strerror(status))
        end do

        !-----------------------
        ! deallocate temporaries
        !-----------------------
        deallocate (buffer)
        deallocate (ix)
        deallocate (nix)
      endif

      !-----------
      ! close file
      !-----------
      if (dace% lpio) then
        if (abs(rm_old) == 3) n_veri(i) = fb% n_veri
        if (.not. preopened) then
          call close_fdbk   (fb)
          call cleanup_fdbk (fb)
        end if
      endif
    end do

    !---------
    ! clean up
    !---------
    if (rm_old == 3) rm_old = -3
    call destruct (on)
    call destruct (oi)
    call destruct (of)
    call destruct (or)
    call destruct (op)
    call destruct (dec_info)
    deallocate    (dec_info)
    call destruct (dec_infb)
    deallocate    (dec_infb)

    if (.not. preopened) &
       deallocate (fb)

  end subroutine add_veri_ens
!------------------------------------------------------------------------------
  function fdbk_meta_atm (fc, grid, ensm, vv) result (meta)
  !-------------------------------------------------
  ! derive feedback file verification meta data from
  ! atmospheric grid and state
  !-------------------------------------------------
  type(t_atm)  ,intent(in)                   :: fc   ! atmospheric state
  type(t_grid) ,intent(in) ,target ,optional :: grid ! grid meta data
  integer      ,intent(in)         ,optional :: ensm ! special meaning
  integer      ,intent(in)         ,optional :: vv   ! nominal forecast time
  type(t_fdbk_meta)                          :: meta ! verification meta data

    type(t_grid) ,pointer :: g
    !-----------------------------------------------
    ! meta data derived from the model grid or state
    !-----------------------------------------------
    character(len=5)          :: c_hhhmm    ! forecast time character repr.
    integer                   :: i_hhhmm    ! forecast time integer   repr.
    integer                   :: vt         ! runtype integer representation
    character(len=40)         :: description
    real(sp)                  :: dx, dy     ! grid spacing
    integer                   :: nx, ny, nz ! grid point numbers
    character(len=8)          :: cm         ! model
    integer                   :: ke         ! ensemble member

    g => fc% grid
    if (present (grid)) g =>  grid
    if (.not. associated (g)) call finish('fdbk_meta_atm','grid must be present')
    !--------------------------------
    ! set run_time, class,
    ! change run_type for first guess
    !--------------------------------
    i_hhhmm = ihhhmm (fc% time - fc% ref_time)
    if (present (vv)) i_hhhmm = vv
    c_hhhmm = char5  (i_hhhmm)
    if (c_hhhmm(1:1)=='0') c_hhhmm(1:1)=' '
    select case (fc% runtype)
    case ('analysis','nudging')
      vt      = VT_ANALYSIS
    case ('init_ana')
      vt      = VT_INIT_ANA
    case ('forecast')
      vt      = VT_FORECAST
    case default
      call finish ('fdbk_meta_atm','unsupported runtype: '//trim (fc% runtype))
    end select
    if (vt == VT_FORECAST .and. fc% runclass == RC_ASS .and.&
        (g% gridtype == DWD6_ICOSAHEDRON .or.         &
         g% gridtype == DWD6_ICON        .or.         &
         prefix_out  == "mof"                )        ) then
      vt          = VT_FIRSTGUESS
      description = 'first guess vv='//c_hhhmm
    else if (vt == VT_FORECAST) then
      description = fc% runtype//' vv='//c_hhhmm
    else
      description = fc% runtype
    endif
    !------------------------------------
    ! special meaning for ensemble member
    !------------------------------------
    ke = VE_DETERM
    if (present(ensm)) then
      ke = ensm
    endif
    !---------------
    ! grid meta data
    !---------------
    if (g% model == 0) then
      cm = model
    else
      cm = cmodel (g% model)
    endif
    select case (g% gridtype)
    case (WMO6_LATLON, WMO6_ROTLL, WMO6_GAUSSIAN)
      nx = g% nx
      ny = g% ny
      nz = g% nz
      dx = g% di
      dy = g% dj
    case (DWD6_ICOSAHEDRON)
      nx = g% ni
      ny = g% ni
      nz = g% nz
      dx = 63.486_sp/nx
      dy = 63.486_sp/ny
    case (DWD6_ICON)
      nx = g% ni
      ny = g% ni
      nz = g% nz
      dx = 45.416_sp/nx
      dy = 45.416_sp/ny
    case default
      write(0,*) "gridtype =", g% gridtype
      call finish('fdbk_meta_atm','unsupported gridtype')
    end select
    !---------------------------------------------------------
    ! finally set verification meta data derived type variable
    !---------------------------------------------------------
    meta% model         = cm
    meta% run_type      = vt
    meta% run_class     = fc% runclass
    meta% initial_date  = cyyyymmddhhmm (fc% ref_time)
    meta% forecast_time = i_hhhmm
    meta% resolution    = [dx, dy]
    meta% domain_size   = [nx, ny, nz]
    meta% description   = description
    meta% ens_member    = ke
    meta% exp_id        = fc% expid
!   meta% operator_flag

  end function fdbk_meta_atm
!==============================================================================
  subroutine open_veri_files (sync)
    logical, optional      :: sync        ! Update history synchronously?
    !----------------
    ! Local variables
    !----------------
    logical                    :: lsync   ! Local copy of sync
    integer                    :: i
    type(t_fdbk_cntr), pointer :: fb_cntr ! Feedback file container entry
    type(t_fdbk),      pointer :: fb      ! Feedback file handle
    type(t_fdbk)               :: empty   ! Initialization

    lsync = .false.; if (present (sync)) lsync = sync

    if (associated (fb_files)) then
      call finish ("open_veri_files","feedback files already opened")
    end if
    allocate (fb_files(n_source))
    !-------------------------
    ! loop over feedback files
    !-------------------------
    do i = 1, n_source
      fb_cntr         => fb_files(i)
      fb_cntr% opened = .false.
      if (source(i)% filetype /= FT_FEEDBACK) cycle
      fb_cntr% fname  = path_file (output, source(i)% file)
      fb_cntr% fb     = empty
      if (dace% lpio) then
        !----------
        ! open file
        !----------
        write(6,'(a,i3,a,a)') &
          ' opening verification file', i, " : ", trim(source(i)% file)
        fb => fb_cntr% fb
        call open_fdbk_write (fb, fb_cntr% fname)
        if (fb% nc% status /= OPEN_W) &
          call finish ('open_veri_files','cannot open: '//trim(fb_cntr% fname))
        fb_cntr% opened = .true.
        fb_cntr% flush  = lsync
        !---------------
        ! read meta data
        !---------------
        call read_meta (fb)
      end if
    end do
  end subroutine open_veri_files
!==============================================================================
  subroutine close_veri_files ()
    !----------------
    ! Local variables
    !----------------
    integer                    :: i
    type(t_fdbk_cntr), pointer :: fb_cntr ! Feedback file container entry
    type(t_fdbk),      pointer :: fb      ! Feedback file handle

    if (.not. associated (fb_files)) then
      call finish ("close_veri_files","no open feedback files")
    end if
    if (n_source > 0 .and. dace% lpio) &
      write(6,'(/a,/)') ' closing verification files'
    do i = 1, n_source
      fb_cntr => fb_files(i)
      if (fb_cntr% opened) then
        fb => fb_cntr% fb
        if (.not. fb_cntr% flush) &     ! Flush history if not yet done
          call write_history (fb)
        call close_fdbk      (fb)
        call cleanup_fdbk    (fb)
        fb_cntr% opened = .false.
      endif
    end do
    deallocate (fb_files)
    nullify (fb_files)
  end subroutine close_veri_files
!==============================================================================
end MODULE mo_t_veri
