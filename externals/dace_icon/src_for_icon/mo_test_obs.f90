!
!+ Test the integrity of observation data (derived type t_obs and components)
!
MODULE mo_test_obs
!
! Description:
! Procedure to test for the integrity of observation data
! (derived type t_obs and components).
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
! V1_4         2009/03/26 Andreas Rhodin
!  deactivate checks for unused processors
! V1_5         2009/05/25 Harald Anlauf
!  test_spot_idx: slightly more cache friendly on NEC SX-9
! V1_6         2009/06/10 Andreas Rhodin
!  subroutine test_obs: test for pairs of wind observations
! V1_9         2010/04/20 Andreas Rhodin
!  subroutines test_obs_... : pass argument 'verb' to subsequent subroutines
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  skip test on subsequent wind components for obstype RAD
! V1_22        2013-02-13 Andreas Rhodin
!  make test dependent on module-type, not codetype
! V1_28        2014/02/26 Andreas Rhodin
!  improve diagnostic output
! V1_29        2014/04/02 Andreas Rhodin
!  ensure consistent model column input to GPSRO operator
! V1_42        2015-06-08 Andreas Rhodin
!  adaptions for temporal interpolation
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2005-2007
!==============================================================================
  !-------------
  ! Modules used
  !-------------
  use mo_t_obs,      only: t_obs,          &!
                           t_spot,         &!
                           t_index,        &! observation derived types
                           TEMP             ! module type value
  use mo_exception,  only: finish           ! error exit routine
  use mo_mpi_dace,   only: dace,           &! MPI group info
                           p_and            ! MPI and
  use mo_fdbk_tables,only: varno,          &! variable numbers table
                           obstype,        &! observation type code table
                           VN_U, VN_V,     &! wind component codes
                           OT_RAD,         &! codetype value
                           init_fdbk_tables ! initialise tables
  use mo_t_table,    only: name_value       ! find name of table entry
                                            ! (for given value)
  use mo_dec_matrix, only: t_vector         ! vector argument derived type
  use mo_fdbk_tables,only: OT_GPSRO         ! GPS radio occultation report type
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: test_obs  ! test the integrity of observation data
  public :: fix_mcols ! enforce consistent GPSRO arguments (model-columns)

  !-----------
  ! Interfaces
  !-----------
  interface test_obs
    module procedure test_obs_1
    module procedure test_obs_0
    module procedure test_obs_s
  end interface
!==============================================================================
contains
!==============================================================================
  subroutine test_obs_1 (obs, text, verb)
  type (t_obs)     ,intent(in) :: obs(:) ! obs.data type array
  character(len=*) ,intent(in) :: text   ! comment to be printed
  integer          ,intent(in) :: verb   ! <0: do noting; >0 verbose

    integer :: ib

    if (verb<0)                  return
    if (verb>0 .and. dace% lpio) write(6,*) dace% pe, 'test_olev_1 :',text
    do ib=1,size(obs)
      if (verb>1) write(6,*) dace% pe, 'test_olev_1 :',text, ib
      call test_obs_0 (obs(ib), ib, text, verb)
    end do
  end subroutine test_obs_1
!------------------------------------------------------------------------------
  subroutine test_obs_0 (obs, ib, text, verb)
  type (t_obs)     ,intent(in) :: obs  ! obs.data type variable
  integer          ,intent(in) :: ib   ! obs.data box index
  character(len=*) ,intent(in) :: text ! comment to be printed
  integer          ,intent(in) :: verb ! <0: do nothing; >0 verbose

    integer :: is, ns

    if (verb<0)          return
    ns = obs% n_spot
    if (ns==0)           return
    if (obs% pe /= dace% pe) return
    if (verb>0 .and. dace% lpio) write(6,*) dace% pe, 'test_obs_0 :',text
!   !------------------------
!   ! check unique sort order
!   !------------------------
!   if (any (obs% spot (1:obs% n_spot-1)% hd% id      &
!         >= obs% spot (2:obs% n_spot  )% hd% id ))   &
!     call finish('test_obs','no unique sort order')
    !----------------------------
    ! test consistency of indices
    !----------------------------
#ifdef __ICON__
!DIR$ NOINLINE
#endif
    if (test_spot_idx (obs% spot(1:ns)% o, obs% n_obs,'obs',verb)) &
      call finish('test_obs: spot% o',text)
    if (test_spot_idx (obs% spot(1:ns)% i, obs% n_int,'int',verb)) &
      call finish('test_obs: spot% i',text)
!   if (test_spot_idx (obs% spot(1:ns)% l, obs% n_lev,'lev',verb)) &
!     call finish('test_obs: spot% l',text)
    if (test_spot_idx (obs% spot(1:ns)% p, obs% n_par,'par',verb)) &
      call finish('test_obs: spot% p',text)
    !------------
    ! check spots
    !------------
    do is = 1, ns
      if (verb>1) write(6,*) dace% pe, 'test_obs_0 : spot',is,text
      call test_obs_s (obs, obs% spot(is), ib, is, text, verb)
    end do
  end subroutine test_obs_0
!------------------------------------------------------------------------------
  subroutine test_obs_s (obs, spot, ib, is, text, verb)
  type (t_obs)     ,intent(in) :: obs  ! obs. data type variable
  type (t_spot)    ,intent(in) :: spot ! spot data type variable
  integer          ,intent(in) :: ib   ! obs.data box index
  integer          ,intent(in) :: is   ! spot index
  character(len=*) ,intent(in) :: text ! comment to be printed
  integer          ,intent(in) :: verb ! <0: do noting; >0 verbose

    integer :: i,j

    if (verb<0) return
    if (verb>0 .and. dace% lpio) write(6,*) dace% pe, 'test_olev_s :',text
!   !--------------------
!   ! test for valid olev
!   !--------------------
!   select case (spot% hd% modtype)
!   case (TEMP)
!     do i=1, spot% o% n
!       j=i+spot% o% i
!       if (obs% olev(j)==0) then
!         write(0,*)dace% pe,'ERROR: olev=0: ib, is, i,n,j',ib,is,i,spot% o% n,j
!         call finish('test_obs',text)
!       end if
!     end do
!   end select
    !------------------------------------
    ! test for pairs of wind observations
    !------------------------------------
    if (spot% hd% obstype == OT_RAD) return  ! for obstype = channel #
    do i=1, spot% o% n
      j=i+spot% o% i
      select case (obs% varno(j))
      case (VN_U)
        if(i==spot% o% n)         call error ('u,v mismatch (1)')
        if(obs% varno(j+1)/=VN_V) call error ('u,v mismatch (2)')
      case (VN_V)
        if(i==1)                  call error ('u,v mismatch (3)')
        if(obs% varno(j-1)/=VN_U) call error ('u,v mismatch (4)')
        if(obs% body (j-1)% use% state /= obs% body(j)% use% state) &
                                  call error ('u,v mismatch (5)')
      end select
    end do

  contains

    subroutine error (text)
    character(len=*) ,intent(in) :: text
      integer :: k
      call init_fdbk_tables
      write(6,*) dace% pe,'ERROR: test_spot_idx: ',text
      write(6,*) dace% pe,'       statid  = ',spot% ident, spot%statid
      write(6,*) dace% pe,'       obstype = ',spot% hd% obstype
      write(6,*) dace% pe,'       obstype = ',spot% hd% obstype, &
                           name_value (obstype, spot% hd% obstype)
      do k = spot% o% i + 1, spot% o% i + spot% o% n
        write(6,*) dace% pe,'       varno   = ',obs% varno(k),           &
                         name_value (varno,     obs% varno(k)),          &
                                                obs% body(k)% lev_typ,   &
                         name_value (varno, int(obs% body(k)% lev_typ)) ,&
                                                obs% body(k)% use% state,&
                                                obs% body(k)% plev
      end do
      call finish ('test_obs_s',text)
    end subroutine error

  end subroutine test_obs_s
!------------------------------------------------------------------------------
  function test_spot_idx (idx, m, test, verb) result (fault)
  logical                      :: fault
  type (t_index)   ,intent(in) :: idx (:) ! indices to be tested
  integer          ,intent(in) :: m       ! max.no.of components
  character(len=*) ,intent(in) :: test    ! tested array
  integer          ,intent(in) :: verb    ! <0: do noting; >0 verbose

    logical :: used (m)
    integer :: i, i1, in
    fault = .false.
    used  = .false.
    do i=1,size(idx)
      if (verb>1) write(6,*) dace% pe, 'test_spot_idx :',idx(i),m,test
      if (idx(i)% i == -1) cycle
      i1 = idx(i)% i + 1
      in = idx(i)% i + idx(i)% n
      if (idx(i)% i < -1) call error ('idx% i < -1')
      if (fault) return
      if (idx(i)% n < 0)  call error ('idx% n < 0')
      if (fault) return
      if (in > m) call error ('idx% n > m')
      if (fault) return
      if (any(used(i1:in))) call error ('double used index')
      if (fault) return
      used(i1:in) = .true.
    end do

  contains

    subroutine error (text)
    character (len=*) :: text
      write(0,*) dace% pe,'ERROR: test_spot_idx: ',text
      write(0,*) dace% pe,'       test  = ',test
      write(0,*) dace% pe,'       spot  =',i
      write(0,*) dace% pe,'       i,n   =',idx(i)% i,idx(i)% n
      write(0,*) dace% pe,'       i1,in =',i1,in
      write(0,*) dace% pe,'       m     =',m
      fault = .true.
    end subroutine error

  end function test_spot_idx
!==============================================================================
  subroutine fix_mcols (obs, x, task, same, comment)
  type (t_obs)     ,intent(in)            :: obs(:)  ! observation meta data
  type (t_vector)  ,intent(inout)         :: x       ! argument to test / fix
  integer          ,intent(in)            :: task    ! 1=fix 2=test 3=test+print
  logical          ,intent(out) ,optional :: same    ! result for task=2,3
  character(len=*) ,intent(in)  ,optional :: comment ! comment for printout
  !--------------------------------------------------------------------
  ! The radio occultation operator implementation (by Michael Gorbunov)
  ! uses an intermediate data structure to pass arguments to the
  ! observation operator. If distinct occultations (within a 3dvar-box)
  ! require the same model columns these columns are passed by
  ! the same variable. For this reason it is not possible to pass
  ! different values of these input variables to different
  ! observation operators, allthough different elements are
  ! provided by the 'interpolation space' variable in the 3D-Var-PSAS.
  !
  ! This routine checks   (task=2) if 'interpolation space' entries are
  ! consistent or ensures (task=1) that they are consisting (by
  ! replicating the content of the first column to subsequent ones.
  !--------------------------------------------------------------------
    integer :: ib          ! box index
    integer :: is          ! spot index
    integer :: m_spt       ! max number of columns used by an occultation
    integer :: n_gpsro     ! number of GPSRO reports in a box
    integer :: ig          ! index of GPSRO report in a box
    integer :: i1, in      ! range of column entries to copy
    integer :: j1, jn      ! range of column entries to copy to
    integer :: m           ! number of input variables / column
    integer :: i, j        ! report and column index of previous report
    integer :: k, l, n     ! indices
    character(len=64) :: c ! comment line
    integer ,allocatable :: imcol (:,:,:) ! index array for cross-references

    !---------------------------
    ! process optional arguments
    !---------------------------
    c = ''; if (present (comment)) c = comment
    if (present(same)) same = .true.
    !----------------------------------------
    ! loop over boxes with GPSRO reports only
    !----------------------------------------
    do ib = 1, size (obs)
      if (obs(ib)% pe /= dace% pe) cycle
      !+++++++++++++++++++++++++++++++++++++++++++++
      ! endless compilation with cray compiler 8.2.1
      !+++++++++++++++++++++++++++++++++++++++++++++
!     n_gpsro = count (obs(ib)% spot(:)% hd% obstype == OT_GPSRO)
      n_gpsro = 0
      do i = 1,obs(ib)% n_spot
        if (obs(ib)% spot(i)% hd% obstype == OT_GPSRO) n_gpsro = n_gpsro + 1
       end do

      if (n_gpsro == 0) cycle
      !-------------------------------
      ! allocate cross-reference array
      !-------------------------------
      m_spt = maxval (obs(ib)% spot(:)% n_spt,                 &
                 mask=obs(ib)% spot(:)% hd% obstype == OT_GPSRO)
      allocate (imcol (3, m_spt, n_gpsro))
      !--------------------------------------
      ! loop over reports, GPSRO reports only
      !--------------------------------------
      imcol = 0
      ig    = 0
      m     = 0
      do is = 1, obs(ib)% n_spot
        if (obs(ib)% spot(is)% hd% obstype /= OT_GPSRO) cycle
        ig = ig + 1
        !------------------------
        ! set up cross-references
        !------------------------
        imcol (1,1:obs(ib)% spot(is)% n_spt,ig) = obs(ib)% spot(is)% imcol% imc(1)
        imcol (2,1                         ,ig) = is
        !---------------------------------------------------------
        ! 1st report: determine number of input variables / column
        !---------------------------------------------------------
        if (m==0) then
          n = obs(ib)% spot(is)% i% n
          l = obs(ib)% spot(is)% n_spt
          m = n / l
        endif
        !---------------------------------
        ! loop over columns of this report
        !---------------------------------
lk:     do k = 1, obs(ib)% spot(is)% n_spt
          !------------------------------------------------
          ! loop over columns of all previous GPSRO reports
          !------------------------------------------------
          do i = 1, ig - 1
            do j = 1, m_spt
              !-----------------------------------------------
              ! update cross-references if same column is used
              !-----------------------------------------------
              if (imcol (1,j,i) == obs(ib)% spot(is)% imcol(k)% imc(1)) then
                imcol (2,2,ig) = imcol (2,1,i)
                imcol (3,k,ig) = j
                cycle lk
              endif
            end do
          end do
        end do lk
      end do
      !------------------------------
      ! again loop over GPSRO reports
      !------------------------------
      if (task > 0) then
        do ig = 1, n_gpsro
          is = imcol (2,1,ig)
          !---------------------------------
          ! loop over columns (doubles only)
          !---------------------------------
          do k = 1, obs(ib)% spot(is)% n_spt
            if (imcol (3,k,ig) == 0) cycle
            !------------------------------
            ! check for correct column size
            !------------------------------
            n = obs(ib)% spot(is)% i% n
            l = obs(ib)% spot(is)% n_spt
            if (m * l /= n) call finish ('fix_mcols','m * l /= n')
            !-------------------------------
            ! set indices of previous report
            !-------------------------------
            i = imcol (2,2,ig)
            j = imcol (3,k,ig)
            i1 = obs(ib)% spot(is)% i% i + (k-1) * m + 1
            in = obs(ib)% spot(is)% i% i +  k    * m
            j1 = obs(ib)% spot(i )% i% i + (j-1) * m + 1
            jn = obs(ib)% spot(i )% i% i +  j    * m
            !------------------------------------
            ! actually test or fix column entries
            !------------------------------------
            select case (task)
            case (1)
              x% s(ib)% x(i1:in) = x% s(ib)% x(j1:jn)
            case (2,3)
              if (all (x% s(ib)% x(i1:in) == x% s(ib)% x(j1:jn))) then
                if (task==3) &
                  write(6,*) '### GPSRO column check: ==',ib,ig,i,k,j,trim(c)
              else
                if (present(same)) same = .false.
                if (task==3) &
                  write(6,*) '### GPSRO column check: /=',ib,ig,i,k,j,trim(c)
              endif
            end select
          end do
        end do
      endif
      !----------------------------
      ! deallocate cross-references
      !----------------------------
      deallocate (imcol)
    end do
    if (present(same)) same = p_and (same)

  end subroutine fix_mcols
!==============================================================================
end module mo_test_obs
