!
!+ Routines to write atmospheric fields to GRADS file
!
MODULE mo_grads_atm
!
! Description:
!   Routines to write atmospheric fields to GRADS file
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_23        2013-03-26 Andreas Rhodin
!  moved routines 'to_grads' from mo_memory to new module mo_grads_atm
!  allow to write GRADS files for distributed lat/lon fields
! V1_26        2013/06/27 Andreas Rhodin
!  new subroutine from_grads: read atmospheric state from GRADS file
! V1_31        2014-08-21 Andreas Rhodin
!  to_grads: fix for globally allocated fields, new optional parameter iostat
! V1_42        2015-06-08 Andreas Rhodin
!  minor cleanup
! V1_45        2015-12-15 Harald Anlauf
!  Cleanup towards TR15581/F2003 compatibility
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================
#include "tr15581.incf"
!-------------
! Modules used
!-------------
  use mo_kind,       only : wp             ! working precision kind parameter
  use mo_memory,     only : t_m,          &! atmospheric field derived type
                            allocate       ! allocate pointer component
  use mo_mpi_dace,   only : dace,         &! DACE communication info
                            p_bcast        ! generic MPI broadcast routine
  use mo_grads,      only : t_ctl,        &! GRADS metadata derived type
                            write_var,    &! write variable  to   GRADS file
                            read_var,     &! read  variable  from GRADS file
                            record         ! get record number in GRADS file
  use mo_atm_decomp, only : t_atm_dec      ! atm. field decomposition info
  use mo_atm_transp, only : gather_multi, &! gather  multi level field
                            scatter_multi  ! scatter multi level field
  implicit none
!----------------
! Public entities
!----------------
  private
  public :: to_grads               ! write atmospheric field to   GRADS file
  public :: from_grads             ! read  atmospheric field from GRADS file

!-----------
! Interfaces
!-----------
  interface to_grads
    module procedure m0_to_grads
    module procedure m_to_grads
  end interface to_grads

  interface from_grads
    module procedure m0_from_grads
    module procedure m_from_grads
  end interface from_grads

contains
!==============================================================================
  subroutine m_from_grads (ctl, m, dec, name, t, yrev)
  !------------------------------------------------
  ! read list of atmospheric fields from GRADS file
  !------------------------------------------------
  type (t_ctl)         ,intent(in)           :: ctl     ! GRADS meta data
  type (t_m)           ,intent(inout)        :: m(:)    ! list of atm. fields
  type (t_atm_dec)     ,intent(in) ,optional :: dec     ! decomposition info
  character(len=*)     ,intent(in) ,optional :: name    ! pass name to GRADS
  integer              ,intent(in) ,optional :: t       ! time slice
  logical              ,intent(in) ,optional :: yrev    ! latitudes flipped
    integer :: i
    do i = 1,size(m)
      call from_grads (ctl, m(i), dec, name, t=t, yrev=yrev)
    end do
  end subroutine m_from_grads
!------------------------------------------------------------------------------
  subroutine m0_from_grads (ctl, m0, dec, name, t, yrev)
  !---------------------------------------
  ! read atmospheric field from GRADS file
  !---------------------------------------
  type (t_ctl)         ,intent(in)           :: ctl     ! GRADS meta data
  type (t_m)           ,intent(inout) TARGET :: m0      ! atmospheric field
  type (t_atm_dec)     ,intent(in) ,optional :: dec     ! decomposition info
  character(len=*)     ,intent(in) ,optional :: name    ! pass name to GRADS
  integer              ,intent(in) ,optional :: t       ! time slice
  logical              ,intent(in) ,optional :: yrev    ! latitudes flipped

    !----------------
    ! local variables
    !----------------
    logical               :: lpar          ! parallel mode, read from 1 pe only
    logical               :: lpio          ! read from this pe
    integer               :: irec          ! >0 if field is present in file
    integer               :: lb(2), ub(2)  ! bounds of full field
    character(len=32)     :: nam           ! name to search for in file
    real(wp) ,ALLOCATABLE                 &!
                   TARGET :: ptr (:,:,:,:) ! temporary to store full field
    real(wp),     pointer :: p   (:,:,:,:) ! auxiliary pointer to full field

    !------------------------
    ! check for parallel mode
    !------------------------
    lpio = .true.
    lpar = .false.
    if (present (dec)) then
      if (dec% nproc1 * dec% nproc2 /= 1) then
        lpar = .true.
        lpio = (dace% lpio)
      endif
    endif

    !--------------
    ! read variable
    !--------------
    if (lpio) then
      !-------------------------------
      ! check for presence of variable
      !-------------------------------
      if (present (name)) then
        nam = trim(name)//'.'//m0%i% name
      else
        nam = m0%i% name
      endif
      irec = record (ctl, nam, t)

      !--------------
      ! read if found
      !--------------
      if (irec/=0) then
        lb (1) = dec% ilim1 (0)
        ub (1) = dec% ilim1 (dec% nproc1) - 1
        lb (2) = dec% ilim2 (0)
        ub (2) = dec% ilim2 (dec% nproc2) - 1
        allocate (ptr (      lb(1) :       ub(1), &
                             lb(2) :       ub(2), &
                       m0%i% lb(3) : m0%i% ub(3), &
                       m0%i% lb(4) : m0%i% ub(4) ))
        call read_var (ctl, ptr(:,:,:,1), nam,&
                       t=t, yrev=yrev         )
      endif
    endif

    !----------------------------
    ! distribute result and store
    !----------------------------
    if (lpar) then
      call p_bcast (irec, dace% pio)
      if (irec/=0) then
        call allocate (m0) ! ,nlev)
        p => ptr
        call scatter_multi (p, m0% ptr (:,:,:,:), dec, dace% pio)
        if (lpio) deallocate (ptr)
      endif
    else
      if (irec/=0) then
        if (ALLOCATED (m0% ptr)) deallocate (m0% ptr)
        MOVE_ALLOC (ptr, m0% ptr)   ! m0% ptr => ptr
        m0%i% alloc = .true.
      endif
    endif

  end subroutine m0_from_grads
!==============================================================================
  subroutine m_to_grads (ctl, m, dec, name, t, comment, yrev, iostat)
  !-----------------------------------------------
  ! write list of atmospheric fields to GRADS file
  !-----------------------------------------------
  type (t_ctl)     ,intent(inout)         :: ctl     ! GRADS meta data
  type (t_m)       ,intent(in)            :: m(:)    ! list of atm. fields
  type (t_atm_dec) ,intent(in)  ,optional :: dec     ! decomposition info
  character(len=*) ,intent(in)  ,optional :: name    ! pass name to GRADS
  integer          ,intent(in)  ,optional :: t       ! time slice
  character(len=*) ,intent(in)  ,optional :: comment ! data set description
  logical          ,intent(in)  ,optional :: yrev    ! latitudes flipped
  integer          ,intent(out) ,optional :: iostat  ! I/O error status
    integer :: i
    if (present(iostat)) iostat = 0
    do i = 1,size(m)
      call to_grads (ctl, m(i), dec, name, t=t, comment=comment, &
                     yrev=yrev, iostat=iostat)
      if (present (iostat)) then
        if (iostat /= 0) exit
      endif
    end do
  end subroutine m_to_grads
!------------------------------------------------------------------------------
  subroutine m0_to_grads (ctl, m0, dec, name, t, comment, yrev, iostat)
  !--------------------------------------
  ! write atmospheric field to GRADS file
  !--------------------------------------
  type (t_ctl)     ,intent(inout)         :: ctl     ! GRADS meta data
  type (t_m)       ,intent(in)     TARGET :: m0      ! atmospheric field
  type (t_atm_dec) ,intent(in)  ,optional :: dec     ! decomposition info
  character(len=*) ,intent(in)  ,optional :: name    ! pass name to GRADS
  integer          ,intent(in)  ,optional :: t       ! time slice
  character(len=*) ,intent(in)  ,optional :: comment ! data set description
  logical          ,intent(in)  ,optional :: yrev    ! latitudes flipped
  integer          ,intent(out) ,optional :: iostat  ! I/O error status

    real(wp)  _POINTER :: ptr  (:,:,:,:)
!   real(wp) , pointer :: lptr (:,:,:,:)
    logical            :: lpio
    logical            :: lpara
    integer            :: lb(2), ub(2)

    if (present (iostat)) iostat = 0
    if (m0%i% alloc) then
!   if (m0%i% alloc .and. (m0%i% rep == 'gg' .or. m0%i% rep == 'ff')) then

      ptr  => m0% ptr(:,:,:,:)
      lpio  = .true.
      lpara = .false.
      !--------------------------------------------------
      ! gather field from multiple processors if required
      !--------------------------------------------------
      if (present (dec)) then
        if (dec% nproc1 * dec% nproc2 /= 1) then
          lpio  = (dace% lpio)
          lpara = .true.
          lb (1) = dec% ilim1 (0)
          ub (1) = dec% ilim1 (dec% nproc1) - 1
          lb (2) = dec% ilim2 (0)
          ub (2) = dec% ilim2 (dec% nproc2) - 1
          if (lb (1) /= m0% i% lb(1) .or. &
              lb (2) /= m0% i% lb(2) .or. &
              ub (1) /= m0% i% ub(1) .or. &
              ub (2) /= m0% i% ub(2)      ) then
            if (lpio) then
              allocate (ptr (      lb(1) :       ub(1), &
                                   lb(2) :       ub(2), &
                             m0%i% lb(3) : m0%i% ub(3), &
                             m0%i% lb(4) : m0%i% ub(4) ))
            endif
!           lptr => m0% ptr (dec% ilim1 (dec% my_num1): dec% ilim1 (dec% my_num1+1)-1,&
!                            dec% ilim2 (dec% my_num2): dec% ilim2 (dec% my_num2+1)-1,&
!                                      :, :                                         )
!           call gather_multi (ptr,    lptr (:,:,:,:), dec, dace% pio)
            call gather_multi (ptr, m0% ptr (:,:,:,:), dec, dace% pio)
          endif
        endif
      endif
      !---------------
      ! write to GRADS
      !---------------
      if (lpio) then
        if (present (name)) then
          call write_var (ctl, ptr(:,:,:,1), trim(name)//'.'//m0%i% name, &
                          t=t, comment=comment, yrev=yrev, iostat=iostat  )
        else
          call write_var (ctl, ptr(:,:,:,1), m0%i% name,                &
                          t=t, comment=comment, yrev=yrev, iostat=iostat)
        endif
        if (.not.associated (ptr, m0% ptr)) deallocate (ptr)
      endif
      if (present (iostat)) then
        if (lpara) call p_bcast (iostat, dace% pio)
      endif
    endif

  end subroutine m0_to_grads
!==============================================================================
end module mo_grads_atm
