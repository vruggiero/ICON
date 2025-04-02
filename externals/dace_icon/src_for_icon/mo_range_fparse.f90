!+ generalization of formula parser
!
module mo_range_fparse
!
! Description:
!   Genralization of mo_fparse, which allows different functions for different
!   ranges (of the input variables).
!
! Current Code Owner: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    email: robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2019/07/24 Robin Faulwetter
!  Adapted fparser to DACE
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Robin Faulwetter           2019   initial version
!==============================================================================

  use kind_parameters,   only: wp

  use mo_exception,      only: finish

  use mo_dace_string,    only: tolower

  use data_constants,    only: pi

  use mo_mpi_dace,       only: p_bcast,             &
                               dace

  use mo_fparse,         only: init,                &!
                               evaluate,            &!
                               destruct,            &!
                               t_fparse,            &!
                               deg2rad,             &!
                               invalid,             &!
                               lv,                  &!
                               lt,                  &!
                               follows_bracket,     &!
                               p_bcast,             &!
                               err_msg,             &!
                               assignment(=)


  implicit none
  private

  public :: t_range_fparse
  public :: init
  public :: evaluate
  public :: destruct
  public :: construct
  public :: p_bcast
  public :: assignment(=)

  integer,          parameter :: mxlen = 1000
  integer,          parameter :: mxrng = 10
  real(kind=wp),    parameter :: pi_2 = 0.5_wp * pi
  character(len=6), parameter :: c_prev = 'y_prev' ! Reserved variable name, represents
                                                   ! the result from the previous range

  ! Type with all information about a "range", i.e. a region in the space spanned
  ! by all variables. It contains the function to be used there as well as the
  ! boundaries
  type t_range
    character(len=lt)          :: func      = ''
    type(t_fparse)             :: fp
    integer,           pointer :: i_v(:)    => null() ! index of var in t_range_fparse%vnames
    character(len=lt)          :: rstr      = ''
    integer                    :: n         = 0
    character(len=lv), pointer :: x_name(:) => null() ! name of the dimension/variable
    integer,           pointer :: i_x(:)    => null() ! index of x in t_range_fparse%vnames
    real(kind=wp),     pointer :: x_s(:)    => null() ! start          of range
    real(kind=wp),     pointer :: x_e(:)    => null() ! end            of range
    real(kind=wp),     pointer :: x_b(:)    => null() ! boundary width of range
  end type t_range

  ! Main type of this module. Contains a set of regions and a list of all the required
  ! variables
  type t_range_fparse
    character(len=mxlen)       :: func      =  ''
    integer                    :: nvar      =  0      ! #required variables
    character(len=lv), pointer :: vnames(:) => null() ! name of the dimension/variable
    !TODO: what is the purpose of the following ???
    real(kind=wp),     pointer :: v(:)      => null() ! variable values (used in evaluation)
    integer                    :: nr        =  0      ! #ranges
    type(t_range),     pointer :: r(:)      => null() ! ranges
  end type t_range_fparse


  interface construct
    module procedure construct_range
    module procedure construct_range_fparse
  end interface

  interface destruct
    module procedure destruct_range
    module procedure destruct_range_fparse
  end interface

  interface init
    module procedure init_range_fparse
  end interface

  interface evaluate
    module procedure evaluate_range_fparse
  end interface

  interface p_bcast
    module procedure bcast_range
    module procedure bcast_range_fparse
  end interface

  interface assignment(=)
    module procedure assign_range
    module procedure assign_range_fparse
  end interface

#if (defined (__GFORTRAN__) && (__GNUC__ >= 10)) || defined (NAGFOR)
    !----------------------------------------------
    ! include interfaces for external bcast routine
    !----------------------------------------------
    interface
       subroutine p_bcast_derivedtype (buffer, count, source, comm)
         type(*) ,intent(inout)     :: buffer       ! variable to bcast
         integer ,intent(in)        :: count        ! len(byte) of variable
         integer ,intent(in)        :: source       ! source processor index
         integer ,intent(in)        :: comm         ! communicator
       end subroutine p_bcast_derivedtype
    end interface
#endif

contains

  subroutine init_range_fparse(func, variablenames, trf, status, used_vnames)
    character(len=*),     intent(in)            :: func
    character(len=*),     intent(in)            :: variablenames(:)
    type(t_range_fparse), intent(out), target   :: trf
    integer,              intent(out)           :: status
    logical,              intent(in),  optional :: used_vnames

    logical                        :: used_vnames_
    integer                        :: n_y
    ! Split into different function for different ranges
    integer                        :: i_s  (mxrng)  ! start index of function
    integer                        :: i_e  (mxrng)  ! end   index of function
    integer                        :: i_r_s(mxrng)  ! start index of range definition
    integer                        :: i_r_e(mxrng)  ! end index of range definition
    integer                        :: lenf, lenr
    integer                        :: nr, ir
    integer                        :: stat
    integer                        :: i, j, k, i1, i2, iop
    logical                        :: mask(mxrng)
    ! Analyze a single range-string
    character(len=37), parameter   :: vchars = 'abcdefghijklmnopqrstuvwxyz0123456789_'
    character(len=26), parameter   :: chars  = 'abcdefghijklmnopqrstuvwxyz'
    character(len=3),  parameter   :: rop    = '=:;' ! operators in range string
    character(len=lv), allocatable :: varnames(:)
    type(t_range),     pointer     :: pr => null()
    integer                        :: nword, nv, lv_
    character(len=lv)              :: word(100)
    integer                        :: i_read   ! -1:none 0:x_name ; 1: x_s; 2:x_e ; 3:x_b
    character(len=8), parameter    :: c_read(3) = (/'start   ', 'end     ', 'boundary'/)
    real(kind=wp)                  :: xdum

    status = 13

    if (present(used_vnames)) then
      used_vnames_ = used_vnames
    else
      used_vnames_ = .false.
    end if

    if (len_trim(func) > mxlen) call finish('init_range_fparse', 'func "'//trim(func)//'" too long.')
    ! trf%func = tolower(trim(adjustl(func)))
    ! Workaround necessary for CRAY compiler:
    trf%func = trim(adjustl(func))
    trf%func = tolower(trf%func)
    lenf = len_trim(trf%func)

    i_s = 0
    i_e = 0
    i_r_s = 0
    i_r_e = 0

    ! Split into functions for different ranges
    nr = 1
    i_s(nr) = 1
    i = 0
    do while (i < lenf - 6)
      i = i + 1
      if (i > 1) then
        if (     trf%func(i-1:i-1) /= ' ') cycle
      end if
      if (verify(trf%func(i  :i  ), 'rR') > 0) cycle
      if (verify(trf%func(i  :i  ), 'rR') > 0) cycle
      if (verify(trf%func(i+1:i+1), 'aA') > 0) cycle
      if (verify(trf%func(i+2:i+2), 'nN') > 0) cycle
      if (verify(trf%func(i+3:i+3), 'gG') > 0) cycle
      if (verify(trf%func(i+4:i+4), 'eE') > 0) cycle
      i1 = follows_bracket(trf%func,lenf,i+5,set='([{')
      if (i1 <= 0) call finish('init_range_fparse', 'missing "(" after "range"')
      i2 = scan(trf%func(i1+1:),')]}')
      if (i2 <= 0) call finish('init_range_fparse', 'missing ")" after "range"')
      i2 = i1 + i2
      if (i-1 > i_s(nr)) then
        if (trim(trf%func(i_s(nr):i-1)) /= '') then
          i_e(nr) = i-1
          nr = nr + 1
        end if
      end if
      i_r_s(nr) = i1 + 1
      i_r_e(nr) = i2 - 1
      i_s  (nr) = i2 + 1
      i = i2 + 1
    end do
    i_e(nr) = lenf

    mask(1:nr) = (i_e(1:nr)-i_s(1:nr)+1 > 0)
    i = count(mask(1:nr))
    i_s(1:i) = pack(i_s(1:nr), mask=mask(1:nr))
    i_e(1:i) = pack(i_e(1:nr), mask=mask(1:nr))
    nr = i

    allocate(trf%r(nr))
    trf%nr = nr

    ! Fill individual t_range entries
    do ir = 1, nr
      pr => trf%r(ir)
      ! Function to be evaluated
      pr%func = trim(trf%func(i_s(ir):i_e(ir)))
      nv = size(variablenames)
      allocate(varnames(nv+1))
      do i = 1, nv
        if (len_trim(variablenames(i)) > lv) call finish('init_range_fparse', &
                     'too long variable name "'//trim(variablenames(i))//'".')
        if (trim(variablenames(i)) == c_prev) call finish('init_range_fparse', &
                     'variable name "'//c_prev//'" not allowed. It is reserved &
                     &for the result from the previous range.')
        varnames(i) = trim(variablenames(i))
      end do
      varnames(nv+1) = c_prev
      call init(pr%func, varnames, stat, tfp=pr%fp, used_vnames=.true.)
      if (stat /= 0) call finish('init_range_fparse', &
           'failed to initialize function "'//trim(pr%func)//'": '//&
           trim(err_msg(stat)))
      deallocate(varnames)
      ! range(s)
      if (i_r_s(ir) > 0) then
        pr%rstr = trim(trf%func(i_r_s(ir):i_r_e(ir)))
        lenr = len_trim(pr%rstr)
        if (pr%rstr(lenr:lenr) /= ';') then
          ! Add ';' at the end in order to facilitate the analysis of pr%rstr below
          pr%rstr = trim(pr%rstr)//';'
          lenr = lenr + 1
        end if
        ! count number of variables/constraints, set up arrays
        nv = 0
        do i = 1, lenr
          if (pr%rstr(i:i) == '=') nv = nv + 1
        end do
        allocate(pr%x_name(nv), pr%x_s(nv), pr%x_e(nv), pr%x_b(nv))
        pr%x_name = ''
        pr%x_s    = -huge(0._wp)
        pr%x_e    =  huge(0._wp)
        pr%x_b    =  0._wp
        ! Determine details
        i1 = 1
        i = 1
        nv = 0
        i_read = -1
        do i = 1, lenr
          iop = index(rop, pr%rstr(i:i))
          select case (iop)
          case(1)
            nv = nv + 1
            if (i_read /= -1) call finish('init_range_fparse', &
                   'unable to interpret "'//trim(pr%rstr)//'". Missing ";"?')
            read(pr%rstr(i1:i-1),*,iostat=stat) pr%x_name(nv)
            if (stat /= 0 .or. pr%x_name(nv)=='') call finish('init_range_fparse', &
                 'failed to determine variable name from "'//trim(pr%rstr)//'".')
            i_read = 0
            i1 = i+1
          case(2,3)
            i_read = i_read + 1
            if (i > i1) then
              read(pr%rstr(i1:i-1),*,iostat=stat) xdum
              if (stat /= 0) call finish('init_range_fparse', &
                   'failed to determine range '//c_read(i_read)//' from "'//trim(pr%rstr)//'".')
            else
              xdum = invalid
            end if
            select case(i_read)
            case(:0)
              call finish('init_range_fparse', &
                   'missing variable name in "'//trim(pr%rstr)//'".')
            case(1)
              if (xdum /= invalid) pr%x_s(nv) = xdum
            case(2)
              if (xdum /= invalid) pr%x_e(nv) = xdum
            case(3)
              if (xdum /= invalid) then
                pr%x_b(nv) = xdum
                pr%x_b(nv) = 0.5_wp * pr%x_b(nv) ! Half of the width is more
                                                 ! convenient in the evaluation
              end if
            case(4:)
              call finish('init_range_fparse', &
                   'junk at the end of "'//trim(pr%rstr)//'".')
            end select
            if (iop == 3) i_read = -1
            i1 = i+1
          end select
        end do
        pr%n = nv
      else
        pr%n = 0
      end if ! i_r_s(ir) > 0
    end do

    if (used_vnames_) then
      ! Produce the list of all required variables
      nv = 0
      do ir = 1, trf%nr
        pr => trf%r(ir)
        ! Required for function evaluation
        do j = 1, pr%fp%nvar
          if (pr%fp%varnames(j) /= c_prev .and. &
              .not.any(word(1:nv) == pr%fp%varnames(j))) then
            nv = nv + 1
            word(nv) = pr%fp%varnames(j)
          end if
        end do
        ! Required for region evaluation
        do j = 1, pr%n
          if (.not.any(word(1:nv) == pr%x_name(j))) then
            nv = nv + 1
            word(nv) = pr%x_name(j)
          end if
        end do
      end do
      allocate(trf%vnames(nv))
      do i = 1, nv
        trf%vnames(i) = trim(word(i))
      end do
    else
      ! Use the supplied list of variablenames
      nv = size(variablenames)
      allocate(trf%vnames(nv))
      do i = 1, nv
        trf%vnames(i) = trim(variablenames(i))
      end do
    end if
    trf%nvar = nv

    ! Determine indices pointing to the required variables
    n_y = 0
    do ir = 1, trf%nr
      pr => trf%r(ir)
      ! Required for function evaluation
      allocate(pr%i_v(pr%fp%nvar))
      do j = 1, pr%fp%nvar
        if (trim(pr%fp%varnames(j)) == c_prev) then
          n_y = 1
          pr%i_v(j) = nv + 1
        else
          do k = 1, nv
            if (trf%vnames(k) == pr%fp%varnames(j)) then
              pr%i_v(j) = k
              EXIT
            endif
          end do
        end if
      end do
      ! Required for region evaluation
      allocate(pr%i_x(pr%n))
      do j = 1, pr%n
        do k = 1, nv
          if (trf%vnames(k) == pr%x_name(j)) then
            pr%i_x(j) = k
            EXIT
          endif
        end do
      end do
    end do
    allocate(trf%v(nv + n_y))

    status = 0

  end subroutine init_range_fparse


  function evaluate_range_fparse(vars, trf) result(answer)
    real(kind=wp)                       :: answer
    real(kind=wp),        intent(in)    :: vars(:)
    type(t_range_fparse), intent(inout) :: trf

    real(kind=wp), parameter :: eps = epsilon(0._wp)
    real(kind=wp)            :: fact, y_, y_old, v
    type(t_range), pointer   :: pr => null()
    integer                  :: ir, ix
    answer = invalid
    if (trf%nvar /= size(vars)) RETURN
    trf%v(1:trf%nvar) = vars(:)

    y_old = 0._wp
    range_loop: do ir = 1, trf%nr
      pr => trf%r(ir)
      fact = 1._wp
      x_loop: do ix = 1, pr%n
        v = vars(pr%i_x(ix))
        if (v < pr%x_s(ix)-pr%x_b(ix) .or. &
            v > pr%x_e(ix)+pr%x_b(ix)) CYCLE range_loop
        if (pr%x_b(ix) > eps) then
          if (v < pr%x_s(ix)+pr%x_b(ix)) then
            ! On start boundary
            fact = fact * boundary_smoother(v-pr%x_s(ix))
          end if
          if (v > pr%x_e(ix)-pr%x_b(ix)) then
            ! On end boundary
            fact = fact * boundary_smoother(pr%x_e(ix)-v)
          end if
        end if
      end do x_loop
      if (any(pr%i_v(:) == trf%nvar+1)) trf%v(trf%nvar+1) = y_old  ! Result from previous ranges, variable "y"
      y_ = evaluate(trf%v(pr%i_v(:)), pr%fp)
      answer = fact * y_ + (1._wp -fact) * y_old
      y_old = answer
    end do range_loop

  contains

    function boundary_smoother(x) result(f)
      real(kind=wp) :: f
      real(kind=wp), intent(in) :: x
      f = 0.5_wp * (1._wp + sin(x/pr%x_b(ix) * pi_2))
    end function boundary_smoother

  end function evaluate_range_fparse


  elemental subroutine construct_range(tr)
    type(t_range), intent(out) :: tr
  end subroutine construct_range


  elemental subroutine destruct_range(tr)
    type(t_range), intent(inout) :: tr
    call destruct(tr%fp)
    if (associated(tr%i_v   )) deallocate(tr%i_v   )
    if (associated(tr%i_x   )) deallocate(tr%i_x   )
    if (associated(tr%x_name)) deallocate(tr%x_name)
    if (associated(tr%x_s   )) deallocate(tr%x_s   )
    if (associated(tr%x_e   )) deallocate(tr%x_e   )
    if (associated(tr%x_b   )) deallocate(tr%x_b   )
  end subroutine destruct_range


  elemental subroutine construct_range_fparse(trf)
    type(t_range_fparse), intent(out) :: trf
  end subroutine construct_range_fparse


  elemental subroutine destruct_range_fparse(trf)
    type(t_range_fparse), intent(inout) :: trf

    if (associated(trf%r     )) then
      call destruct(trf%r(:))
      deallocate(trf%r)
    end if
    if (associated(trf%vnames)) deallocate(trf%vnames)
    if (associated(trf%v     )) deallocate(trf%v     )
  end subroutine destruct_range_fparse


  subroutine bcast_range (r, source, comm)
    type (t_range)   ,intent(inout) :: r
    integer          ,intent(in)    :: source
    integer ,optional,intent(in)    :: comm

    integer :: lcom
    integer :: count = 0

    lcom = dace% comm ;if(present(comm)) lcom = comm

    if (count == 0) count = size(transfer(r,(/' '/)))

    if (dace% pe /= source) call destruct (r)
    call p_bcast_derivedtype (r, count, source, lcom)

    call p_bcast (r%fp, source, lcom)

    if (dace% pe /= source) then
      if (associated(r% i_v   )) allocate(r% i_v   (r%fp%nvar))
      if (associated(r% x_name)) allocate(r% x_name(r%n      ))
      if (associated(r% i_x   )) allocate(r% i_x   (r%n      ))
      if (associated(r% x_s   )) allocate(r% x_s   (r%n      ))
      if (associated(r% x_e   )) allocate(r% x_e   (r%n      ))
      if (associated(r% x_b   )) allocate(r% x_b   (r%n      ))
    end if
    if (associated(r% i_v   )) call p_bcast(r% i_v   (1:r%fp%nvar), source, lcom)
    if (associated(r% x_name)) call p_bcast(r% x_name(1:r%n      ), source, lcom)
    if (associated(r% i_x   )) call p_bcast(r% i_x   (1:r%n      ), source, lcom)
    if (associated(r% x_s   )) call p_bcast(r% x_s   (1:r%n      ), source, lcom)
    if (associated(r% x_e   )) call p_bcast(r% x_e   (1:r%n      ), source, lcom)
    if (associated(r% x_b   )) call p_bcast(r% x_b   (1:r%n      ), source, lcom)

  end subroutine bcast_range


  subroutine bcast_range_fparse(rf, source, comm)
    type (t_range_fparse), intent(inout) :: rf
    integer,               intent(in)    :: source
    integer, optional,     intent(in)    :: comm

    integer :: lcom
    integer :: count = 0
    integer :: i

    lcom = dace% comm ;if(present(comm)) lcom = comm

    if (count == 0) count = size(transfer(rf,(/' '/)))

    if (dace% pe /= source) call destruct (rf)
    call p_bcast_derivedtype (rf, count, source, lcom)

    if (dace% pe /= source) then
      if (associated(rf% vnames)) allocate(rf% vnames(rf%nvar))
      if (associated(rf% v     )) allocate(rf% v     (rf%nvar))
      if (associated(rf% r     )) allocate(rf% r     (rf%nr  ))
    end if
    if (associated(rf% vnames)) call p_bcast(rf% vnames(1:rf%nvar), source, lcom)
    if (associated(rf% v     )) call p_bcast(rf% v     (1:rf%nvar), source, lcom)
    if (associated(rf% r     )) then
      do i = 1, rf%nr
        call p_bcast(rf% r(i), source, lcom)
      end do
    end if

  end subroutine bcast_range_fparse

  elemental subroutine assign_range(y, x)
    type(t_range), intent(inout) :: y
    type(t_range), intent(in)    :: x

    call destruct(y)
    y%func = x%func
    y%fp   = x%fp    ! see assign_fparse@mo_fparser.f90
    y%rstr = x%rstr
    y%n    = x%n
    if (associated(x%i_v)) then
      allocate(y%i_v(size(x%i_v)))
      y%i_v = x%i_v
    end if
    if (associated(x%x_name)) then
      allocate(y%x_name(size(x%x_name)))
      y%x_name = x%x_name
    end if
    if (associated(x%i_x)) then
      allocate(y%i_x(size(x%i_x)))
      y%i_x = x%i_x
    end if
    if (associated(x%x_s)) then
      allocate(y%x_s(size(x%x_s)))
      y%x_s = x%x_s
    end if
    if (associated(x%x_e)) then
      allocate(y%x_e(size(x%x_e)))
      y%x_e = x%x_e
    end if
    if (associated(x%x_b)) then
      allocate(y%x_b(size(x%x_b)))
      y%x_b = x%x_b
    end if
  end subroutine assign_range

  elemental subroutine assign_range_fparse(y, x)
    type(t_range_fparse), intent(inout) :: y
    type(t_range_fparse), intent(in)    :: x

    call destruct(y)
    y%func = x%func
    y%nvar = x%nvar
    y%nr   = x%nr
    if (associated(x%vnames)) then
      allocate(y%vnames(size(x%vnames)))
      y%vnames = x%vnames
    end if
    if (associated(x%v)) then
      allocate(y%v(size(x%v)))
      y%v = x%v
    end if
    if (associated(x%r)) then
      allocate(y%r(size(x%r)))
      y%r = x%r   ! See above: assign_range
    end if
  end subroutine assign_range_fparse

end module mo_range_fparse
