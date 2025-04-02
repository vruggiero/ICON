!
!+ generic formula parser
!
module mo_fparse
!
! Description:
!   Generic formula parser. Original code by Ivomar Brito Soares,
!   https://github.com/ivomarb/Fortran-Expression-Evaluator
!   http://www.labfit.net/functionparser.htm
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
! Ivomar Brito Soares et al. ????   original code
! Robin Faulwetter           2019   adaptation to DACE and implementation of t_fparse
!==============================================================================

  !TODO: check whether everything is deconstructed
  !Allow to get only variable names with replacements and without full init (use in calc_tovs_cloud@mo_cloud)
  
  use kind_parameters,   only: wp

  use data_constants,    only: pi

  use mo_dace_string,    only: split2, tolower

  use mo_mpi_dace,       only: p_bcast, &
                               dace

  use mo_namelist,       only: position_nml,     &! routine to position nml group
                               nnml,             &! namelist fortran unit number
                               POSITIONED         ! position_nml: OK    return flag

  use mo_exception,      only: finish

  implicit none

  private

  public :: init
  public :: evaluate
  public :: evaluate_save
  public :: destruct
  public :: construct
  public :: read_nml_fparse_shortcuts
  public :: replace_shortcuts
  public :: t_fparse
  public :: deg2rad
  public :: invalid
  public :: lv
  public :: lt
  public :: follows_bracket
  public :: p_bcast
  public :: err_msg
  public :: assignment(=)

  real(kind=wp),     parameter :: deg2rad = pi/180._wp
  real(kind=wp),     parameter :: invalid = -huge(1._wp)
  integer,           parameter :: lv      = 30
  integer,           parameter :: lt      = 1000

  character(len=11), parameter :: numbers = '.0123456789'
  character(len=26), parameter :: chars = 'abcdefghijklmnopqrstuvwxyz'
  character(len=7),  parameter :: operators = '+-*/^()'
  

  ! Contains the results of the init subroutine, i.e. all the temporary information
  ! required to evaluate the function
  type t_fparse
    character(len=lt)          :: func          =  ''     ! original function expression
    integer                    :: nvar          =  0
    character(len=lv), pointer :: varnames(:)   => null() ! variablenames, that occur
    integer                    :: ioperations   =  0      ! number of operations
    integer,           pointer :: operations(:) => null() ! operations
    integer                    :: numberk       =  1      ! #numbers in the function
    real(kind=wp),     pointer :: number(:)     => null() ! numbers in the function
    real(kind=wp),     pointer :: pdata(:)      => null() ! temporary data used within the evaluation
  end type t_fparse

  type(t_fparse),save,target      :: tmp
  type(t_fparse),     pointer     :: tmp_ => null()

  integer                         :: nvar
  character(len=lv),  pointer     :: varnames(:)
  integer                         :: nvar_used
  character(len=lv),  pointer     :: varnames_used(:)
  character(len=lv),  pointer     :: varnames_used_(:)
  character(len=lt),  allocatable :: stokens(:)
  integer,            pointer     :: operations(:)
  integer                         :: ntokens       =  0
  character,          pointer     :: opaddsub(:)   => null()
  integer                         :: isaddsub      =  1
  character,          pointer     :: opmuldiv(:)   => null()
  integer                         :: ismuldiv      =  1
  character(len=lt)               :: toke          =  ''
  integer                         :: itoke         =  1
  integer                         :: ioperations   =  1
  integer                         :: numberk       =  1
  real(kind=wp),      pointer     :: pdata(:)      => null()
  real(kind=wp),      pointer     :: number(:)     => null()
  integer                         :: status_fparse =  0
  logical                         :: used_vnames_  = .false.
  logical,            pointer     :: use_var(:)    => null()
  character(len=lt)               :: err_msg_add   =  ''

  ! Shortcuts
  ! Chortcuts might be defined in fparse_shortcut namelists.
  ! Function strings are scanned for shortcut names and names are
  ! replaced by the corresponding expressions.
  integer, parameter :: len_name    = 100
  integer, parameter :: len_expr    = 1000
  integer, parameter :: mx_shortcut = 100
  type t_fparse_shortcut
    character(len=len_name) :: name = ''
    character(len=len_expr) :: expr = ''
  end type t_fparse_shortcut
  type(t_fparse_shortcut) :: shortcuts(mx_shortcut)
  integer                 :: n_shortcut = 0

  
  interface construct
    module procedure construct_fparse
  end interface

  interface destruct
    module procedure destruct_fparse1
    module procedure destruct_fparse2
  end interface

  interface init
    module procedure init
  end interface

  interface evaluate
    module procedure evaluate_fparse
  end interface

  interface p_bcast
    module procedure bcast_fparse
  end interface

  interface assignment(=)
    module procedure assign_fparse
  end interface

contains

  subroutine read_nml_fparse_shortcuts
    character(len=*), parameter :: proc = 'read_nml_fparse_shortcut'
    character(len=100)          :: msg  = ''
    logical                 :: lfirst
    integer                 :: ierr, istat, inml
    ! namelist
    character(len=len_name) :: name = ''
    character(len=len_expr) :: expr = ''
    namelist /FPARSE_SHORTCUT/ name, expr
    
    lfirst = .true.
    inml   = 0
    nml_loop: do
      inml = inml + 1
      name = ''
      expr = ''
      if (dace% lpio) then
        call position_nml ('FPARSE_SHORTCUT' ,lrewind=lfirst ,status=ierr)
        if (ierr == POSITIONED) then
          write(msg, '("namelist ",I2)') inml
#if defined(__ibm__)
          read (nnml ,nml=FPARSE_SHORTCUT, iostat=istat)
          if (istat /= 0) call finish (proc, 'ERROR in /FPARSE_SHORTCUT/ '//trim(msg))
#else
          read (nnml ,nml=FPARSE_SHORTCUT)
#endif
          name = tolower(name)
          expr = tolower(expr)
          if (lfirst) then
            write(*,*)
            write(*,'(1x,A)') 'read /FPARSE_SHORTCUT/ namelists:'
          end if
          write(*,'(3x,A)') trim(msg)//':'
          write(*,'(5x,"name = ",A)') trim(name)
          write(*,'(5x,"expr = ",A)') trim(expr)
          istat = 0
          if (trim(name) == '') then
            istat = 1
            write(*,'(A)') '*** WARNING: invalid namelist: empty name!'
          endif
          if (trim(expr) == '') then
            istat = 1
            write(*,'(A)') '*** WARNING: invalid exprlist: empty expr!'
          endif
        end if
      end if
      lfirst = .false.
      call p_bcast(ierr, dace% pio)
      call p_bcast(istat,dace% pio)
      if (ierr  /= POSITIONED) exit  nml_loop
      if (istat /= 0         ) cycle nml_loop
      
      call p_bcast(name, dace% pio)
      call p_bcast(expr, dace% pio)
      if (trim(name) == trim(expr)) call finish(proc, 'Invalid recursive replacement')
      
      n_shortcut = n_shortcut + 1
      if (n_shortcut > mx_shortcut) call finish(proc, 'Too many shortcuts: reduce number of &
           &FPARSE_SHORTCUT namelists or increase mx_shortcut!')
      shortcuts(n_shortcut) = t_fparse_shortcut(name, expr)
    end do nml_loop
  end subroutine read_nml_fparse_shortcuts

  
  subroutine replace_shortcuts(func)
    character(len=*), intent(inout) :: func
    character(len=*), parameter     :: proc = 'replace_shortcuts'
    character(len=lt)               :: f1, f2
    integer                         :: i, is, ie, j
    integer                         :: lf, ltr, ldiff

    lf  = len(func)
    ltr = len_trim(func)
    i   = 1
    do while (i <= ltr)
      !It's a variable, or function  name
      if (index(chars,func(i:i)) /= 0) then
        is = i
        do while (i < ltr)
          if (index(operators, func(i+1:i+1)) /= 0 .or. func(i+1:i+1) == ' ') exit
          i = i + 1
        end do
        ie = i
        do j = 1, n_shortcut
          if (func(is:ie) == trim(shortcuts(j)%name)) then
            ldiff = len_trim(shortcuts(j)%expr) - (ie-is+1)
            if (ltr + ldiff > lf) then
              write(0,'("func=",A)') trim(func)
              write(0,'("shortcut%name=",A)') trim(shortcuts(j)%name)
              write(0,'("shortcut%expr=",A)') trim(shortcuts(j)%expr)
              write(0,'("is,ie=",2(1x,I4))') is, ie
              write(0,'("ltr,ldiff,lf=",3(1x,I4))') ltr,ldiff,lf
              call finish(proc, 'length with replacement will increase maximum length')
            end if
            if (is > 1) then
              f1 = func(1:is-1)
            else
              f1 = ''
            end if
            if (ie < ltr) then
              f2 = func(ie+1:ltr)
            else
              f2 = ''
            end if
            ltr = ltr + ldiff
            if (ltr > len(func)) call finish(proc,'function string length exceeded')
            func = trim(f1)//trim(shortcuts(j)%expr)//trim(f2)
            !i  = i  + ldiff
            i  = i  + 1
          end if
        end do
      else if (index(numbers,func(i:i)) /= 0) then
        do while (i < ltr)
          if ( (index(operators, func(i+1:i+1)) /= 0 &
                 .and..not.(verify(func(i:i),'EeDd')==0 .and. verify(func(i+1:i+1),'+-')==0)) &
               .or. func(i+1:i+1) == ' ') exit
          i = i + 1
        end do
      end if
      i = i + 1
    end do

  end subroutine replace_shortcuts
  
  
  subroutine init (func_, variablenames, status, tfp, used_vnames)
    !
    !        This subroutine shifts all characters of the function
    !        expression to lowercase and converts exponents signals ** to ^
    !
    character(len=*), intent(inout)                  :: func_
    character(len=*), intent(in)                     :: variablenames(:)
    integer,          intent(out)                    :: status
    type(t_fparse),   intent(out),  optional, target :: tfp
    logical,          intent(in),   optional         :: used_vnames

    character(len=lt)                :: func = ''
    integer                          :: i

    err_msg_add   =  ''
    status = 99

    if (present(used_vnames)) then
      used_vnames_ = used_vnames
    else
      used_vnames_ = .false.
    end if

    if (len(func_) > lt) then
      status = -4
      return
    end if
    func = tolower(func_)

    call replace_shortcuts(func)
    
    status = check_err(func)
    if (status /= 0) then
      write(0,*) 'func='//trim(func)
      RETURN
    end if
    call convert(func)
    call blanks(func)
    call recog_variables (func, variablenames)
    status = status_fparse
    if (status /= 0)  write(0,*) 'func='//trim(func)

    if (allocated(stokens)) then
      deallocate(stokens)
    end if
    if (associated(opaddsub)) then
      deallocate(opaddsub)
    end if
    if (associated(opmuldiv)) then
      deallocate(opmuldiv)
    end if

    if (present(tfp)) then
      tmp_ => tfp
    else
      tmp_ => tmp
    end if
    i = min(len_trim(func),len(tmp_%func))
    tmp_%func(1:i)   = func(1:i)
    if (used_vnames_) then
      tmp_%nvar      =  nvar_used     ; nvar_used      =  0
      tmp_%varnames  => varnames_used ; varnames_used  => null()
      nvar = 0
      deallocate(varnames)
    else
      tmp_%nvar      =  nvar        ; nvar        =  0
      tmp_%varnames  => varnames    ; varnames    => null()
    endif
    tmp_%ioperations =  ioperations ; ioperations =  0
    tmp_%operations  => operations  ; operations  => null()
    tmp_%numberk     =  numberk     ; numberk     =  0
    tmp_%number      => number      ; number      => null()
    tmp_%pdata       => pdata       ; pdata       => null()

  end subroutine init


  subroutine recog_variables (func, variablenames)
    !
    !       This subroutine recognizes the variables and set their values
    !
    character(len=*), intent(in)    :: variablenames(:)
    character(len=*), intent(inout) :: func
    integer :: i
    nvar = size(variablenames)
    allocate(varnames(nvar))
    if (used_vnames_) then
      nvar_used = 0
      allocate(varnames_used(nvar))
    end if
    do i = 1, nvar
      varnames(i) = trim(variablenames(i))
    end do
    call tokens_analyzer (func)

  end subroutine recog_variables


  subroutine tokens_analyzer (func)
    !
    !       This subroutine scans the func string storing its basic elements
    !
    character(len=*), intent(in) :: func
    integer                      :: k = 1, i = 1
    integer                      :: lenf, k0
    integer                      :: irightbrackets = 1, ileftbrackets = 1
    logical                      :: status = .true.

    k = 1
    i = 1
    irightbrackets = 1
    ileftbrackets = 1
    ntokens = 0
    status = .true.
    lenf = len_trim(func)

    do while (k <= lenf)
      !It's a variable, or function  name
      if (index(chars,func(k:k)) /= 0) then
        do while (k < lenf)
          if (index(operators, func(k+1:k+1)) /= 0 .or. func(k+1:k+1) == ' ') exit
          k = k + 1
        end do
        k = k + 1
        ntokens = ntokens + 1

        !It's a number
      else if (index(numbers,func(k:k)) /= 0) then
        do while (k < lenf)
          if ( (index(operators, func(k+1:k+1)) /= 0 &
                 .and..not.(verify(func(k:k),'EeDd')==0 .and. verify(func(k+1:k+1),'+-')==0)) &
               .or. func(k+1:k+1) == ' ') exit
          k = k + 1
        end do
        k = k + 1
        ntokens = ntokens + 1

        !It's an operator or delimitator
      else
        k = k + 1
        ntokens = ntokens + 1
      end if
    end do

    allocate(stokens(ntokens))

    k = 1
    i = 1

    do while (k <= len_trim(func))
      !It's a variable, or function  name
      if (index(chars,func(k:k)) /= 0) then
        k0 = k
        do while (k < lenf)
          if (index(operators, func(k+1:k+1)) /= 0 .or. func(k+1:k+1) == ' ') exit
          k = k + 1
        end do
        stokens(i) = func(k0:k)
        k = k + 1
        i = i + 1

        !It's a number
      else if (index(numbers,func(k:k)) /= 0) then
        k0 = k
        do while (k < lenf)
          if ( (index(operators, func(k+1:k+1)) /= 0 &
                 .and..not.(verify(func(k:k),'EeDd')==0 .and. verify(func(k+1:k+1),'+-')==0)) &
               .or. func(k+1:k+1) == ' ') exit
          k = k + 1
        end do
        stokens(i) = func(k0:k)
        k = k + 1
        i = i + 1

        !It's an operator or delimitator
      else
        stokens(i) = func(k:k)
        if(stokens(i) == '(')then
          irightbrackets = irightbrackets + 1
        else if(stokens(i) == ')') then
          ileftbrackets = ileftbrackets + 1
        end if
        i = i + 1
        k = k + 1
      end if
    end do

    if (irightbrackets /= ileftbrackets) then
      status_fparse = -1
      return
    end if

    itoke = 1
    isaddsub = 1
    ismuldiv = 1
    ioperations = 1
    numberk = 1
    toke = stokens(itoke)
    allocate(opaddsub(2))
    allocate(opmuldiv(2))
    allocate(number(ntokens))
    allocate(pdata(ntokens))
    allocate(operations(ntokens))

    call add_sub()
    ioperations = ioperations - 1

  end subroutine tokens_analyzer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !The following subroutines call themselves recursively
  !to build the expression to be parsed based on an algorithm
  !called Recursive Descendent Parsing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_sub ()
    !
    ! Enter description here
    !

    call mul_div ()

    do while (trim(toke) == '+' .or. trim(toke) == '-')
      opaddsub(isaddsub) = trim(toke)
      isaddsub = isaddsub + 1
      itoke = itoke + 1
      toke = stokens(itoke)
      call mul_div()

      selectcase(opaddsub(isaddsub-1))
      case('+')
        isaddsub = isaddsub - 1
        operations(ioperations) = 3
        ioperations = ioperations + 1

      case('-')
        isaddsub = isaddsub - 1
        operations(ioperations) = 4
        ioperations = ioperations + 1
      end select
    end do

  end subroutine add_sub

  subroutine mul_div ()
    !
    ! Enter description here
    !

    call unary()

    do while (trim(toke) == '*' .or. trim(toke) == '/')
      opmuldiv(ismuldiv) = trim(toke)
      ismuldiv = ismuldiv + 1
      itoke = itoke + 1
      toke = stokens(itoke)
      call unary()

      selectcase(opmuldiv(ismuldiv-1))
      case('*')
        ismuldiv = ismuldiv - 1
        operations(ioperations) = 5
        ioperations = ioperations + 1
      case('/')
        ismuldiv = ismuldiv - 1
        operations(ioperations) = 6
        ioperations = ioperations + 1
      end select
    end do

  end subroutine mul_div

  subroutine unary()
    !
    ! Enter description here
    !

    if (trim(toke) == '-') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call pow()
      operations(ioperations) = 2
      ioperations = ioperations + 1
    else if (trim(toke) == '+') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call pow()
    else
      call pow()
    end if

  end subroutine unary

  subroutine pow ()
    !
    ! Enter description here
    !

    call functions()

    if (trim(toke) == '^') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call functions()
      operations(ioperations) = 7
      ioperations = ioperations + 1
    end if

  end subroutine pow

  subroutine functions ()
    !
    ! Enter description here
    !
    if (trim(toke) == 'sin') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 8
      ioperations = ioperations + 1

    else if(trim(toke) == 'cos') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 9
      ioperations = ioperations + 1

    else if(trim(toke) == 'tan') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 10
      ioperations = ioperations + 1

    else if(trim(toke) == 'asin') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 11
      ioperations = ioperations + 1

    else if(trim(toke) == 'acos') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 12
      ioperations = ioperations + 1

    else if(trim(toke) == 'atan') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 13
      ioperations = ioperations + 1

    else if(trim(toke) == 'sinh') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 14
      ioperations = ioperations + 1

    else if(trim(toke) == 'cosh') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 15
      ioperations = ioperations + 1

    else if(trim(toke) == 'tanh') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 16
      ioperations = ioperations + 1

    else if(trim(toke) == 'sind') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 17
      ioperations = ioperations + 1

    else if(trim(toke) == 'cosd') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 18
      ioperations = ioperations + 1

    else if(trim(toke) == 'tand') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 19
      ioperations = ioperations + 1

    else if (trim(toke) == 'log') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 20
      ioperations = ioperations + 1

    else if (trim(toke) == 'log10') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 21
      ioperations = ioperations + 1

    else if (trim(toke) == 'nint') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 22
      ioperations = ioperations + 1

    else if (trim(toke) == 'anint') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 23
      ioperations = ioperations + 1

    else if (trim(toke) == 'aint') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 24
      ioperations = ioperations + 1

    else if (trim(toke) == 'exp') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 25
      ioperations = ioperations + 1

    else if (trim(toke) == 'sqrt') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 26
      ioperations = ioperations + 1

    else if (trim(toke) == 'abs') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 27
      ioperations = ioperations + 1

    else if (trim(toke) == 'floor') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call brackets()
      operations(ioperations) = 28
      ioperations = ioperations + 1

    else
      call brackets()

    end if

  end subroutine functions

  subroutine brackets()
    !
    ! Enter description here
    !
    if (trim(toke) == '(') then
      itoke = itoke + 1
      toke = stokens(itoke)
      call add_sub()
      if (trim(toke) /= ')') then
        status_fparse =  -2
        return
      end if
      if (itoke < ntokens) then
        itoke = itoke + 1
        toke = stokens(itoke)
      end if
      if (trim(toke) == '(') then
        status_fparse =  -3
        return
      end if
    else
      call recog_vars ()
    end if

  end subroutine brackets

  subroutine recog_vars ()
    !
    ! Enter description here
    !

    integer          :: i, j, i_used
    integer          :: ierror

    !Expression has an error
    if (index(operators, trim(toke)) /= 0) then
      err_msg_add = ': "'//trim(toke)//'"'
      status_fparse = -6
      return
    end if

    do i = 1, nvar
      !It's a variable
!      if (trim(toke) == varnames(i)) then
      if (match_var(varnames(i), trim(toke))) then
        ! print*,trim(toke)//' is variable '//trim(varnames(i))
        if (used_vnames_) then
          i_used = -1
          do j = 1, nvar_used
            if (trim(toke) == trim(varnames_used(j))) then
              i_used = j
              exit
            end if
          end do
          if (i_used <= 0) then
            nvar_used = nvar_used + 1
            if (size(varnames_used) < nvar_used) then
              allocate(varnames_used_(nvar_used + 5))
              varnames_used_(1:nvar_used-1) = varnames_used(1:nvar_used-1)
              deallocate(varnames_used)
              varnames_used => varnames_used_
              varnames_used_ => null()
            end if
            varnames_used(nvar_used) = trim(toke)
            i_used = nvar_used
          end if
          !print*,'new used variable',nvar_used,trim(toke)
          operations(ioperations) = 28+i_used
        else
          operations(ioperations) = 28+i
        end if
        ioperations = ioperations + 1
        if (itoke < ntokens) then
          itoke = itoke + 1
          toke = stokens(itoke)
        end if
        return
      end if
    end do
    !print*,trim(toke)//' is not a variable'

    !It's a number
    toke = trim(toke)
    read(toke, *, iostat = ierror) number(numberk)
    if (ierror /= 0) then
      err_msg_add = ': "'//trim(toke)//'"'
      status_fparse = -7
      return
    else
      operations(ioperations) = 1
      ioperations = ioperations + 1
      if (itoke < ntokens) then
        itoke = itoke + 1
        toke = stokens(itoke)
      end if
      numberk = numberk + 1
    end if

  end subroutine recog_vars


  function evaluate_fparse (vars, tfp) result (answer)
    !
    !       This function will evaluate the expression supplied
    !

    real(kind = wp)                                :: answer
    real(kind = wp), intent(in)                    :: vars(:)
    type(t_fparse),  intent(in),  optional, target :: tfp

    integer                      :: st
    integer                      :: dt
    integer                      :: i

    if (present(tfp)) then
      tmp_ => tfp
    else
      tmp_ => tmp
    end if

    answer = invalid
    if (tmp_%ioperations <= 0) RETURN
    if (size(vars) /= tmp_%nvar) RETURN

    st = 0
    dt = 1
    do i = 1, tmp_%ioperations
      select case(tmp_%operations(i))
      case (1)
        st = st + 1
        tmp_%pdata(st) = tmp_%number(dt)
        dt = dt + 1
      case (2)
        tmp_%pdata(st) =     - tmp_%pdata(st)
      case (3)
        tmp_%pdata(st-1) = tmp_%pdata(st-1) + tmp_%pdata(st)
        st = st - 1
      case (4)
        tmp_%pdata(st-1) = tmp_%pdata(st-1) - tmp_%pdata(st)
        st = st - 1
      case (5)
        tmp_%pdata(st-1) = tmp_%pdata(st-1) * tmp_%pdata(st)
        st = st - 1
      case (6)
        tmp_%pdata(st-1) = tmp_%pdata(st-1) / tmp_%pdata(st)
        st = st - 1
      case (7)
        tmp_%pdata(st-1) = tmp_%pdata(st-1) ** tmp_%pdata(st)
        st = st - 1
      case (8)
        tmp_%pdata(st) = sin(tmp_%pdata(st))
      case (9)
        tmp_%pdata(st) = cos(tmp_%pdata(st))
      case (10)
        tmp_%pdata(st) = tan(tmp_%pdata(st))
      case (11)
        tmp_%pdata(st) = asin(tmp_%pdata(st))
      case (12)
        tmp_%pdata(st) = acos(tmp_%pdata(st))
      case (13)
        tmp_%pdata(st) = atan(tmp_%pdata(st))
      case (14)
        tmp_%pdata(st) = sinh(tmp_%pdata(st))
      case (15)
        tmp_%pdata(st) = cosh(tmp_%pdata(st))
      case (16)
        tmp_%pdata(st) = tanh(tmp_%pdata(st))
      case (17)
        tmp_%pdata(st) = sin(tmp_%pdata(st)*deg2rad)
      case (18)
        tmp_%pdata(st) = cos(tmp_%pdata(st)*deg2rad)
      case (19)
        tmp_%pdata(st) = tan(tmp_%pdata(st)*deg2rad)
      case (20)
        tmp_%pdata(st) = log(tmp_%pdata(st))
      case (21)
        tmp_%pdata(st) = log10(tmp_%pdata(st))
      case (22)
        tmp_%pdata(st) = nint(tmp_%pdata(st))
      case (23)
        tmp_%pdata(st) = anint(tmp_%pdata(st))
      case (24)
        tmp_%pdata(st) = aint(tmp_%pdata(st))
      case (25)
        tmp_%pdata(st) = exp(tmp_%pdata(st))
      case (26)
        tmp_%pdata(st) = sqrt(tmp_%pdata(st))
      case (27)
        tmp_%pdata(st) = abs(tmp_%pdata(st))
      case (28)
        tmp_%pdata(st) = floor(tmp_%pdata(st))
      case default
        st = st + 1
        tmp_%pdata(st) = vars(tmp_%operations(i)-28)
      end select
    end do

    answer = tmp_%pdata(1)

  end function evaluate_fparse


  function evaluate_save (vars, tfp) result (answer)
    !
    !       This function will evaluate the expression supplied
    !
    real(kind = wp)                                :: answer
    real(kind = wp), intent(in)                    :: vars(:)
    type(t_fparse),  intent(in),  optional, target :: tfp

    integer                      :: st
    integer                      :: dt
    integer                      :: i

    if (present(tfp)) then
      tmp_ => tfp
    else
      tmp_ => tmp
    end if

    answer = invalid
    if (tmp_%ioperations <= 0) RETURN
    if (size(vars) /= tmp_%nvar) RETURN

    st = 0
    dt = 1

    do i = 1, tmp_%ioperations
      select case(tmp_%operations(i))
      case (1)
        st = st + 1
        tmp_%pdata(st) = tmp_%number(dt)
        dt = dt + 1
      case (2)
        tmp_%pdata(st) =     - tmp_%pdata(st)
      case (3)
        tmp_%pdata(st-1) = tmp_%pdata(st-1) + tmp_%pdata(st)
        st = st - 1
      case (4)
        tmp_%pdata(st-1) = tmp_%pdata(st-1) - tmp_%pdata(st)
        st = st - 1
      case (5)
        tmp_%pdata(st-1) = tmp_%pdata(st-1) * tmp_%pdata(st)
        st = st - 1
      case (6)
        if(abs(tmp_%pdata(st)) < 1.0e-30) then
          answer = invalid
          return
        end if
        tmp_%pdata(st-1) = tmp_%pdata(st-1) / tmp_%pdata(st)
        st = st - 1
      case (7)
        if(tmp_%pdata(st-1) < 0.0 .AND. (tmp_%pdata(st)-int(tmp_%pdata(st))) /= 0.0) then
          answer = invalid
          return
        end if
        if(tmp_%pdata(st)*log(abs(tmp_%pdata(st-1)+1.0E-15)) > 65.0) then
          answer = invalid
          return
        end if
        tmp_%pdata(st-1) = tmp_%pdata(st-1) ** tmp_%pdata(st)
        st = st - 1
      case (8)
        tmp_%pdata(st) = sin(tmp_%pdata(st))
      case (9)
        tmp_%pdata(st) = cos(tmp_%pdata(st))
      case (10)
        if((abs(tmp_%pdata(st)) > 89.99*3.141593/180. .and. abs(tmp_%pdata(st)) < 90.01*3.141593/180)&
             .or. (abs(tmp_%pdata(st)) > 269.99*3.141593/180. .and. abs(tmp_%pdata(st)) < 270.01*3.141593/180)) then
          answer = invalid
          return
        end if
        tmp_%pdata(st) = tan(tmp_%pdata(st))
      case (11)
        if(abs(tmp_%pdata(st)) > 1.0) then
          answer = invalid
          return
        end if
        tmp_%pdata(st) = asin(tmp_%pdata(st))
      case (12)
        if(abs(tmp_%pdata(st)) > 1.0) then
          answer = invalid
          return
        end if
        tmp_%pdata(st) = acos(tmp_%pdata(st))
      case (13)
        if(abs(tmp_%pdata(st)) > 1.0e+10) then
          answer = invalid
          return
        end if
        tmp_%pdata(st) = atan(tmp_%pdata(st))
      case (14)
        if(tmp_%pdata(st) > 60) then
          answer = invalid
          return
        end if
        tmp_%pdata(st) = sinh(tmp_%pdata(st))
      case (15)
        if(tmp_%pdata(st) > 60) then
          answer = invalid
          return
        end if
        tmp_%pdata(st) = cosh(tmp_%pdata(st))
      case (16)
        tmp_%pdata(st) = tanh(tmp_%pdata(st))
      case (17)
        tmp_%pdata(st) = sin(tmp_%pdata(st)*deg2rad)
      case (18)
        tmp_%pdata(st) = cos(tmp_%pdata(st)*deg2rad)
      case (19)
        tmp_%pdata(st) = tan(tmp_%pdata(st)*deg2rad)
      case (20)
        if(tmp_%pdata(st) <= 1.0e-15) then
          answer = invalid
          return
        end if
        tmp_%pdata(st) = log(tmp_%pdata(st))
      case (21)
        if(tmp_%pdata(st) <= 1.0e-15) then
          answer = invalid
          return
        end if
        tmp_%pdata(st) = log10(tmp_%pdata(st))
      case (22)
        tmp_%pdata(st) = nint(tmp_%pdata(st))
      case (23)
        tmp_%pdata(st) = anint(tmp_%pdata(st))
      case (24)
        tmp_%pdata(st) = aint(tmp_%pdata(st))
      case (25)
        if(tmp_%pdata(st) > 55) then
          answer = invalid
          return
        end if
        tmp_%pdata(st) = exp(tmp_%pdata(st))
      case (26)
        if(tmp_%pdata(st) < 0) then
          answer = invalid
          return
        end if
        tmp_%pdata(st) = sqrt(tmp_%pdata(st))
      case (27)
        tmp_%pdata(st) = abs(tmp_%pdata(st))
      case (28)
        tmp_%pdata(st) = floor(tmp_%pdata(st))
      case default
        st = st + 1
        tmp_%pdata(st) = vars(tmp_%operations(i)-28)
      end select
      if(abs(tmp_%pdata(st)) > 1.0d+60) then
        answer = invalid
        return
      end if

    end do

    answer = tmp_%pdata(1)

  end function evaluate_save

  elemental subroutine destruct_fparse1(tfp)
    type(t_fparse),  intent(inout) :: tfp
    if (associated(tfp%varnames  )) deallocate(tfp%varnames  )
    if (associated(tfp%operations)) deallocate(tfp%operations)
    if (associated(tfp%number    )) deallocate(tfp%number    )
    if (associated(tfp%pdata     )) deallocate(tfp%pdata     )
    call construct(tfp)
  end subroutine destruct_fparse1

  subroutine destruct_fparse2()
    if (associated(tmp%varnames  )) deallocate(tmp%varnames  )
    if (associated(tmp%operations)) deallocate(tmp%operations)
    if (associated(tmp%number    )) deallocate(tmp%number    )
    if (associated(tmp%pdata     )) deallocate(tmp%pdata     )
    call construct(tmp)
  end subroutine destruct_fparse2

  elemental subroutine construct_fparse(tmp)
    type(t_fparse),  intent(out) :: tmp
  end subroutine construct_fparse

  !> This subroutine removes unnecessary blank spaces
  recursive subroutine blanks(func)
    character(len=*), intent(inout) :: func
    integer                         :: k

    func = adjustl(func)
    k = index(trim(func), ' ')
    if (k /= 0) then
      func = func(:k-1) // func(k+1:)
      call blanks(func)
    end if

  end subroutine blanks

  !> Checks for errors
  function check_err(func) result(status)
    integer                          :: status
    character (len=*), intent(inout) :: func

    character(len=5),      parameter :: funcs5( 3) = (/'log10','anint','floor'/)
    integer,               parameter :: stat5 ( 3) = (/     9 ,    10 ,    11 /)
    character(len=4),      parameter :: funcs4(12) = (/'asin','acos','atan','sinh','cosh','tanh',&
                                                       'sind','cosd','tand','nint','aint','sqrt'/)
    integer,               parameter :: stat4 (12) = (/   12 ,   13 ,   14 ,   15 ,   16 ,   17 ,&
                                                          18 ,   19 ,   20 ,   21 ,   22 ,   23 /)
    character(len=5),      parameter :: funcs3 (6) = (/'sin','cos','tan','log','exp','abs'/)
    integer,               parameter :: stat3  (6) = (/  24 ,  25 ,  26 ,  27 ,  28 ,  29 /)
    
    integer,               parameter :: n_variav = 37
    character(n_variav),   parameter :: variav = '0123456789abcdefghijklmnopqrstuvwxyz_'
    integer                          :: i, j, nchar
    logical                          :: l_okay

    
    
    status = 0

    nchar = len_trim(func)

    if(func(nchar:nchar) == '-' .or. func(nchar:nchar) == '+' .or. func(nchar:nchar) == '/' .or. func(nchar:nchar) == '*') then
      status=1
      return
    end if

    if(func(1:1) == '*' .or. func(1:1) == '/') then
      status=2
      return
    end if

    do i = 1, nchar-1
      if(func(i:i+1) == '--' .or. func(i:i+1) == '-+' .or. func(i:i+1) == '-/' .or. func(i:i+1) == '-*') status=3
      if(func(i:i+1) == '+-' .or. func(i:i+1) == '++' .or. func(i:i+1) == '+/' .or. func(i:i+1) == '+*') status=4
      if(func(i:i+1) == '*-' .or. func(i:i+1) == '*+' .or. func(i:i+1) == '*/')                          status=5
      if(func(i:i+1) == '/-' .or. func(i:i+1) == '/+' .or. func(i:i+1) == '//' .or. func(i:i+1) == '/*') status=6
    end do
    if(status /= 0) return

    do i = 1, nchar-1
      do j = 1, n_variav
        if(func(i:i+1) == ')'//variav(j:j)) status=7
      end do
    end do
    if(status /= 0) return

    do i = 1, nchar-1
      do j = 1, n_variav
        if (func(i:i) == '0' .or. func(i:i) == 'n' .or. func(i:i) == 's' .or. &
            func(i:i) == 'h' .or. func(i:i) == 'd' .or. func(i:i) == 'g' .or. &
            func(i:i) == 't' .or. func(i:i) == 'p' .or. func(i:i) == 'r') then
          ! do not test, can be one of the defined functions
        else
          if(func(i:i+1) == variav(j:j)//'(') status=8
        end if
      end do
    end do
    if(status /= 0) return

    if(nchar >= 5) then
      do i = 1,nchar-4
        do j = 1, size(funcs5)
          if(func(i:i+4) == funcs5(j)) then
            if (.not.follows_valid(func(i+5:))) then ; status=stat5(j) ; return ; endif
          end if
        end do
      end do
    end if
    if(nchar >= 4) then
      do i = 1,nchar-3
        do j = 1, size(funcs4)
          if(func(i:i+3) == funcs4(j)) then
            if (.not.follows_valid(func(i+4:))) then ; status=stat4(j) ; return ; endif
          end if
        end do
      end do
    end if
    if(nchar >= 3) then
      do i = 1,nchar-2
        do j = 1, size(funcs3)
          if(func(i:i+2) == funcs3(j)) then
            if (.not.follows_valid(func(i+3:))) then ; status=stat3(j) ; return ; endif
          end if
        end do
      end do
    end if

  end function check_err


  function err_msg(status) result(msg)
    character(len=100)  :: msg
    integer, intent(in) :: status

    select case(status)
    case(-7)
      msg = 'Invalid number'
    case(-6)
      msg = 'Expected a variable, but got an operator'
    case(-5)
      msg = 'Expanded function too long. Recompile with increased "lt" in mo_fparser.'
    case(-4)
      msg = 'Function too long. Recompile with increased "lt" in mo_fparser.'
    case(-3)
      msg = 'Invalid left bracket'
    case(-2)
      msg = 'Missing right bracket'
    case(-1)
      msg = 'Number of left and right brackets not equal'
    case( 1)
      msg = 'Trailing operator at the end of function'
    case( 2)
      msg = 'Function starts with operator "*" or "/"'
    case( 3)
      msg = 'Operator follows "-"'
    case( 4)
      msg = 'Operator follows "+"'
    case( 5)
      msg = 'Operator follows "*"'
    case( 6)
      msg = 'Operator follows "/"'
    case( 7)
      msg = 'Missing operator after ")"'
    case( 8)
      msg = 'No operator or function in front of "("'
    case( 9)
      msg = 'No bracket after function "log10"'
    case(10)
      msg = 'No bracket after function "anint"'
    case(11)
      msg = 'No bracket after function "floor"'
    case(12)
      msg = 'No bracket after function "asin"'
    case(13)
      msg = 'No bracket after function "acos"'
    case(14)
      msg = 'No bracket after function "atan"'
    case(15)
      msg = 'No bracket after function "sinh"'
    case(16)
      msg = 'No bracket after function "cosh"'
    case(17)
      msg = 'No bracket after function "tanh"'
    case(18)
      msg = 'No bracket after function "sind"'
    case(19)
      msg = 'No bracket after function "cosd"'
    case(20)
      msg = 'No bracket after function "tand"'
    case(21)
      msg = 'No bracket after function "nint"'
    case(22)
      msg = 'No bracket after function "aint"'
    case(23)
      msg = 'No bracket after function "sqrt"'
    case(24)
      msg = 'No bracket or "h" or "d" after function "sin"'
    case(25)
      msg = 'No bracket or "h" or "d" after function "cos"'
    case(26)
      msg = 'No bracket or "h" or "d" after function "tan"'
    case(27)
      msg = 'No bracket after function "log"'
    case(28)
      msg = 'No bracket after function "exp"'
    case(29)
      msg = 'No bracket after function "abs"'
    case default
      msg = 'Unknown error'
    end select

    msg = trim(msg)//trim(err_msg_add)

  end function err_msg


  subroutine convert(func)
    character(len=lt), intent(inout) :: func
    integer :: k, lenf

    k = 1
    lenf = len_trim(func)

    do while (k <= lenf)
      if (func(k:k) == '[' .or. func(k:k) == '{') then
        func(k:k)='('
      elseif (func(k:k) == ']' .or. func(k:k) == '}') then
        func(k:k)=')'
      elseif (func(k:k) == ',') then
        func(k:k)='.'
      elseif (k < lenf) then
        select case(func(k:k+1))
        case('**')
          func = func(:k-1)//'^'//trim(func(k+2:))
          lenf = lenf - 1
        case('pi')
          if (is_word(func,k,k+1)) then
            func = func(:k-1)//'3.14159265359'//trim(func(k+2:))
            lenf = lenf + 11
            k = k + 12
          end if
        case('ln')
          if (is_word(func,k,k+1)) then
            func = func(:k-1)//'log'//trim(func(k+2:))
            lenf = lenf + 1
            k = k + 2
          end if
        end select
      end if
      k = k + 1
    end do

  end subroutine convert


  function follows_bracket(str, nchar, ind, set) result(ii)
    integer :: ii
    character(len=*), intent(in)           :: str
    integer,          intent(in)           :: nchar
    integer,          intent(in)           :: ind
    character(len=*), intent(in), optional :: set
    integer            :: j
    character(len=lt)  :: set_

    ii = -1
    if (present(set)) then
      set_ = set
    else
      set_ = '('
    end if
    do j = ind,nchar
      if(str(j:j) /= ' ') then
        if (scan(str(j:j),trim(set_)) > 0) then
          ii = j
          exit
        else
          return
        end if
      end if
    end do
  end function follows_bracket

  function follows_valid(str) result(lvalid)
    logical                      :: lvalid
    character(len=*), intent(in) :: str
    character(len=*), parameter  :: set = '(0123456789abcdefghijklmnopqrstuvwxyz_'
    integer :: j
    
    lvalid = .false.
    do j = 1, len_trim(str)
      if(str(j:j) /= ' ') then
        if (scan(str(j:j),trim(set)) > 0) lvalid = .true.
        return
      end if
    end do
    
  end function follows_valid

  function match_var(vname ,str) result(l)
    logical :: l
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: str

    character(len=lv) :: v(10)
    integer           :: l_s, l_v, i_s
    integer           :: i, n, stat
    logical           :: l_prev_star

    l = .false.
    call split2(vname, sep='*', array=v, n=n, status=stat)
    if (stat /= 0) RETURN

    if (n > 0) then
      l_prev_star = (vname(1:1) == '*')
      l_s = len_trim(str)
      i_s = 1
      do i = 1, n
        if (i_s > l_s) RETURN ! we are at the end of str
        l_v = len_trim(v(i))
        if (i_s+l_v-1 > l_s) RETURN
        if (l_prev_star) then
          do while (str(i_s:i_s+l_v-1) /= trim(v(i)) .and. i_s+l_v-1 < l_s)
            i_s = i_s + 1
          end do
        end if
        if (str(i_s:i_s+l_v-1) /= trim(v(i))) RETURN
        i_s = i_s + l_v
        l_prev_star = .true.
      end do
      l_v = len_trim(vname)
      if (vname(l_v:l_v) /= '*') then
        ! Check for trailing stuff in str
        if (i_s <= l_s) RETURN ! There is still trailing stuff at the end of str
      end if
      l = .true.
    else
      l = .true.
    end if

  end function match_var


  function is_word(text,i_s,i_e) result(l)
    logical :: l
    character(len=*), intent(in) :: text
    integer,          intent(in) :: i_s
    integer,          intent(in) :: i_e

    l = .true.

    if (i_s > 1) then
      l = l .and. (verify(text(i_s-1:i_s-1), 'abcdefghijklmnopqrstuvwxyz0123456789_') > 0)
    end if
    if (i_e < len_trim(text)) then
      l = l .and. (verify(text(i_e+1:i_e+1), 'abcdefghijklmnopqrstuvwxyz0123456789_') > 0)
    end if

  end function is_word


  subroutine bcast_fparse (fp, source, comm)
    type (t_fparse)   ,intent(inout) :: fp
    integer           ,intent(in)    :: source
    integer ,optional ,intent(in)    :: comm

    integer :: lcom
    integer :: count = 0
    integer :: n

    lcom = dace% comm ;if(present(comm)) lcom = comm

    if (count == 0) count = size(transfer(fp,(/' '/)))

    if (dace% pe /= source) call destruct (fp)
    call p_bcast_derivedtype (fp, count, source, lcom)
    if (dace% pe == source) then
      n = 0
      if (associated(fp%pdata)) n = size(fp%pdata)
    end if
    call p_bcast (n,  source, lcom)

    if (dace% pe /= source) then
      if (associated(fp% varnames  )) allocate(fp% varnames  (fp%nvar       ))
      if (associated(fp% operations)) allocate(fp% operations(fp%ioperations))
      if (associated(fp% number    )) allocate(fp% number    (fp%numberk    ))
      if (associated(fp% pdata     )) allocate(fp% pdata     (n             ))
    end if
    if (associated(fp% varnames  )) call p_bcast(fp% varnames  (1:fp%nvar       ), source, lcom)
    if (associated(fp% operations)) call p_bcast(fp% operations(1:fp%ioperations), source, lcom)
    if (associated(fp% number    )) call p_bcast(fp% number    (1:fp%numberk    ), source, lcom)
    if (associated(fp% pdata     )) call p_bcast(fp% pdata     (1:n             ), source, lcom)

  end subroutine bcast_fparse

  elemental subroutine assign_fparse(y, x)
    type(t_fparse), intent(inout) :: y
    type(t_fparse), intent(in)    :: x

    call destruct(y)
    y%func        = x%func
    y%nvar        = x%nvar
    y%ioperations = x%ioperations
    y%numberk     = x%numberk
    if (associated(x%varnames)) then
      allocate(y%varnames(size(x%varnames)))
      y%varnames = x%varnames
    end if
    if (associated(x%operations)) then
      allocate(y%operations(size(x%operations)))
      y%operations = x%operations
    end if
    if (associated(x%number)) then
      allocate(y%number(size(x%number)))
      y%number = x%number
    end if
    if (associated(x%pdata)) then
      allocate(y%pdata(size(x%pdata)))
      y%pdata = x%pdata
    end if

  end subroutine assign_fparse

end module mo_fparse
