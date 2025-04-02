!
!+ derived type definition to to store domain decomposition information
!
MODULE mo_atm_decomp
!
! Description:
!   Definition of type 't_atm_dec' to store domain decomposition information.
!   To be used as component of 't_grid' (grid definition for
!   atmospheric states).
!   Definition of subroutines to set up the domain decomposition.
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
! V1_7         2009/08/24 Andreas Rhodin
!  make subroutine 'print_decomp' public, change output format
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  remove unused variables
! V1_22        2013-02-13 Harald Anlauf
!  initial support for ICON unstructured grid
! V1_26        2013/06/27 Andreas Rhodin
!  construct_grid: revised handling of decomposition info
! V1_28        2014/02/26 Andreas Rhodin
!  fixes for the case that nproc1/2 exceeds the respective # of grid-points
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Luis Kornblueh  MPI  2002       routines gathered from GME code
! Andreas Rhodin  DWD  2002-2005
!==========================================================================

  !-------------
  ! Modules used
  !-------------
  use mo_mpi_dace, only:  dace,            &! DACE communication info
                          self,            &! self communication info
                          p_max,           &! overloaded mpi_max routine
                          p_barrier,       &! encapsulated MPI_barrier routine
!                         MPI_STATUS_SIZE, &! size of status return variable
                          MPI_INTEGER       ! MPI integer kind parameter
  use mo_exception, only: finish
  IMPLICIT NONE

  !----------------
  ! Public entities
  !----------------
  private
  public :: lmpi
  public :: t_atm_dec      ! domain decomposition data type
  public :: setup_single   ! setup single processor run
  public :: setup_parallel ! setup multi  processor run
  public :: setup_decomp   ! setup decomposition
  public :: setup_com_tri  ! setup communication for icosahedral grid
  public :: setup_com_reg  ! setup communication for regular grid
  public :: setup_com_icon ! setup communication for ICON unstructured grid
  public :: destruct       ! deallocate pointer components of t_atm_dec
  public :: print_decomp   ! print decomposition

  LOGICAL, PARAMETER :: lmpi = .TRUE. ! Wether MPI is used or not

  !---------------------
  ! Data type definition
  !---------------------
  type t_atm_dec
    !---------------------------------------------------------------
    ! parameters set by subroutines 'setup_single', 'setup_parallel'
    !---------------------------------------------------------------
    integer :: comm    ! communicator group for this model instance
    integer :: nproc1  ! Number of PEs in direction 1
    integer :: nproc2  ! Number of PEs in direction 2
    !---------------------------------------------
    ! parameters set up by subroutine 'decomp_tri'
    !---------------------------------------------
    integer :: nproc   ! Number of all processors MPI:=1 ??
    integer :: npe     ! nproc1 * nproc2
    integer :: myproc  ! Number (Rank) of own processor element
    integer :: my_num1 ! Number in direction 1
    integer :: my_num2 ! Number in direction 2

    integer :: nbpe    ! Number of neighbor poleward east
    integer :: nbaw    ! Number of neighbor antipoleward west
    integer :: nbpw    ! Number of neighbor poleward west
    integer :: nbae    ! Number of neighbor antipoleward east

    integer, pointer :: ilim1 (:)   => NULL() ! Limits of dec. in direction 1
    integer, pointer :: ilim2 (:)   => NULL() ! Limits of dec. in direction 2
    integer, pointer :: ilimp (:)   => NULL() ! Limits of dec. for paral.vector
    integer, pointer :: pe_12 (:,:) => NULL() ! Rank of PEs on 2-dim grid
    !--------------------------------------------------
    ! for a better communication in the post-processing
    !--------------------------------------------------
    integer, pointer :: isubpositions (:,:) => NULL() ! limits of every subdom.
    integer          :: imaxnumgp                     ! max.# of gp. in subdom.
    !----------------------------------------------------
    ! parameters set by subroutines 'setup_com'
    !----------------------------------------------------
    integer, pointer :: idx_recv  (:,:) => NULL() ! Indices to put recvd points
    integer, pointer :: idx_send  (:,:) => NULL() ! Indices to get pnts to send
    integer, pointer :: np_recv_s (:)   => NULL() ! start index in array idx_recv
    integer, pointer :: np_recv_1 (:)   => NULL() ! #of points to recv - 1 line
    integer, pointer :: np_recv_2 (:)   => NULL() ! dto - 2 boundary lines
    integer, pointer :: np_send_s (:)   => NULL() ! start index in idx_send
    integer, pointer :: np_send_1 (:)   => NULL() ! #of points to send - 1 line
    integer, pointer :: np_send_2 (:)   => NULL() ! dto - 2 boundary lines
    integer, pointer :: np_recv_t (:)   => NULL() ! startin idx_recv -target PE
                                                  !  (for shmem-communication)
    integer          :: np_recv_tot               ! tot. # of points to receive
    integer          :: np_send_tot               ! tot. # of points to send
    integer          :: np_recv_max               ! max. # of points to receive
    integer          :: np_send_max               ! max. # of points to send
  end type t_atm_dec

  !-----------
  ! Interfaces
  !-----------

  interface destruct
    module procedure destruct_t_atm_dec
  end interface destruct

  ! Neighborhood relationships within the diamonds
  ! mpw, mpe, maw, mae give the numbers of the poleward or antipoleward
  ! west or east neighbors of a diamond

  INTEGER, PARAMETER :: mpw(10) = (/ 5, 1, 2, 3, 4,10, 6, 7, 8, 9 /)
  INTEGER, PARAMETER :: mpe(10) = (/ 2, 3, 4, 5, 1, 7, 8, 9,10, 6 /)
  INTEGER, PARAMETER :: maw(10) = (/10, 6, 7, 8, 9, 1, 2, 3, 4, 5 /)
  INTEGER, PARAMETER :: mae(10) = (/ 6, 7, 8, 9,10, 2, 3, 4, 5, 1 /)

  INTEGER :: mpi_err ! just for convenience

CONTAINS
!------------------------------------------------------------------------------
  subroutine destruct_t_atm_dec (dc)
  type (t_atm_dec) ,intent(inout) :: dc

    if (associated (dc% ilim1        )) deallocate (dc% ilim1        )
    if (associated (dc% ilim2        )) deallocate (dc% ilim2        )
    if (associated (dc% ilimp        )) deallocate (dc% ilimp        )
    if (associated (dc% pe_12        )) deallocate (dc% pe_12        )
    if (associated (dc% isubpositions)) deallocate (dc% isubpositions)
    if (associated (dc% idx_recv     )) deallocate (dc% idx_recv     )
    if (associated (dc% idx_send     )) deallocate (dc% idx_send     )
    if (associated (dc% np_recv_s    )) deallocate (dc% np_recv_s    )
    if (associated (dc% np_recv_1    )) deallocate (dc% np_recv_1    )
    if (associated (dc% np_recv_2    )) deallocate (dc% np_recv_2    )
    if (associated (dc% np_send_s    )) deallocate (dc% np_send_s    )
    if (associated (dc% np_send_1    )) deallocate (dc% np_send_1    )
    if (associated (dc% np_send_2    )) deallocate (dc% np_send_2    )
    if (associated (dc% np_recv_t    )) deallocate (dc% np_recv_t    )

  end subroutine destruct_t_atm_dec
!------------------------------------------------------------------------------
  subroutine print_decomp (dc, unit)
  type (t_atm_dec) ,intent(in)           :: dc
  integer          ,intent(in) ,optional :: unit

    integer :: iu
    character(len=*), parameter :: fint10 = '(a, 8i9/(17x, 8i9))'
    character(len=*), parameter :: fint5  = '(a, 5i9/(17x, 5i9))'
!   character(len=*), parameter :: fint3  = &
!     '(a,3i5,1x,3i5,1x,3i5/(17x, 3i5,1x,3i5,1x,3i5))'

    iu = 6; if (present(unit)) iu = unit

    write(iu,'(    )')
    write(iu,'(a   )') repeat('-',79)
    write(iu,'(    )')
    write(iu,'(a   )') ' type (t_atm_dec) :'
    write(iu,'(    )')
    write(iu,'(a,i0)') ' comm          = ', dc%  comm
    write(iu,'(    )')
    write(iu,'(a,i9)') ' nproc1        = ', dc%  nproc1
    write(iu,'(a,i9)') ' nproc2        = ', dc%  nproc2
    write(iu,'(a,i9)') ' nproc         = ', dc%  nproc
    write(iu,'(a,i9)') ' my_num1       = ', dc%  my_num1
    write(iu,'(a,i9)') ' my_num2       = ', dc%  my_num2
    write(iu,'(    )')
    write(iu,'(a,i9)') ' nbpe          = ', dc%  nbpe
    write(iu,'(a,i9)') ' nbaw          = ', dc%  nbaw
    write(iu,'(a,i9)') ' nbpw          = ', dc%  nbpw
    write(iu,'(a,i9)') ' nbae          = ', dc%  nbae
    write(iu,'(    )')
    write(iu, fint10 ) ' ilim1         = ', dc%  ilim1
    write(iu, fint10 ) ' ilim2         = ', dc%  ilim2
    write(iu, fint10 ) ' ilimp         = ', dc%  ilimp
    write(iu,'(    )')
    write(iu, fint5  ) ' isubpositions = ', dc%  isubpositions
    write(iu,'(a,i9)') ' imaxnumgp     = ', dc%  imaxnumgp
!   write(iu,'(    )')
!   write(iu, fint3  ) ' idx_recv      = ', dc%  idx_recv
!   write(iu, fint3  ) ' idx_send      = ', dc%  idx_send
!   write(iu,'(    )')
!   write(iu, fint10 ) ' np_recv_s     = ', dc%  np_recv_s
!   write(iu, fint10 ) ' np_recv_1     = ', dc%  np_recv_1
!   write(iu, fint10 ) ' np_recv_2     = ', dc%  np_recv_2
!   write(iu, fint10 ) ' np_send_s     = ', dc%  np_send_s
!   write(iu, fint10 ) ' np_send_1     = ', dc%  np_send_1
!   write(iu, fint10 ) ' np_send_2     = ', dc%  np_send_2
!   write(iu, fint10 ) ' np_recv_t     = ', dc%  np_recv_t
!   write(iu,'(    )')
!   write(iu,'(a,i5)') ' np_recv_tot   = ', dc%  np_recv_tot
!   write(iu,'(a,i5)') ' np_send_tot   = ', dc%  np_send_tot
!   write(iu,'(a,i5)') ' np_recv_max   = ', dc%  np_recv_max
!   write(iu,'(a,i5)') ' np_send_max   = ', dc%  np_send_max

  end subroutine print_decomp
!------------------------------------------------------------------------------
  SUBROUTINE setup_single (dc)
  !---------------------------------------------------
  ! setup decomposition table for single processor run
  !---------------------------------------------------
  type (t_atm_dec) ,INTENT(out) :: dc ! decomposition to set

    dc% nproc   =  1
    dc% nproc1  =  1
    dc% nproc2  =  1
    dc% npe     =  1
    dc% myproc  =  0
    dc% my_num1 =  1
    dc% my_num2 =  1
    dc% comm    =  self% comm

  END SUBROUTINE setup_single
!------------------------------------------------------------------------------
  SUBROUTINE setup_parallel (dc, nproc1, nproc2, comm)
  !-------------------------------------------
  ! setup decomposition table for parallel run
  !-------------------------------------------
  type (t_atm_dec)  ,INTENT(out) :: dc     ! decomposition to set
  integer           ,INTENT(in)  :: nproc1 ! number of PEs in direction 1
  integer           ,INTENT(in)  :: nproc2 ! number of PEs in direction 2
  integer           ,INTENT(in)  :: comm   ! communicator group

    dc% nproc   = 1
    dc% nproc1  = nproc1
    dc% nproc2  = nproc2
    dc% npe     = nproc1 * nproc2
    dc% my_num1 = -1
    dc% my_num2 = -1
    dc% comm    = comm
    if      (dc% comm == self% comm ) then
      dc% myproc  = self% pe
    else if (dc% comm == dace% comm ) then
      dc% myproc  = dace% pe
    else
      call finish('setup_parallel','general case not implemented')
    endif

  END SUBROUTINE setup_parallel
!------------------------------------------------------------------------------
  SUBROUTINE setup_decomp (dc, lb, ub, ilim1, ilim2)
    !--------------------------------------------
    ! setup domain decomposition for general grid
    !--------------------------------------------
    type(t_atm_dec)  ,INTENT(inout) :: dc
    INTEGER          ,INTENT(inout) :: lb   (4)
    INTEGER          ,INTENT(inout) :: ub   (4)
    INTEGER,optional ,INTENT(in)    :: ilim1(0:)
    INTEGER,optional ,INTENT(in)    :: ilim2(0:)

    INTEGER :: ipes, irows, irem, i, j
    INTEGER :: nx, ny, ni
    !-----------------------------------------------
    ! for determining the indices of every subdomain
    !-----------------------------------------------
    INTEGER :: ilimits(5)

    nx = ub(1) - lb(1) + 1
    ny = ub(2) - lb(2) + 1
    ni = ny - 1            ! ? ni for parallel vector ?

    !---------------------------------------------------------------
    ! Allocate isubpositions: is needed in parallel and seq. version
    !---------------------------------------------------------------
    ALLOCATE (dc% isubpositions(5,0:dc% nproc1*dc% nproc2-1))

    IF (lmpi) THEN
      !======================
      ! Decomposition for MPI
      !======================
      if (dc% npe > 1) CALL p_barrier(dc% comm)
      !----------------------------
      ! Check nproc, nproc1, nproc2
      !----------------------------
      IF (dc% nproc /= 1) THEN
        IF (dc% myproc == 0) THEN
          WRITE (0,*) 'nproc = ', dc% nproc
          WRITE (0,*) 'nproc MUST be 1 for MPI-Version'
          call finish('setup_decomp','nproc MUST be 1 for MPI-Version')
        ENDIF
      ENDIF

      ipes = 1
#ifndef NOMPI
      if (dc% npe > 1) CALL MPI_Comm_size(dc% comm, ipes, mpi_err)
#endif

      IF (dc% nproc1*dc% nproc2 /= ipes) THEN
        IF (dc% myproc == 0) THEN
          WRITE (0,*) 'nproc1 = ',dc% nproc1,' dc% nproc2 = ',dc% nproc2,  &
               ' but # of MPI processors = ',ipes
          call finish('setup_decomp','invalid no. of MPI processors')
        ENDIF
      ENDIF
      !--------------------------------------
      ! Allocate and set decomposition limits
      !--------------------------------------
      ALLOCATE (dc% ilim1(0:dc% nproc1))
      ALLOCATE (dc% ilim2(0:dc% nproc2))
      ALLOCATE (dc% ilimp(0:1))
      ALLOCATE (dc% pe_12(dc% nproc1,dc% nproc2))

      irows = nx / dc% nproc1
      irem  = MOD(nx, dc% nproc1)

      dc% ilim1(0) = lb(1)
      DO i=1,irem
        dc% ilim1(i) = dc% ilim1(i-1) + (irows+1)
      ENDDO
      DO i=irem+1,dc% nproc1
        dc% ilim1(i) = dc% ilim1(i-1) + irows
      ENDDO

      irows = ny /dc% nproc2
      irem  =  MOD(ny ,dc% nproc2)

      dc% ilim2(0) = 1
      DO i=1,irem
        dc% ilim2(i) = dc% ilim2(i-1) + (irows+1)
      ENDDO
      DO i=irem+1,dc% nproc2
        dc% ilim2(i) = dc% ilim2(i-1) + irows
      ENDDO

      if (present (ilim1) .or. present (ilim2)) then
         if (.not. present (ilim1) .or. .not. present (ilim2)) &
              call finish ('setup_decomp','ilim1,ilim2 must both be present')
         if (size (ilim1) /= dc% nproc1+1 .or. size (ilim2) /= dc% nproc2+1) then
            write(*,*) "size (ilim1), nproc1 =", size (ilim1), dc% nproc1
            write(*,*) "size (ilim2), nproc2 =", size (ilim2), dc% nproc2
            call finish ('setup_decomp','sizes of ilim1,ilim2')
         end if
         if (any (ilim1(0:dc% nproc1-1) > ilim1(1:dc% nproc1))) then
            write(*,*) 'ilim1=', ilim1
            call finish ('setup_decomp','ilim1')
         end if
         if (any (ilim2(0:dc% nproc2-1) > ilim2(1:dc% nproc2))) then
            write(*,*) 'ilim2=', ilim2
            call finish ('setup_decomp','ilim2')
         end if
         dc% ilim1 = ilim1
         dc% ilim2 = ilim2
      end if

!      IF (dc% myproc == 0) THEN
!        WRITE (0,*) 'ilim1: ',dc% ilim1
!        WRITE (0,*) 'ilim2: ',dc% ilim2
!      ENDIF
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! The following relationships only hold for dc% comm == MPI_COMM_WORLD
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !------------------
      ! Set ng1s ... ng2e
      !------------------
      dc% my_num1 = mod(dc% myproc,dc% nproc1)
      dc% my_num2 = dc% myproc/dc% nproc1

      lb(1) = dc% ilim1(dc% my_num1)
      ub(1) = dc% ilim1(dc% my_num1+1) - 1

      lb(2) = dc% ilim2(dc% my_num2)
      ub(2) = dc% ilim2(dc% my_num2+1) - 1

      dc% ilimp(0) = lb(2)
      dc% ilimp(1) = ub(2)+1
      !--------------------------------
      ! set pe_12 (rank of PEs on grid)
      !--------------------------------
      do i=1, dc% nproc1
        do j=1, dc% nproc2
          dc% pe_12(i,j) = (j-1) * dc% nproc1 + (i-1)
        end do
      end do
      !--------------
      ! Set neighbors
      !--------------
      dc% nbpe                                  = dc% myproc-1;
      IF (dc% my_num1 == 0) dc% nbpe            = -1
      dc% nbaw                                  = dc% myproc+1
      IF (dc% my_num1 == dc% nproc1-1) dc% nbaw = -1
      dc% nbpw                                  = dc% myproc-dc% nproc1
      IF (dc% my_num2 == 0) dc% nbpw            = -1
      dc% nbae                                  = dc% myproc+dc% nproc1
      IF (dc% my_num2 == dc% nproc2-1) dc% nbae = -1
      !----------
      ! my limits
      !----------
      ilimits(1) = lb(1)
      ilimits(2) = ub(1)
      ilimits(3) = lb(2)
      ilimits(4) = ub(2)
      ilimits(5) = (ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)*10

      if (dc% npe > 1) then
#ifndef NOMPI
        CALL MPI_ALLGATHER (ilimits,       5, MPI_INTEGER,  &
           dc% isubpositions, 5, MPI_INTEGER,   &
           dc% comm, mpi_err)
#endif
      else
        dc% isubpositions(:,0) = ilimits
      endif

      ! Determine maximal size of a subdomain

!      IF (dc% myproc == 0) THEN
        dc% imaxnumgp = MAXVAL (dc% isubpositions(5,:))
!      ENDIF
!
!      ! broadcast to all other PEs
!
!      CALL p_bcast (dc% imaxnumgp, 0)
!
!      CALL p_barrier (dc% comm)

    ELSE

      ! Decomposition for Parallel Vector

      ! Check nproc1, nproc2

      IF (dc% nproc1/=1 .OR. dc% nproc2/=1) THEN
        WRITE (0,*) 'nproc1 = ',dc% nproc1,' dc% nproc2 = ',dc% nproc2,  &
             ' both must be 1 for PVP version!'
        call finish('setup_decomp','nproc1,nproc2 must be 1 for PVP version')
      ENDIF

      if (present (ilim1) .or. present (ilim2)) then
         call finish ('setup_decomp',&
                      'ilim1,ilim2 argument unsupported for Parallel Vector')
      end if

      ! Allocate and set decomposition limits

      ALLOCATE (dc% ilim1(0:1))
      ALLOCATE (dc% ilim2(0:1))
      ALLOCATE (dc% ilimp(0:dc% nproc))

      dc% ilim1(0) = 0
      dc% ilim1(1) = ni+1

      dc% ilim2(0) = 1
      dc% ilim2(1) = ni+2

      irows = (ni+1)/dc% nproc
      irem  = MOD(ni+1,dc% nproc)

      dc% ilimp(0) = 1
      DO i=1,irem
        dc% ilimp(i) = dc% ilimp(i-1) + (irows+1)
      ENDDO
      DO i=irem+1,dc% nproc
        dc% ilimp(i) = dc% ilimp(i-1) + irows
      ENDDO

      ! Set ng1s ... ng2e

      dc% my_num1 = 1
      dc% my_num2 = 1

      ! Set isubpositions and dc% imaxnumgp for the sequentiell case

      dc% isubpositions(1,0) = lb(1)
      dc% isubpositions(2,0) = ub(1)
      dc% isubpositions(3,0) = lb(2)
      dc% isubpositions(4,0) = ub(2)
      dc% isubpositions(5,0) = (ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)*10
      dc% imaxnumgp = dc% isubpositions(5,0)

      ! Set neighbors

      dc% nbpe = -1
      dc% nbaw = -1
      dc% nbpw = -1
      dc% nbae = -1

    ENDIF

  END SUBROUTINE setup_decomp
!------------------------------------------------------------------------------
  subroutine print_com (marr, unit)
  integer          ,intent(in)           :: marr (:,:,:,:)
  integer          ,intent(in) ,optional :: unit
    integer :: j
    integer :: iu

    iu = 6; if (present(unit)) iu = unit
    write (iu,'()')
    write (iu,'(a,4i5)') ' size  (marr) =',size(marr)
    write (iu,'(a,4i5)') ' lbound(marr) =',lbound(marr)
    write (iu,'(a,4i5)') ' ubound(marr) =',lbound(marr)
    if (size(marr,2)<=500) then
      write (iu,'(/a/)')' pe:'
      do j = lbound(marr,3), ubound(marr,3)
        write (iu,'(500i1)') marr(1,:,j,1)
      end do
      write (iu,'(/a/)')' i :'
      do j = lbound(marr,3), ubound(marr,3)
        write (iu,'(500i1)') mod(marr(2,:,j,1),10)
      end do
      write (iu,'(/a/)')' j :'
      do j = lbound(marr,3), ubound(marr,3)
        write (iu,'(500i1)') mod(marr(3,:,j,1),10)
      end do
      write (iu,'(/a/)')' d :'
      do j = lbound(marr,3), ubound(marr,3)
        write (iu,'(500i1)') marr(4,:,j,1)
      end do
    endif
  end subroutine print_com
!------------------------------------------------------------------------------
  SUBROUTINE setup_com_reg (dc, marr, lbg, ubg)
  type (t_atm_dec) ,INTENT(inout) :: dc
  INTEGER ,pointer                :: marr(:,:,:,:) ! (4,i,j,d)
  INTEGER          ,INTENT(in)    :: lbg (4)
  INTEGER          ,INTENT(in)    :: ubg (4)
  !---------------------------------------------
  ! set up communication info for a regular grid
  !---------------------------------------------
    INTEGER :: jp1, jp2, j1, j2, jp

    allocate (marr(4,lbg(1):ubg(1),lbg(2):ubg(2),1))
    marr(:,:,:,:) = -1

    DO jp2 = 1, dc% nproc2
      DO j2 = dc% ilim2(jp2-1), dc% ilim2(jp2)-1
        DO jp1 = 1, dc% nproc1
          jp = (jp2-1)*dc% nproc1 + jp1-1
          DO j1 = dc% ilim1(jp1-1), dc% ilim1(jp1)-1
            marr(1,j1,j2,1) = jp
            marr(2,j1,j2,1) = j1
            marr(3,j1,j2,1) = j2
            marr(4,j1,j2,1) = 1
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  end SUBROUTINE setup_com_reg
!------------------------------------------------------------------------------
  SUBROUTINE setup_com_icon (dc, marr, lbg, ubg)
  type (t_atm_dec) ,INTENT(inout) :: dc
  INTEGER ,pointer                :: marr(:,:,:,:) ! (4,i,j,d)
  INTEGER          ,INTENT(in)    :: lbg (4)
  INTEGER          ,INTENT(in)    :: ubg (4)
  !---------------------------------------------------------
  ! set up communication info for the ICON unstructured grid
  !---------------------------------------------------------
    INTEGER :: jp1, jp2, j1, j2, jp

    if (ubg(2) /= lbg(2)) &
         call finish('setup_com_icon','unimplemented/invalid decomposition')

    allocate (marr(4,lbg(1):ubg(1),lbg(2):ubg(2),1))
    marr(:,:,:,:) = -1

    DO jp2 = 1, dc% nproc2
      DO j2 = dc% ilim2(jp2-1), dc% ilim2(jp2)-1
        DO jp1 = 1, dc% nproc1
          jp = (jp2-1)*dc% nproc1 + jp1-1
          DO j1 = dc% ilim1(jp1-1), dc% ilim1(jp1)-1
            marr(1,j1,j2,1) = jp
            marr(2,j1,j2,1) = j1
            marr(3,j1,j2,1) = j2
            marr(4,j1,j2,1) = 1
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  end SUBROUTINE setup_com_icon
!------------------------------------------------------------------------------
  SUBROUTINE setup_com_tri (dc, marr,                        &
                            kgg1s, kgg1e, kgg2s, kgg2e, knd, &
                            kg1s,  kg1e,  kg2s,  kg2e,  kni  )

    type (t_atm_dec) ,INTENT(inout) :: dc

    INTEGER ,pointer :: marr (:,:,:,:) ! (4,i,j,d)

    INTEGER, INTENT(in) :: kgg1s, kgg1e, kgg2s, kgg2e, knd, &
                           kg1s,  kg1e,  kg2s,  kg2e,  kni

    INTEGER :: jd, jp1, jp2, j1, j2, jp, jb, jr, j
    INTEGER :: mbsize, mb
    INTEGER :: mi1sc, mi1ec, mi2sc, mi2ec

#if 0   /* needed for original code */
#ifndef NOMPI
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
#endif
#endif

    INTEGER :: idx_bound(3,(kg1e-kg1s+5)*(kg2e-kg2s+5) )

    INTEGER, ALLOCATABLE :: idx_help(:,:)

#ifndef NOMPI
    INTEGER              :: np                          ! Total PEs this comm.
    INTEGER, ALLOCATABLE :: rcount(:), scount(:)        ! Counts
    INTEGER, ALLOCATABLE :: rdispl(:), sdispl(:)        ! Displacements
#endif

    ! Section 0:
    !
    ! Allocate arrays

    allocate(marr(4,kgg1s-2:kgg1e+2,kgg2s-2:kgg2e+2,1:knd))

    ALLOCATE(dc% np_recv_s(0:dc% nproc1*dc% nproc2-1), &
             dc% np_recv_1(0:dc% nproc1*dc% nproc2-1), &
             dc% np_recv_2(0:dc% nproc1*dc% nproc2-1), &
             dc% np_send_s(0:dc% nproc1*dc% nproc2-1), &
             dc% np_send_1(0:dc% nproc1*dc% nproc2-1), &
             dc% np_send_2(0:dc% nproc1*dc% nproc2-1), &
             dc% np_recv_t(0:dc% nproc1*dc% nproc2-1))

    ! Section 1:
    !
    ! Set up the array marr for all complete ( = core + extended)
    ! diamonds so that it can tell us later exactly for every point,
    ! from which processor and which diamond at which location
    ! a particular point in the extension area has to be picked.
    !
    ! marr(1,.,.,.) processor number
    ! marr(2,.,.,.) j1-index
    ! marr(3,.,.,.) j2-index
    ! marr(4,.,.,.) diamond number

    ! Set array marr to -1 (just for debugging)

    marr(:,:,:,:) = -1

    ! Set up marr at the core of the diamonds

#ifdef DEBUG
    WRITE (0,'(a,i0)') 'setup_com: nproc1 = ', dc% nproc1
    WRITE (0,'(a,i0)') 'setup_com: nproc2 = ', dc% nproc2

    WRITE (0,'(a,64i4)') 'setup_com: Index limits in direction 1 = ', &
         dc% ilim1(:)
    WRITE (0,'(a,64i4)') 'setup_com: Index limits in direction 2 = ', &
         dc% ilim2(:)
#endif

    DO jd = 1, knd
      DO jp2 = 1, dc% nproc2
        DO j2 = dc% ilim2(jp2-1), MIN(dc% ilim2(jp2)-1,kni)
          DO jp1 = 1, dc% nproc1
            DO j1 = MAX(dc% ilim1(jp1-1),1), dc% ilim1(jp1)-1
              jp = (jp2-1)*dc% nproc1 + jp1-1
              marr(1,j1,j2,jd) = jp
              marr(2,j1,j2,jd) = j1
              marr(3,j1,j2,jd) = j2
              marr(4,j1,j2,jd) = jd
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Set the pole

    DO jd =  1, knd
      marr(1,0,1,jd) = 0 ! Pole is always owned by proc 0
      marr(2,0,1,jd) = 0
      marr(3,0,1,jd) = 1
      IF(jd <= 5) THEN
        marr(4,0,1,jd) = 1 ! Use always diamond 1 for north pole
      ELSE
        marr(4,0,1,jd) = 6 ! Use always diamond 6 for south pole
      ENDIF
    ENDDO

    ! Set boundaries

    DO jd = 1, knd

      ! poleward west

      DO j= 1, kni
        marr(:,j  , 0,jd) = marr(:,1,j,mpw(jd))
        marr(:,j+1,-1,jd) = marr(:,2,j,mpw(jd))
      ENDDO

      ! poleward east

      DO j= 1, kni
        marr(:, 0,j+1,jd) = marr(:,j,1,mpe(jd))
        marr(:,-1,j+2,jd) = marr(:,j,2,mpe(jd))
        marr(:,-2,j+3,jd) = marr(:,j,3,mpe(jd))
      ENDDO

      ! anti-poleward west

      DO j= 1, kni
        marr(:,kni+1,kni-j+1,jd) = marr(:,j,kni  ,maw(jd))
        marr(:,kni+2,kni-j+1,jd) = marr(:,j,kni-1,maw(jd))
      ENDDO

      ! anti-poleward east

      DO j= 1, kni
        marr(:,kni-j+1,kni+1,jd) = marr(:,kni  ,j,mae(jd))
        marr(:,kni-j+1,kni+2,jd) = marr(:,kni-1,j,mae(jd))
        marr(:,kni-j+1,kni+3,jd) = marr(:,kni-2,j,mae(jd))
      ENDDO

    ENDDO

    ! Set special points

    DO jd = 1, knd

      ! Pole

      marr(:,-1, 2,jd) = marr(:,1,1,mpe(mpe(jd)))
      marr(:,-2, 3,jd) = marr(:,2,1,mpe(mpe(jd)))
      marr(:,-2, 2,jd) = marr(:,1,2,mpe(mpe(jd)))

      marr(:, 0, 0,jd) = marr(:,1,1,mpw(mpw(jd)))
      marr(:, 0,-1,jd) = marr(:,2,1,mpw(mpw(jd)))
      marr(:, 1,-1,jd) = marr(:,1,2,mpw(mpw(jd)))

      marr(:,-1, 0,jd) = marr(:, 0, 0,jd) ! Undefined
      marr(:,-1, 1,jd) = marr(:, 0, 0,jd) ! Mirror

      marr(:,-1,-1,jd) = marr(:, 0,-1,jd) ! Undefined
      marr(:,-2,-1,jd) = marr(:, 0,-1,jd) ! Undefined
      marr(:,-2, 0,jd) = marr(:, 0,-1,jd) ! Undefined
      marr(:,-2, 1,jd) = marr(:, 0,-1,jd) ! Mirror

      ! West Corner

      marr(:,kni+1, 0,jd) = marr(:,kni,0,jd) ! Mirror

      marr(:,kni+2,-1,jd) = marr(:,kni+1,-1,jd) ! Undefined
      marr(:,kni+2, 0,jd) = marr(:,kni+1,-1,jd) ! Mirror

      ! Antipole

      marr(:,kni+1,kni+1,jd) = marr(:,kni  ,kni+2,jd) ! Mirror
      marr(:,kni+1,kni+2,jd) = marr(:,kni  ,kni+2,jd) ! Undefined

      marr(:,kni+2,kni+1,jd) = marr(:,kni  ,kni+3,jd) ! Mirror
      marr(:,kni+2,kni+2,jd) = marr(:,kni  ,kni+3,jd) ! Undefined
      marr(:,kni+2,kni+3,jd) = marr(:,kni  ,kni+3,jd) ! Undefined
      marr(:,kni+1,kni+3,jd) = marr(:,kni  ,kni+3,jd) ! Undefined

      ! East Corner

      marr(:, 0,kni+2,jd) = marr(:,-1,kni+2,jd) ! Copy
      marr(:,-1,kni+2,jd) = marr(:,-1,kni+1,jd) ! Mirror

      marr(:, 0,kni+3,jd) = marr(:,-2,kni+3,jd) ! Copy
      marr(:,-2,kni+3,jd) = marr(:,-2,kni+2,jd) ! Undefined
      marr(:,-1,kni+3,jd) = marr(:,-2,kni+2,jd) ! Mirror

    ENDDO

#ifdef DEBUG
    if (dace% peio) then
      DO jp = 0, 0 ! dc% nproc1*dc% nproc2-1
        WRITE (0,'(a,i0)') 'setup_com: status of PE: ', jp
        DO jd = 1, 1 ! knd
          WRITE (0,'(a,i0)') 'setup_com: - for diamond ', jd
          DO j1 = kgg1s-2, kgg1e+2
            WRITE (0,'(200i1)') marr(1,j1,:,jd)
          ENDDO
        ENDDO
      ENDDO
    endif
#endif

    ! Section 2:
    !
    ! Gather all our boundary points we need during extension
    ! in array idx_bound.
    !
    ! This is done only for 1 diamond since all other diamonds
    ! are identical
    !
    ! idx_bound(1,.) = j1-index of boundary point
    ! idx_bound(2,.) = j2-index of boundary point
    ! idx_bound(3,.) == 1 if this point is needed in a 1 line exchange
    ! idx_bound(3,.) == 2 if this point is needed in a 2 line exchange
    ! idx_bound(3,.) == 3 these points are not really needed
    !                     they are exchanged in a 2 line exchange
    !                     for compatibility purposes only

    ! Set the computational boundaries

    mi1sc = MAX(kg1s,1)
    mi1ec = kg1e
    mi2sc = kg2s
    mi2ec = MIN(kg2e,kni)

    mb = 0

    DO j1 = kg1s-2, kg1e+2
      DO j2 = kg2s-2, kg2e+2

        IF (j1 < mi1sc .OR. j1 > mi1ec .OR. j2 < mi2sc .OR. j2 > mi2ec) THEN

          ! This is a boundary point

          mb = mb+1
          idx_bound(1,mb) = j1
          idx_bound(2,mb) = j2

          IF (j1 < mi1sc-2 .OR. j1 > mi1ec+2 .OR. &
              j2 < mi2sc-2 .OR. j2 > mi2ec+2) THEN
            idx_bound(3,mb) = 3
          ELSE IF (j1 < mi1sc-1 .OR. j1 > mi1ec+1 .OR.   &
                   j2 < mi2sc-1 .OR. j2 > mi2ec+1) THEN
            idx_bound(3,mb) = 2
          ELSE
            idx_bound(3,mb) = 1
          ENDIF

        ENDIF

      ENDDO
    ENDDO

    mbsize = mb

    ! The processor which owns the pole (always 0) needs
    ! some special treatment since it needs some more outer
    ! points in order to calculate the stencils at the pole
    ! correctly

    IF (dc% myproc == 0) THEN
      DO j = 1, mbsize
        IF ( idx_bound(1,j) == -1 .AND. idx_bound(2,j) <= 2 )   &
             idx_bound(3,j) = 1 ! was previously set to 2
        IF ( idx_bound(1,j) == -2 .AND. idx_bound(2,j) <= 3 )   &
             idx_bound(3,j) = 2 ! was previously set to 3
      ENDDO
    ENDIF

#ifdef DEBUG
    DO jp = 1,mbsize
      WRITE (0,'(a,3i5)') 'setup_com: index bounds = ', &
           idx_bound(1:3,jp)
    ENDDO
#endif

    ! Section 3:
    !
    ! Figure out which of our boundary points we have to receive from
    ! which processor, sort them accordingly by processor and
    ! set the receive arrays.

    ! Allocate idx_recv and idx_help

    ALLOCATE ( dc% idx_recv(3,10*mbsize) )
    ALLOCATE ( idx_help(3,10*mbsize) )

    ! Sort boundary points by processor to receive them from

    mb = 0

    DO jp = 0, dc% nproc1*dc% nproc2-1   ! Loop over all processors

      dc% np_recv_s(jp) = mb+1

      DO jr = 1, 2 ! Loop over the 2 boundary line cases
        DO jd = 1, knd ! Loop over the diamonds
          DO jb = 1, mbsize
            j1 = idx_bound(1,jb)
            j2 = idx_bound(2,jb)

            IF (marr(1,j1,j2,jd) /= jp) CYCLE

            IF ( (jr==1 .AND. idx_bound(3,jb)==1) .OR.   &
                 (jr==2 .AND. idx_bound(3,jb)/=1) ) THEN
              mb = mb+1
              dc% idx_recv(1,mb) = j1
              dc% idx_recv(2,mb) = j2
              dc% idx_recv(3,mb) = jd

              idx_help(1,mb) = marr(2,j1,j2,jd)
              idx_help(2,mb) = marr(3,j1,j2,jd)
              idx_help(3,mb) = marr(4,j1,j2,jd)
            ENDIF
          ENDDO
        ENDDO

        IF (jr==1) THEN
          dc% np_recv_1(jp) = mb - dc% np_recv_s(jp) + 1
        ELSE
          dc% np_recv_2(jp) = mb - dc% np_recv_s(jp) + 1
        ENDIF

      ENDDO
    ENDDO

    dc% np_recv_tot = mb

    ! Internal check:

    IF (dc% np_recv_tot /= 10*mbsize) THEN
      WRITE (0,*) 'set_boundaries: Internal error!'
      call finish ('setup_com','set_boundaries: Internal error')
    ENDIF

    ! Section 4:
    !
    ! Now every processor knows, which and how many boundary points
    ! it has to receive from which other processor.
    ! It does not yet know, however, how many and which interior points
    ! it has to send to others.
    !
    ! To figure that out, every processor first sends to all others
    ! how many points it needs and then sends the indices of the points

    ! Send how many points we need

    IF ( dc% nproc1*dc% nproc2 > 1) THEN ! no communication if only 1 processor

#if 0   /* original code */

      DO jp = 0, dc% nproc1*dc% nproc2-1   ! Loop over all processors

        ! dc% np_recv_x array -> dc% np_send_x array-element of
        !                    the corresponding processor

#ifndef NOMPI
        CALL MPI_Scatter(dc% np_recv_1, 1, MPI_INTEGER, dc% np_send_1(jp), 1,   &
             MPI_INTEGER, jp, dc% comm, mpi_err)
        CALL MPI_Scatter(dc% np_recv_2, 1, MPI_INTEGER, dc% np_send_2(jp), 1,   &
             MPI_INTEGER, jp, dc% comm, mpi_err)
        CALL MPI_Scatter(dc% np_recv_s, 1, MPI_INTEGER, dc% np_recv_t(jp), 1,   &
             MPI_INTEGER, jp, dc% comm, mpi_err)
#endif

      ENDDO

#else   /* faster code using MPI collectives */

#ifndef NOMPI
      CALL MPI_Alltoall (dc% np_recv_1, 1, MPI_INTEGER,                  &
                         dc% np_send_1, 1, MPI_INTEGER, dc% comm, mpi_err)
      CALL MPI_Alltoall (dc% np_recv_2, 1, MPI_INTEGER,                  &
                         dc% np_send_2, 1, MPI_INTEGER, dc% comm, mpi_err)
      CALL MPI_Alltoall (dc% np_recv_s, 1, MPI_INTEGER,                  &
                         dc% np_recv_t, 1, MPI_INTEGER, dc% comm, mpi_err)
#endif

#endif  /* */

    ELSE

      dc% np_send_1(0) = dc% np_recv_1(0)
      dc% np_send_2(0) = dc% np_recv_2(0)
      dc% np_recv_t(0) = dc% np_recv_s(0)

    ENDIF

    ! Setup dc% np_send_s, count total number of points to send

    mb = 0

    DO jp = 0, dc% nproc1*dc% nproc2-1
      dc% np_send_s(jp) = mb + 1
      mb = mb + dc% np_send_2(jp)
    ENDDO

    dc% np_send_tot = mb

    ! Allocate idx_send

    ALLOCATE ( dc% idx_send(3,dc% np_send_tot) )

    ! Get the indices of the points we have to send

#if 0   /* original code */

    DO jp = 0, dc% nproc1*dc% nproc2-1   ! Loop over all processors

      IF (jp == dc% myproc) THEN

        ! It is our turn to send which points we need

        DO jp1 = 0, dc% nproc1*dc% nproc2-1

          IF ( jp1 == dc% myproc ) THEN

            ! Don't send to ourself, just copy

            DO jb = 1, dc% np_recv_2(dc% myproc)
              j1 = dc% np_send_s(dc% myproc) + jb - 1
              j2 = dc% np_recv_s(dc% myproc) + jb - 1
              dc% idx_send(1,j1) = idx_help(1,j2)
              dc% idx_send(2,j1) = idx_help(2,j2)
              dc% idx_send(3,j1) = idx_help(3,j2)
            ENDDO

          ELSE

            ! Send what we want to receive

#ifndef NOMPI
            IF (dc% np_recv_2(jp1) > 0)                           &
                 CALL MPI_Send(idx_help(1,dc% np_recv_s(jp1)),    &
                      3*dc% np_recv_2(jp1), MPI_INTEGER, jp1, 3,  &
                      dc% comm, mpi_err)
#endif
          ENDIF

        ENDDO

      ELSE

        ! It is our turn to receive which points are needed from
        ! processor jp

#ifndef NOMPI
        IF (dc% np_send_2(jp) > 0)                          &
            CALL MPI_Recv(dc% idx_send(1,dc% np_send_s(jp)),    &
                 3*dc% np_send_2(jp), MPI_INTEGER, jp, 3,   &
                 dc% comm, mpi_status, mpi_err)
#endif

      ENDIF

#ifndef NOMPI
      IF (lmpi) CALL p_barrier(dc% comm)
#endif

    ENDDO

#else   /* faster code using MPI collectives */

#ifndef NOMPI
    np = dc% nproc1 * dc% nproc2
    IF (np > 1) THEN                           ! skip if only 1 processor

       ALLOCATE (rcount(0:np-1), rdispl(0:np-1), &
                 scount(0:np-1), sdispl(0:np-1)  )

       scount(:) = 3 *  dc% np_recv_2(:)        ! Send counts
       rcount(:) = 3 *  dc% np_send_2(:)        ! Recv counts
       sdispl(:) = 3 * (dc% np_recv_s(:) - 1)   ! Sendbuffer displacements
       rdispl(:) = 3 * (dc% np_send_s(:) - 1)   ! Recvbuffer displacements

       if (sdispl(0) /= 0 .or. rdispl(0) /= 0) &
            call finish ('setup_com_tri','Internal error deriving displacements')

       CALL MPI_Alltoallv (    idx_help, scount, sdispl, MPI_INTEGER, &
                           dc% idx_send, rcount, rdispl, MPI_INTEGER, &
                           dc% comm, mpi_err                          )

    END IF
#endif

#endif

    ! Section 5:
    !
    ! Final work

    ! We are all done, deallocate idx_help

    DEALLOCATE ( idx_help )

    ! Calculate maximum of points to send/receive over all processors

    IF (lmpi) THEN
      dc% np_recv_max = p_max (dc% np_recv_tot, dc% comm)
      dc% np_send_max = p_max (dc% np_send_tot, dc% comm)
    ELSE
      dc% np_recv_max = dc% np_recv_tot
      dc% np_send_max = dc% np_send_tot
    ENDIF

!if (dace% peio) call print_decomp (dc)
!if (dace% peio) call print_com    (marr)

  END SUBROUTINE setup_com_tri
!------------------------------------------------------------------------------
END MODULE mo_atm_decomp
