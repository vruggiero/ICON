!
!+ Sparse matrix conversion and matrix-vector products from SPARSKIT
!
MODULE mo_sparskit
!
! Description:
!   This module contains selected sparse matrix conversion routines
!   and reference implementations of sparse matrix-vector products
!   from SPARSKIT:
!   http://www-users.cs.umn.edu/~saad/software/SPARSKIT/sparskit.html
!
!   Copyrights of the original code:
!
!   Copyright (C) 2005, the Regents of the University of Minnesota
!
!   SPARSKIT is  free software; you  can redistribute it and/or  modify it
!   under the terms of the  GNU Lesser General Public License as published
!   by the  Free Software Foundation [version  2.1 of the  License, or any
!   later version.]
!  ----------------------------------------------------------------------
!   The original Fortran77 code was converted using TO_F90 by Alan Miller
!   and manually fixed.
!
!   The subset of routines retained is described below.
!
! Current Code Owner: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_5         2009/05/25 Harald Anlauf
!  amuxj_vec_*: optimized for SX-9 by Rudolph Fischer and Jens-Olaf Beismann (NEC)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_19        2012-04-16 Harald Anlauf
!  optimize for sxf90 rev.430+
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Yousef Saad (saad@cs.umn.edu)  Original code, 1986(?)-2005.
! Harald Anlauf, DWD 2008/02     Adaptions to Fortran 95
!============================================================================

  !----------------------------------------------------------------------c
  !                          S P A R S K I T                             c
  !----------------------------------------------------------------------c
  !          BASIC MATRIX-VECTOR OPERATIONS - MATVEC MODULE              c
  !         Matrix-vector Multiplications and Triang. Solves             c
  !----------------------------------------------------------------------c
  ! 1) Matrix-vector products:                                           c
  !---------------------------                                           c
  ! amux  : A times a vector. Compressed Sparse Row (CSR) format.        c
  ! amuxms: A times a vector. Modified Compressed Sparse Row format.     c
  ! atmux : Transp(A) times a vector. CSR format.                        c
  ! atmuxr: Transp(A) times a vector. CSR format. A rectangular.         c
  ! amuxj : A times a vector. Jagged Diagonal (JAD) format.              c
  !----------------------------------------------------------------------c

  !----------------------------------------------------------------------c
  !                    FORMAT CONVERSION MODULE                          c
  !----------------------------------------------------------------------c
  ! csrdns  : converts a row-stored sparse matrix into the dense format. c
  ! dnscsr  : converts a dense matrix to a sparse storage format.        c
  ! coocsr  : converts coordinate to  to csr format                      c
  ! coicsr  : in-place conversion of coordinate to csr format            c
  ! csrcoo  : converts compressed sparse row to coordinate.              c
  ! csrssr  : converts compressed sparse row to symmetric sparse row     c
  ! ssrcsr  : converts symmetric sparse row to compressed sparse row     c
  ! csrmsr  : converts compressed sparse row format to modified sparse   c
  !           row format                                                 c
  ! msrcsr  : converts modified sparse row format to compressed sparse   c
  !           row format.                                                c
  ! csrcsc  : converts compressed sparse row format to compressed sparse c
  !           column format (transposition)                              c
  ! csrcsc2 : rectangular version of csrcsc                              c
  ! csrbnd  : converts a compressed sparse row format into a banded      c
  !           format (linpack style).                                    c
  ! bndcsr  : converts a banded format (linpack style) into a compressed c
  !           sparse row storage.                                        c
  ! csrjad  : converts the csr format into the jagged diagonal format    c
  ! jadcsr  : converts the jagged-diagonal format into the csr format    c
  ! csorted : Checks if matrix in CSR format is sorted by columns        c
  !--------- miscalleneous additions not involving the csr format--------c
  ! dcsort  : sorting routine used by crsjad                             c
  !----------------------------------------------------------------------c

  !----------------------------------------------------------------------c
  !                     UNARY SUBROUTINES MODULE                         c
  !----------------------------------------------------------------------c
  ! dvperm : permutes a real vector (in-place)                           c
  ! getbwd : returns the bandwidth information on a matrix.              c
  !----------------------------------------------------------------------c

  !----------------------------------------------------------------------c
  !                        INPUT-OUTPUT MODULE                           c
  !----------------------------------------------------------------------c
  !  dump   : outputs matrix rows in a simple format (debugging purposes)c
  !----------------------------------------------------------------------c

  !============================================================================
  !-------------
  ! Modules used
  !-------------
  use mo_kind, only: wp, sp
  !============================================================================
  implicit none
  !----------------
  ! Public entities
  !----------------
  public                                ! Everything in this module is public

  !===========
  ! Interfaces
  !===========
  interface csrjad
     module procedure csrjad            ! Convert CSR -> JAD
     module procedure csrjad_r4
  end interface

  interface jadcsr
     module procedure jadcsr            ! Convert JAD -> CSR
     module procedure jadcsr_r4
  end interface

  interface csrcsc
     module procedure csrcsc            ! Convert CSR <-> CSC
     module procedure csrcsc_r4
  end interface

  !--------------------------
  ! Matvec JAD * x = Perm * y
  !--------------------------
  interface amuxj
!    module procedure amuxj_r8          ! Original version (single diagonal)
!    module procedure amuxj_r4_r8       ! ..real(4) matrix, real(8) vectors
!    module procedure amuxj_grp_r8      ! Grouped version (up to 4 diagonals)
!    module procedure amuxj_grp_r4_r8   ! ..real(4) matrix, real(8) vectors
     module procedure amuxj_vec_r8      ! Optimized for NEC SX-9
     module procedure amuxj_vec_r4_r8   ! ..real(4) matrix, real(8) vectors
  end interface

  !============================================================================
contains
  !============================================================================
  subroutine amux (n, x, y, a, ja, ia)
    integer,  intent(in)    :: n        ! row dimension of A
    real(wp), intent(in)    :: x(:)
    real(wp), intent(out)   :: y(:)
    real(wp), intent(in)    :: a(:)
    integer,  intent(in)    :: ja(:)
    integer,  intent(in)    :: ia(:)
    !-----------------------------------------------------------------------
    !         A times a vector
    !-----------------------------------------------------------------------
    ! multiplies a matrix by a vector using the dot product form
    ! Matrix A is stored in compressed sparse row storage.
    !
    ! on entry:
    !----------
    ! n     = row dimension of A
    ! x     = real array of length equal to the column dimension of
    !         the A matrix.
    ! a, ja,
    !    ia = input matrix in compressed sparse row format.
    !
    ! on return:
    !-----------
    ! y     = real array of length n, containing the product y=Ax
    !-----------------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------------
    integer  :: i, k
    real(wp) :: t
    !-----------------------------------------------------------------------
    do  i = 1,n

       !     compute the inner product of row i with vector x

       t = 0.0d0
       do  k=ia(i), ia(i+1)-1
          t = t + a(k)*x(ja(k))
       end do

       !     store result in y(i)

       y(i) = t
    end do
    !---------end-of-amux---------------------------------------------------
    !-----------------------------------------------------------------------
  end subroutine amux
  !-----------------------------------------------------------------------
  subroutine amuxms (n, x, y, a, ja)
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: x(:)
    real(wp), intent(out) :: y(:)
    real(wp), intent(in)  :: a(:)
    integer,  intent(in)  :: ja(:)
    !-----------------------------------------------------------------------
    !         A times a vector in MSR format
    !-----------------------------------------------------------------------
    ! multiplies a matrix by a vector using the dot product form
    ! Matrix A is stored in Modified Sparse Row storage.
    !
    ! on entry:
    !----------
    ! n     = row dimension of A
    ! x     = real array of length equal to the column dimension of
    !         the A matrix.
    ! a, ja,= input matrix in modified compressed sparse row format.
    !
    ! on return:
    !-----------
    ! y     = real array of length n, containing the product y=Ax
    !-----------------------------------------------------------------------
    ! local variables

    integer :: i, k
    !-----------------------------------------------------------------------
    do  i=1, n
       y(i) = a(i)*x(i)
    end do
    do  i = 1,n

       !     compute the inner product of row i with vector x

       do  k=ja(i), ja(i+1)-1
          y(i) = y(i) + a(k) *x(ja(k))
       end do
    end do
    !---------end-of-amuxm--------------------------------------------------
    !-----------------------------------------------------------------------
  end subroutine amuxms
  !-----------------------------------------------------------------------
  subroutine atmux (n, x, y, a, ja, ia)
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: x(:)
    real(wp), intent(out) :: y(:)
    real(wp), intent(in)  :: a(:)
    integer,  intent(in)  :: ja(:)
    integer,  intent(in)  :: ia(:)
    !-----------------------------------------------------------------------
    !         transp( A ) times a vector
    !-----------------------------------------------------------------------
    ! multiplies the transpose of a matrix by a vector when the original
    ! matrix is stored in compressed sparse row storage. Can also be
    ! viewed as the product of a matrix by a vector when the original
    ! matrix is stored in the compressed sparse column format.
    !-----------------------------------------------------------------------
    !
    ! on entry:
    !----------
    ! n     = row dimension of A
    ! x     = real array of length equal to the column dimension of
    !         the A matrix.
    ! a, ja,
    !    ia = input matrix in compressed sparse row format.
    !
    ! on return:
    !-----------
    ! y     = real array of length n, containing the product y=transp(A)*x
    !-----------------------------------------------------------------------
    !     local variables

    integer :: i, k
    !-----------------------------------------------------------------------
    y(1:n) = 0          ! zero out output vector

    ! loop over the rows

    do  i = 1,n
!NEC$ ivdep
       do  k=ia(i), ia(i+1)-1
          y(ja(k)) = y(ja(k)) + x(i)*a(k)
       end do
    end do
    !-------------end-of-atmux----------------------------------------------
    !-----------------------------------------------------------------------
  end subroutine atmux
  !-----------------------------------------------------------------------
  subroutine atmuxr (m, n, x, y, a, ja, ia)
    integer,  intent(in)  :: m
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: x(:)
    real(wp), intent(out) :: y(:)
    real(wp), intent(in)  :: a(:)
    integer,  intent(in)  :: ja(:)
    integer,  intent(in)  :: ia(:)
    !-----------------------------------------------------------------------
    !         transp( A ) times a vector, A can be rectangular
    !-----------------------------------------------------------------------
    ! See also atmux.  The essential difference is how the solution vector
    ! is initially zeroed.  If using this to multiply rectangular CSC
    ! matrices by a vector, m number of rows, n is number of columns.
    !-----------------------------------------------------------------------
    !
    ! on entry:
    !----------
    ! m     = column dimension of A
    ! n     = row dimension of A
    ! x     = real array of length equal to the column dimension of
    !         the A matrix.
    ! a, ja,
    !    ia = input matrix in compressed sparse row format.
    !
    ! on return:
    !-----------
    ! y     = real array of length n, containing the product y=transp(A)*x
    !-----------------------------------------------------------------------
    !     local variables

    integer :: i, k
    !-----------------------------------------------------------------------
    y(1:m) = 0          ! zero out output vector

    ! loop over the rows

    do  i = 1,n
!NEC$ ivdep
       do  k=ia(i), ia(i+1)-1
          y(ja(k)) = y(ja(k)) + x(i)*a(k)
       end do
    end do
    !-------------end-of-atmuxr---------------------------------------------
    !-----------------------------------------------------------------------
  end subroutine atmuxr
  !-----------------------------------------------------------------------
  subroutine amuxj_r8 (n, x, y, jdiag, a, ja, ia)
    integer,  intent(in)      :: n        ! Row dimension of A
    real(wp), intent(in)      :: x(n)
    real(wp), intent(out)     :: y(n)
    integer,  intent(in)      :: jdiag    ! Number of jagged-diagonals
    real(wp), intent(in)      :: a(:)
    integer,  intent(in)      :: ja(:)
    integer,  intent(in)      :: ia(:)
    !-----------------------------------------------------------------------
    !        A times a vector in Jagged-Diagonal storage format (JAD)
    !-----------------------------------------------------------------------
    ! multiplies a matrix by a vector when the original matrix is stored
    ! in the jagged diagonal storage format.
    !-----------------------------------------------------------------------
    ! on entry:
    !----------
    ! n      = row dimension of A
    ! x      = real array of length equal to the column dimension of
    !          the A matrix.
    ! jdiag  = integer. The number of jagged-diagonals in the data-structure.
    ! a      = real array containing the jagged diagonals of A stored
    !          in succession (in decreasing lengths)
    ! j      = integer array containing the column indices of the
    !          corresponding elements in a.
    ! ia     = integer array containing the lengths of the jagged diagonals
    !
    ! on return:
    !-----------
    ! y      = real array of length n, containing the product y=A*x
    !
    ! Note:
    !-------
    ! Permutation related to the JAD format is not performed.
    ! this can be done by:
    !     call dvperm (n, y, iperm)
    ! after the call to amuxj, where iperm is the permutation produced
    ! by csrjad.
    !-----------------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------------
    integer :: ii, k1, diag_length, j

    y(1:n) = 0
    outer: do  ii=1, jdiag
       k1 = ia(ii) - 1
       diag_length = ia(ii+1) - k1 - 1
!      if (diag_length == 0) exit outer         ! Prevents vectorization!?
       do  j=1, diag_length
          y(j) = y(j) + a(k1+j) * x(ja(k1+j))
       end do
    end do outer
    !----------end-of-amuxj-------------------------------------------------
    !-----------------------------------------------------------------------
  end subroutine amuxj_r8
  !-----------------------------------------------------------------------
  subroutine amuxj_r4_r8 (n, x, y, jdiag, a, ja, ia)
    integer,  intent(in)      :: n        ! Row dimension of A
    real(wp), intent(in)      :: x(n)
    real(wp), intent(out)     :: y(n)
    integer,  intent(in)      :: jdiag    ! Number of jagged-diagonals
    real(sp), intent(in)      :: a(:)
    integer,  intent(in)      :: ja(:)
    integer,  intent(in)      :: ia(:)

    integer :: ii, k1, diag_length, j

    y(1:n) = 0
    outer: do  ii=1, jdiag
       k1 = ia(ii) - 1
       diag_length = ia(ii+1) - k1 - 1
       do  j=1, diag_length
          y(j) = y(j) + a(k1+j) * x(ja(k1+j))
       end do
    end do outer
  end subroutine amuxj_r4_r8
  !-----------------------------------------------------------------------
  subroutine amuxj_grp_r8 (n, x, y, jdiag, a, ja, ia)
    integer,  intent(in)      :: n        ! Row dimension of A
    real(wp), intent(in)      :: x(n)
    real(wp), intent(out)     :: y(n)
    integer,  intent(in)      :: jdiag    ! Number of jagged-diagonals
    real(wp), intent(in)      :: a(:)
    integer,  intent(in)      :: ja(:)
    integer,  intent(in)      :: ia(:)
    !-----------------------------------------------------------------------
    !        A times a vector in Jagged-Diagonal storage format (JAD)
    !-----------------------------------------------------------------------
    ! "Grouped" implementation to simultaneously handle up to 4 diagonals
    ! of same length.
    ! This is supposed to give better performance on vector processors [1].
    !
    ! [1] Sunil R. Tiyyagura, Uwe Küster, Stefan Borowski,
    !     "Performance Improvement of Sparse Matrix Vector Product
    !      on Vector Machines", ICCS '06.
    !-----------------------------------------------------------------------
    ! Machine-dependent tuning parameter:
    ! maximum no. of simultaneous diagonals (range: 1-4)
    !
    integer, parameter :: MAXDIAG = 4
    !-----------------------------------------------------------------------
    ! Local variables:
    !-----------------------------------------------------------------------
    integer            :: ii, j, k        ! Loop indices
    integer            :: diag_len        ! Length of current diagonal(s)
    integer            :: ndiag           ! Number of current diagonals
    integer            :: offset(4)       ! Offsets of current diagonal
    integer            :: jad_len(jdiag)  ! Lengths of jagged diagonals

!CDIR NEIGHBORS
    jad_len(1:jdiag) = ia(2:jdiag+1) - ia(1:jdiag)

!print *, "amuxj_group:jdiag =", jdiag
!print *, "amuxj_group:jad_len:", jad_len(1:jdiag)

    y(1:n) = 0
    ii = 1
    outer: do
       if (ii > jdiag) exit outer
       diag_len = jad_len(ii)                   ! Length of next diagonal
       ndiag = 1
!NEC$ loop_count(3)
       do k = 1, min (jdiag-ii, MAXDIAG-1)      ! Check for further diagonals
          if (jad_len(ii+k) /= diag_len) exit
          ndiag = k+1
       end do
       offset(1:ndiag) = ia(ii:ii+ndiag-1) - 1  ! Offset of diagonals
!print *, "amuxj_group:", ii, diag_len, ndiag, offset(1:ndiag)
       select case (ndiag)
       case (1)          ! Single diagonal
          do j=1, diag_len
             y(j) = y(j) + a(offset(1)+j) * x(ja(offset(1)+j))
          end do
          ii = ii+1
       case (2)          ! Two diagonals of same length
          do j=1, diag_len
             y(j) = y(j) + a(offset(1)+j) * x(ja(offset(1)+j)) &
                         + a(offset(2)+j) * x(ja(offset(2)+j))
          end do
          ii = ii+2
       case (3)          ! Three diagonals of same length
          do j=1, diag_len
             y(j) = y(j) + a(offset(1)+j) * x(ja(offset(1)+j)) &
                         + a(offset(2)+j) * x(ja(offset(2)+j)) &
                         + a(offset(3)+j) * x(ja(offset(3)+j))
          end do
          ii = ii+3
       case default
!      case (4)          ! Four diagonals of same length
          do j=1, diag_len
             y(j) = y(j) + a(offset(1)+j) * x(ja(offset(1)+j)) &
                         + a(offset(2)+j) * x(ja(offset(2)+j)) &
                         + a(offset(3)+j) * x(ja(offset(3)+j)) &
                         + a(offset(4)+j) * x(ja(offset(4)+j))
          end do
          ii = ii+4
       end select
    end do outer
  end subroutine amuxj_grp_r8
  !-----------------------------------------------------------------------
  subroutine amuxj_grp_r4_r8 (n, x, y, jdiag, a, ja, ia)
    integer,  intent(in)      :: n        ! Row dimension of A
    real(wp), intent(in)      :: x(n)
    real(wp), intent(out)     :: y(n)
    integer,  intent(in)      :: jdiag    ! Number of jagged-diagonals
    real(sp), intent(in)      :: a(:)
    integer,  intent(in)      :: ja(:)
    integer,  intent(in)      :: ia(:)
    !-----------------------------------------------------------------------
    !        A times a vector in Jagged-Diagonal storage format (JAD)
    !-----------------------------------------------------------------------
    ! "Grouped" implementation to simultaneously handle up to 4 diagonals.
    ! This is supposed to give better performance on vector processors.
    !-----------------------------------------------------------------------
    ! Machine-dependent tuning parameter:
    ! maximum no. of simultaneous diagonals (range: 1-4)
    !
    integer, parameter :: MAXDIAG = 4
    !-----------------------------------------------------------------------
    ! Local variables:
    !-----------------------------------------------------------------------
    integer            :: ii, j, k        ! Loop indices
    integer            :: diag_len        ! Length of current diagonal(s)
    integer            :: ndiag           ! Number of current diagonals
    integer            :: offset(4)       ! Offsets of current diagonal
    integer            :: jad_len(jdiag)  ! Lengths of jagged diagonals

!CDIR NEIGHBORS
    jad_len(1:jdiag) = ia(2:jdiag+1) - ia(1:jdiag)

!print *, "amuxj_group:jdiag =", jdiag
!print *, "amuxj_group:jad_len:", jad_len(1:jdiag)

    y(1:n) = 0
    ii = 1
    outer: do
       if (ii > jdiag) exit outer
       diag_len = jad_len(ii)                   ! Length of next diagonal
       ndiag = 1
!NEC$ loop_count(3)
       do k = 1, min (jdiag-ii, MAXDIAG-1)      ! Check for further diagonals
          if (jad_len(ii+k) /= diag_len) exit
          ndiag = k+1
       end do
       offset(1:ndiag) = ia(ii:ii+ndiag-1) - 1  ! Offset of diagonals
!print *, "amuxj_group:", ii, diag_len, ndiag, offset(1:ndiag)
       select case (ndiag)
       case (1)          ! Single diagonal
          do j=1, diag_len
             y(j) = y(j) + a(offset(1)+j) * x(ja(offset(1)+j))
          end do
          ii = ii+1
       case (2)          ! Two diagonals of same length
          do j=1, diag_len
             y(j) = y(j) + a(offset(1)+j) * x(ja(offset(1)+j)) &
                         + a(offset(2)+j) * x(ja(offset(2)+j))
          end do
          ii = ii+2
       case (3)          ! Three diagonals of same length
          do j=1, diag_len
             y(j) = y(j) + a(offset(1)+j) * x(ja(offset(1)+j)) &
                         + a(offset(2)+j) * x(ja(offset(2)+j)) &
                         + a(offset(3)+j) * x(ja(offset(3)+j))
          end do
          ii = ii+3
       case default
!      case (4)          ! Four diagonals of same length
          do j=1, diag_len
             y(j) = y(j) + a(offset(1)+j) * x(ja(offset(1)+j)) &
                         + a(offset(2)+j) * x(ja(offset(2)+j)) &
                         + a(offset(3)+j) * x(ja(offset(3)+j)) &
                         + a(offset(4)+j) * x(ja(offset(4)+j))
          end do
          ii = ii+4
       end select
    end do outer
  end subroutine amuxj_grp_r4_r8
  !-----------------------------------------------------------------------
  subroutine amuxj_vec_r8 (n, x, y, jdiag, a, ja, ia)
    integer,  intent(in)      :: n        ! Row dimension of A
    real(wp), intent(in)      :: x(n)
    real(wp), intent(out)     :: y(n)
    integer,  intent(in)      :: jdiag    ! Number of jagged-diagonals
    real(wp), intent(in)      :: a(:)
    integer,  intent(in)      :: ja(:)
    integer,  intent(in)      :: ia(:)
    !-----------------------------------------------------------------------
    !        A times a vector in Jagged-Diagonal storage format (JAD)
    !-----------------------------------------------------------------------
    ! Optimized for NEC SX-9 using vector registers and
    ! Assignable Data Buffers (ADB) by
    ! Rudolf Fischer (NEC) and Jens-Olaf Beismann (NEC).
    !-----------------------------------------------------------------------
    ! Local variables:
    !-----------------------------------------------------------------------
    integer            :: ii, j           ! Loop indices
    integer            :: offset          ! Offsets of current diagonal
    integer            :: jad_len(jdiag)  ! Lengths of jagged diagonals
    integer            :: js              ! Start index for strip-mining
    integer            :: je              ! End   index for strip-mining
    integer            :: l               ! Length of current strip
    integer            :: je_(jdiag)      ! End index of strip-mined diagonal
    integer, parameter :: VL = 256        ! Vector register length
    real(wp)           :: v(VL)
!NEC$ vreg(v)

!CDIR NEIGHBORS
!cdir on_adb(jad_len)
    jad_len(1:jdiag) = ia(2:jdiag+1) - ia(1:jdiag)

    do js = 1, n, VL
       je = min (js + VL-1, n)
       l  = je - js + 1
       v  = 0._wp
!cdir on_adb(jad_len)
!cdir on_adb(je_)
       je_(1:jdiag) = min (je, jad_len(1:jdiag))
       do ii = 1, jdiag
          if (je_(ii) < js) goto 4711
          offset = ia(ii) - 1
!cdir on_adb(x)
!NEC$ shortloop
          do j = js, je_(ii)
             v(j-js+1) = v(j-js+1) + a(offset+j) * x(ja(offset+j))
!            y(j     ) = y(j     ) + a(offset+j) * x(ja(offset+j))
          end do
       end do
 4711  y(js:js+l-1) = v(1:l)
    end do
  end subroutine amuxj_vec_r8
  !-----------------------------------------------------------------------
  subroutine amuxj_vec_r4_r8 (n, x, y, jdiag, a, ja, ia)
    integer,  intent(in)      :: n        ! Row dimension of A
    real(wp), intent(in)      :: x(n)
    real(wp), intent(out)     :: y(n)
    integer,  intent(in)      :: jdiag    ! Number of jagged-diagonals
    real(sp), intent(in)      :: a(:)
    integer,  intent(in)      :: ja(:)
    integer,  intent(in)      :: ia(:)
    !-----------------------------------------------------------------------
    !        A times a vector in Jagged-Diagonal storage format (JAD)
    !-----------------------------------------------------------------------
    ! Optimized for NEC SX-9 using vector registers and
    ! Assignable Data Buffers (ADB) by
    ! Rudolf Fischer (NEC) and Jens-Olaf Beismann (NEC).
    !-----------------------------------------------------------------------
    ! Local variables:
    !-----------------------------------------------------------------------
    integer            :: ii, j           ! Loop indices
    integer            :: offset          ! Offsets of current diagonal
    integer            :: jad_len(jdiag)  ! Lengths of jagged diagonals
    integer            :: js              ! Start index for strip-mining
    integer            :: je              ! End   index for strip-mining
    integer            :: l               ! Length of current strip
    integer            :: je_(jdiag)      ! End index of strip-mined diagonal
    integer, parameter :: VL = 256        ! Vector register length
    real(wp)           :: v(VL)
!NEC$ vreg(v)

!CDIR NEIGHBORS
!cdir on_adb(jad_len)
    jad_len(1:jdiag) = ia(2:jdiag+1) - ia(1:jdiag)

    do js = 1, n, VL
       je = min (js + VL-1, n)
       l  = je - js + 1
       v  = 0._wp
!cdir on_adb(jad_len)
!cdir on_adb(je_)
       je_(1:jdiag) = min (je, jad_len(1:jdiag))
       do ii = 1, jdiag
          if (je_(ii) < js) goto 4711
          offset = ia(ii) - 1
!cdir on_adb(x)
!NEC$ shortloop
          do j = js, je_(ii)
             v(j-js+1) = v(j-js+1) + a(offset+j) * x(ja(offset+j))
!            y(j     ) = y(j     ) + a(offset+j) * x(ja(offset+j))
          end do
       end do
 4711  y(js:js+l-1) = v(1:l)
    end do
  end subroutine amuxj_vec_r4_r8
  !============================================================================

SUBROUTINE csrdns(nrow,ncol,a,ja,ia,dns,ndns,ierr)

INTEGER, INTENT(IN)                      :: nrow
INTEGER, INTENT(IN)                      :: ncol
REAL(WP), INTENT(IN)                     :: a(*)
INTEGER, INTENT(IN)                      :: ja(*)
INTEGER, INTENT(IN)                      :: ia(*)
INTEGER, INTENT(IN)                      :: ndns
REAL(WP), INTENT(OUT)                    :: dns(ndns,*)
INTEGER, INTENT(OUT)                     :: ierr


!-----------------------------------------------------------------------
! Compressed Sparse Row    to    Dense
!-----------------------------------------------------------------------
!
! converts a row-stored sparse matrix into a densely stored one
!
! On entry:
!----------
! nrow = row-dimension of a
! ncol = column dimension of a
! a,
! ja,
! ia    = input matrix in compressed sparse row format.
!         (a=value array, ja=column array, ia=pointer array)
! dns   = array where to store dense matrix
! ndns = first dimension of array dns
!
! on return:
!-----------
! dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*)
!
! ierr  = integer error indicator.
!         ierr .eq. 0  means normal return
!         ierr .eq. i  means that the code has stopped when processing
!         row number i, because it found a column number .gt. ncol.
!-----------------------------------------------------------------------
integer :: i, j, k

ierr = 0
DO  i=1, nrow
  DO  j=1,ncol
    dns(i,j) = 0.0D0
  END DO
END DO

DO  i=1,nrow
  DO  k=ia(i),ia(i+1)-1
    j = ja(k)
    IF (j > ncol) THEN
      ierr = i
      RETURN
    END IF
    dns(i,j) = a(k)
  END DO
END DO
RETURN
!---- end of csrdns ----------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrdns
!-----------------------------------------------------------------------

SUBROUTINE dnscsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)

INTEGER, INTENT(IN)                      :: nrow
INTEGER, INTENT(IN)                      :: ncol
INTEGER, INTENT(IN)                      :: nzmax
INTEGER, INTENT(IN)                      :: ndns
REAL(WP), INTENT(IN)                     :: dns(ndns,*)
REAL(WP), INTENT(OUT)                    :: a(*)
INTEGER, INTENT(OUT)                     :: ja(*)
INTEGER, INTENT(OUT)                     :: ia(*)
INTEGER, INTENT(OUT)                     :: ierr

!-----------------------------------------------------------------------
! Dense  to    Compressed Row Sparse
!-----------------------------------------------------------------------

! converts a densely stored matrix into a row orientied
! compactly sparse matrix. ( reverse of csrdns )
! Note: this routine does not check whether an element
! is small. It considers that a(i,j) is zero if it is exactly
! equal to zero: see test below.
!-----------------------------------------------------------------------
! on entry:
!---------

! nrow = row-dimension of a
! ncol = column dimension of a
! nzmax = maximum number of nonzero elements allowed. This
!         should be set to be the lengths of the arrays a and ja.
! dns   = input nrow x ncol (dense) matrix.
! ndns = first dimension of dns.

! on return:
!----------

! a, ja, ia = value, column, pointer  arrays for output matrix

! ierr = integer error indicator:
!         ierr .eq. 0 means normal retur
!         ierr .eq. i means that the the code stopped while
!         processing row number i, because there was no space left in
!         a, and ja (as defined by parameter nzmax).
!-----------------------------------------------------------------------
integer :: i, j, next

ierr = 0
next = 1
ia(1) = 1
DO  i=1,nrow
  DO  j=1, ncol
    IF (dns(i,j) == 0.0D0) CYCLE
    IF (next > nzmax) THEN
      ierr = i
      RETURN
    END IF
    ja(next) = j
    a(next) = dns(i,j)
    next = next+1
  END DO
  ia(i+1) = next
END DO
RETURN
!---- end of dnscsr ----------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE dnscsr
!-----------------------------------------------------------------------

SUBROUTINE coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)                      :: nrow
INTEGER, INTENT(IN)                      :: nnz
REAL(WP), INTENT(IN)                       :: a(*)
INTEGER, INTENT(IN)                      :: ir(*)
INTEGER, INTENT(IN)                      :: jc(*)
REAL(WP), INTENT(OUT)                      :: ao(*)
INTEGER, INTENT(OUT)                     :: jao(*)
INTEGER, INTENT(OUT)                     :: iao(*)
REAL(WP)  x

!-----------------------------------------------------------------------
!  Coordinate     to   Compressed Sparse Row
!-----------------------------------------------------------------------
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.

! on entry:
!---------
! nrow = dimension of the matrix
! nnz = number of nonzero elements in matrix
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
!    the elements, ir(k) = its row number and jc(k) = its column
!   number. The order of the elements is arbitrary.

! on return:
!-----------
! ir  is destroyed

! ao, jao, iao = matrix in general sparse matrix format with ao
!  continung the real values, jao containing the column indices,
! and iao being the pointer to the beginning of the row,
! in arrays ao, jao.

! Notes:
!------ This routine is NOT in place.  See coicsr

!------------------------------------------------------------------------
integer :: i, j, k, k0, iad

DO  k=1,nrow+1
  iao(k) = 0
END DO
! determine row-lengths.
DO  k=1, nnz
  iao(ir(k)) = iao(ir(k))+1
END DO
! starting position of each row..
k = 1
DO  j=1,nrow+1
  k0 = iao(j)
  iao(j) = k
  k = k+k0
END DO
! go through the structure  once more. Fill in output matrix.
DO  k=1, nnz
  i = ir(k)
  j = jc(k)
  x = a(k)
  iad = iao(i)
  ao(iad) =  x
  jao(iad) = j
  iao(i) = iad+1
END DO
! shift back iao
DO  j=nrow,1,-1
  iao(j+1) = iao(j)
END DO
iao(1) = 1
RETURN
!------------- end of coocsr -------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE coocsr
!-----------------------------------------------------------------------

SUBROUTINE coicsr (n,nnz,job,a,ja,ia,iwk)

INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: nnz
INTEGER, INTENT(IN OUT)                  :: job
REAL(WP), INTENT(IN OUT)                   :: a(*)
INTEGER, INTENT(IN OUT)                  :: ja(nnz)
INTEGER, INTENT(IN OUT)                  :: ia(nnz)
INTEGER, INTENT(OUT)                     :: iwk(n+1)


!------------------------------------------------------------------------
! IN-PLACE coo-csr conversion routine.
!------------------------------------------------------------------------
! this subroutine converts a matrix stored in coordinate format into
! the csr format. The conversion is done in place in that the arrays
! a,ja,ia of the result are overwritten onto the original arrays.
!------------------------------------------------------------------------
! on entry:
!---------
! n = integer. row dimension of A.
! nnz = integer. number of nonzero elements in A.
! job   = integer. Job indicator. when job=1, the real values in a are
!         filled. Otherwise a is not touched and the structure of the
!         array only (i.e. ja, ia)  is obtained.
! a = real array of size nnz (number of nonzero elements in A)
!         containing the nonzero elements
! ja = integer array of length nnz containing the column positions
!    of the corresponding elements in a.
! ia = integer array of length nnz containing the row positions
!    of the corresponding elements in a.
! iwk = integer work array of length n+1
! on return:
!----------
! a
! ja
! ia = contains the compressed sparse row data structure for the
!         resulting matrix.
! Note:
!-------
!         the entries of the output matrix are not sorted (the column
!         indices in each are not in increasing order) use coocsr
!         if you want them sorted.
!----------------------------------------------------------------------c
!  Coded by Y. Saad, Sep. 26 1989                                      c
!----------------------------------------------------------------------c
integer  :: i, j, k, init, ipos, inext, jnext
REAL(WP) :: t,tnext
LOGICAL  :: values
!-----------------------------------------------------------------------
values = (job == 1)
! find pointer array for resulting matrix.
DO  i=1,n+1
   iwk(i) = 0
END DO
DO  k=1,nnz
   i = ia(k)
   iwk(i+1) = iwk(i+1)+1
END DO
!------------------------------------------------------------------------
iwk(1) = 1
DO  i=2,n
   iwk(i) = iwk(i-1) + iwk(i)
END DO

!     loop for a cycle in chasing process.

init = 1
k = 0
do
   IF (values) t = a(init)
   i = ia(init)
   j = ja(init)
   ia(init) = -1
   !------------------------------------------------------------------------
6  k = k+1
   !     current row number is i.  determine  where to go.
   ipos = iwk(i)
   !     save the chased element.
   IF (values) tnext = a(ipos)
   inext = ia(ipos)
   jnext = ja(ipos)
   !     then occupy its location.
   IF (values) a(ipos)  = t
   ja(ipos) = j
   !     update pointer information for next element to come in row i.
   iwk(i) = ipos+1
   !     determine  next element to be chased,
   IF (ia(ipos) < 0) GO TO 65
   t = tnext
   i = inext
   j = jnext
   ia(ipos) = -1
   IF (k < nnz) GO TO 6
   exit
65 init = init+1
   IF (init > nnz) exit
   IF (ia(init) < 0) GO TO 65
   !     restart chasing --
end do
DO  i=1,n
   ia(i+1) = iwk(i)
END DO
ia(1) = 1
RETURN
!----------------- end of coicsr ----------------------------------------
!------------------------------------------------------------------------
END SUBROUTINE coicsr
!-----------------------------------------------------------------------
SUBROUTINE csrcoo (nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)
!-----------------------------------------------------------------------
INTEGER, INTENT(IN)                      :: nrow
INTEGER, INTENT(IN OUT)                  :: job
INTEGER, INTENT(IN)                      :: nzmax
REAL(WP), INTENT(IN)                       :: a(*)
INTEGER, INTENT(IN)                      :: ja(*)
INTEGER, INTENT(IN)                      :: ia(nrow+1)
INTEGER, INTENT(OUT)                     :: nnz
REAL(WP), INTENT(OUT)                      :: ao(*)
INTEGER, INTENT(OUT)                     :: ir(*)
INTEGER, INTENT(OUT)                     :: jc(*)
INTEGER, INTENT(OUT)                     :: ierr

!-----------------------------------------------------------------------
!  Compressed Sparse Row      to      Coordinate
!-----------------------------------------------------------------------
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.

! on entry:
!---------
! nrow = dimension of the matrix.
! job   = integer serving as a job indicator.
!         if job = 1 fill in only the array ir, ignore jc, and ao.
!         if job = 2 fill in ir, and jc but not ao
!         if job = 3 fill in everything.
!         The reason why these options are provided is that on return
!         ao and jc are the same as a, ja. So when job = 3, a and ja are
!         simply copied into ao, jc.  When job=2, only jc and ir are
!         returned. With job=1 only the array ir is returned. Moreover,
!         the algorithm is in place:
!      call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr)
!         will write the output matrix in coordinate format on a, ja,ia.

! a,
! ja,
! ia    = matrix in compressed sparse row format.
! nzmax = length of space available in ao, ir, jc.
!         the code will stop immediatly if the number of
!         nonzero elements found in input matrix exceeds nzmax.

! on return:
!-----------
! ao, ir, jc = matrix in coordinate format.

! nnz        = number of nonzero elements in matrix.
! ierr       = integer error indicator.
!         ierr .eq. 0 means normal retur
!         ierr .eq. 1 means that the the code stopped
!         because there was no space in ao, ir, jc
!         (according to the value of  nzmax).

! NOTES: 1)This routine is PARTIALLY in place: csrcoo can be called with
!         ao being the same array as as a, and jc the same array as ja.
!         but ir CANNOT be the same as ia.
!         2) note the order in the output arrays,
!------------------------------------------------------------------------
integer :: i, k, k1, k2

ierr = 0
nnz = ia(nrow+1)-1
IF (nnz > nzmax) THEN
  ierr = 1
  RETURN
END IF
!------------------------------------------------------------------------
SELECT CASE ( job )
  CASE (    1)
    GO TO 3
  CASE (    2)
    GO TO 2
  CASE (    3)
    GO TO 1
END SELECT
1    DO  k=1,nnz
  ao(k) = a(k)
END DO
2    DO  k=1,nnz
  jc(k) = ja(k)
END DO

!     copy backward to allow for in-place processing.

3    DO  i=nrow,1,-1
  k1 = ia(i+1)-1
  k2 = ia(i)
  DO  k=k1,k2,-1
    ir(k) = i
  END DO
END DO
RETURN
!------------- end-of-csrcoo -------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrcoo
!-----------------------------------------------------------------------

SUBROUTINE csrssr (nrow,a,ja,ia,nzmax,ao,jao,iao,ierr)

INTEGER, INTENT(IN)                      :: nrow
REAL(WP), INTENT(IN)                       :: a(*)
INTEGER, INTENT(IN)                      :: ja(*)
INTEGER, INTENT(IN)                      :: ia(*)
INTEGER, INTENT(IN)                      :: nzmax
REAL(WP), INTENT(OUT)                      :: ao(*)
INTEGER, INTENT(OUT)                     :: jao(*)
INTEGER, INTENT(OUT)                     :: iao(*)
INTEGER, INTENT(OUT)                     :: ierr
REAL(WP)  t

!-----------------------------------------------------------------------
! Compressed Sparse Row     to     Symmetric Sparse Row
!-----------------------------------------------------------------------
! this subroutine extracts the lower triangular part of a matrix.
! It can used as a means for converting a symmetric matrix for
! which all the entries are stored in sparse format into one
! in which only the lower part is stored. The routine is in place in
! that the output matrix ao, jao, iao can be overwritten on
! the  input matrix  a, ja, ia if desired. Csrssr has been coded to
! put the diagonal elements of the matrix in the last position in
! each row (i.e. in position  ao(ia(i+1)-1   of ao and jao)
!-----------------------------------------------------------------------
! On entry
!-----------
! nrow  = dimension of the matrix a.
! a, ja,
!    ia = matrix stored in compressed row sparse format

! nzmax = length of arrays ao,  and jao.

! On return:
!-----------
! ao, jao,
!     iao = lower part of input matrix (a,ja,ia) stored in compressed sparse
!          row format format.

! ierr   = integer error indicator.
!          ierr .eq. 0  means normal return
!          ierr .eq. i  means that the code has stopped when processing
!          row number i, because there is not enough space in ao, jao
!          (according to the value of nzmax)

!-----------------------------------------------------------------------
integer :: i, k, ko, kold, kdiag

ierr = 0
ko = 0
!-----------------------------------------------------------------------
DO   i=1, nrow
  kold = ko
  kdiag = 0
  DO  k = ia(i), ia(i+1) -1
    IF (ja(k)  > i) CYCLE
    ko = ko+1
    IF (ko > nzmax) THEN
      ierr = i
      RETURN
    END IF
    ao(ko) = a(k)
    jao(ko) = ja(k)
    IF (ja(k)  == i) kdiag = ko
  END DO
  IF (kdiag == 0 .OR. kdiag == ko) GO TO 72

!     exchange

  t = ao(kdiag)
  ao(kdiag) = ao(ko)
  ao(ko) = t

  k = jao(kdiag)
  jao(kdiag) = jao(ko)
  jao(ko) = k
  72      iao(i) = kold+1
END DO
!     redefine iao(n+1)
iao(nrow+1) = ko+1
RETURN
!--------- end of csrssr -----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrssr
!-----------------------------------------------------------------------

SUBROUTINE ssrcsr(job, value2, nrow, a, ja, ia, nzmax,  &
    ao, jao, iao, indu, iwk, ierr)
!     .. Scalar Arguments ..

INTEGER, INTENT(IN)                      :: job
INTEGER, INTENT(IN)                      :: value2
INTEGER, INTENT(IN)                      :: nrow
REAL(WP), INTENT(IN)                       :: a(*)
INTEGER, INTENT(IN)                      :: ja(*)
INTEGER, INTENT(IN)                      :: ia(nrow+1)
INTEGER, INTENT(IN)                      :: nzmax
REAL(WP), INTENT(OUT)                      :: ao(nzmax)
INTEGER, INTENT(OUT)                     :: jao(nzmax)
INTEGER, INTENT(OUT)                     :: iao(nrow+1)
INTEGER, INTENT(OUT)                     :: indu(nrow)
INTEGER, INTENT(OUT)                     :: iwk(nrow+1)
INTEGER, INTENT(OUT)                     :: ierr

!     ..
!     .. Array Arguments ..


!     ..
!-----------------------------------------------------------------------
!     Symmetric Sparse Row to Compressed Sparse Row format
!-----------------------------------------------------------------------
!     This subroutine converts a given matrix in SSR format to regular
!     CSR format by computing Ao = A + A' - diag(A), where A' is A
!     transpose.

!     Typically this routine is used to expand the SSR matrix of
!     Harwell Boeing matrices, or to obtain a symmetrized graph of
!     unsymmetric matrices.

!     This routine is inplace, i.e., (Ao,jao,iao) may be same as
!     (a,ja,ia).

!     It is possible to input an arbitrary CSR matrix to this routine,
!     since there is no syntactical difference between CSR and SSR
!     format. It also removes duplicate entries and perform a partial
!     ordering. The output matrix has an order of lower half, main
!     diagonal and upper half after the partial ordering.
!-----------------------------------------------------------------------
! on entry:
!---------

! job   = options
!         0 -- duplicate entries are not removed. If the input matrix is
!              SSR (not an arbitary CSR) matrix, no duplicate entry should
!              arise from this routine.
!         1 -- eliminate duplicate entries, zero entries.
!         2 -- eliminate duplicate entries and perform partial ordering.
!         3 -- eliminate duplicate entries, sort the entries in the
!              increasing order of clumn indices.

! value2= will the values of A be copied?
!         0 -- only expand the graph (a, ao are not touched)
!         1 -- expand the matrix with the values.

! nrow  = column dimension of inout matrix
! a,
! ia,
! ja    = matrix in compressed sparse row format.

! nzmax = size of arrays ao and jao. SSRCSR will abort if the storage
!          provided in ao, jao is not sufficient to store A. See ierr.

! on return:
!----------
! ao, jao, iao
!       = output matrix in compressed sparse row format. The resulting
!         matrix is symmetric and is equal to A+A'-D. ao, jao, iao,
!         can be the same as a, ja, ia in the calling sequence.

! indu  = integer array of length nrow. INDU will contain pointers
!         to the beginning of upper traigular part if job > 1.
!         Otherwise it is also used as a work array (size nrow).

! iwk   = integer work space (size nrow+1).

! ierr  = integer. Serving as error message. If the length of the arrays
!         ao, jao exceeds nzmax, ierr returns the minimum value
!         needed for nzmax. otherwise ierr=0 (normal return).

!-----------------------------------------------------------------------
!     .. Local Scalars ..
INTEGER :: i, ipos, j, k, kfirst, klast, ko, kosav, nnz
REAL(WP)             tmp
!     ..
!     .. Executable Statements ..
ierr = 0
DO  i = 1, nrow
  indu(i) = 0
  iwk(i) = 0
END DO
iwk(nrow+1) = 0

!     .. compute number of elements in each row of (A'-D)
!     put result in iwk(i+1)  for row i.

DO  i = 1, nrow
  DO  k = ia(i), ia(i+1) - 1
    j = ja(k)
    IF (j /= i) iwk(j+1) = iwk(j+1) + 1
  END DO
END DO

!     .. find addresses of first elements of ouput matrix. result in iwk

iwk(1) = 1
DO  i = 1, nrow
  indu(i) = iwk(i) + ia(i+1) - ia(i)
  iwk(i+1) = iwk(i+1) + indu(i)
  indu(i) = indu(i) - 1
END DO
!.....Have we been given enough storage in ao, jao ?
nnz = iwk(nrow+1) - 1
IF (nnz > nzmax) THEN
  ierr = nnz
  RETURN
END IF

!     .. copy the existing matrix (backwards).

kosav = iwk(nrow+1)
DO  i = nrow, 1, -1
  klast = ia(i+1) - 1
  kfirst = ia(i)
  iao(i+1) = kosav
  kosav = iwk(i)
  ko = iwk(i) - kfirst
  iwk(i) = ko + klast + 1
  DO  k = klast, kfirst, -1
    IF (value2 /= 0) ao(k+ko) = a(k)
    jao(k+ko) = ja(k)
  END DO
END DO
iao(1) = 1

!     now copy (A'-D). Go through the structure of ao, jao, iao
!     that has already been copied. iwk(i) is the address
!     of the next free location in row i for ao, jao.

DO  i = 1, nrow
  DO  k = iao(i), indu(i)
    j = jao(k)
    IF (j /= i) THEN
      ipos = iwk(j)
      IF (value2 /= 0) ao(ipos) = ao(k)
      jao(ipos) = i
      iwk(j) = ipos + 1
    END IF
  END DO
END DO
IF (job <= 0) RETURN

!     .. eliminate duplicate entries --
!     array INDU is used as marker for existing indices, it is also the
!     location of the entry.
!     IWK is used to stored the old IAO array.
!     matrix is copied to squeeze out the space taken by the duplicated
!     entries.

DO  i = 1, nrow
  indu(i) = 0
  iwk(i) = iao(i)
END DO
iwk(nrow+1) = iao(nrow+1)
k = 1
DO  i = 1, nrow
  iao(i) = k
  ipos = iwk(i)
  klast = iwk(i+1)
  100     IF (ipos < klast) THEN
    j = jao(ipos)
    IF (indu(j) == 0) THEN
!     .. new entry ..
      IF (value2 /= 0) THEN
        IF (ao(ipos) /= 0.0D0) THEN
          indu(j) = k
          jao(k) = jao(ipos)
          ao(k) = ao(ipos)
          k = k + 1
        END IF
      ELSE
        indu(j) = k
        jao(k) = jao(ipos)
        k = k + 1
      END IF
    ELSE IF (value2 /= 0) THEN
!     .. duplicate entry ..
      ao(indu(j)) = ao(indu(j)) + ao(ipos)
    END IF
    ipos = ipos + 1
    GO TO 100
  END IF
!     .. remove marks before working on the next row ..
  DO  ipos = iao(i), k - 1
    indu(jao(ipos)) = 0
  END DO
END DO
iao(nrow+1) = k
IF (job <= 1) RETURN

!     .. partial ordering ..
!     split the matrix into strict upper/lower triangular
!     parts, INDU points to the the beginning of the strict upper part.

DO  i = 1, nrow
  klast = iao(i+1) - 1
  kfirst = iao(i)
  130     IF (klast > kfirst) THEN
    IF (jao(klast) < i .AND. jao(kfirst) >= i) THEN
!     .. swap klast with kfirst ..
      j = jao(klast)
      jao(klast) = jao(kfirst)
      jao(kfirst) = j
      IF (value2 /= 0) THEN
        tmp = ao(klast)
        ao(klast) = ao(kfirst)
        ao(kfirst) = tmp
      END IF
    END IF
    IF (jao(klast) >= i) klast = klast - 1
    IF (jao(kfirst) < i) kfirst = kfirst + 1
    GO TO 130
  END IF

  IF (jao(klast) < i) THEN
    indu(i) = klast + 1
  ELSE
    indu(i) = klast
  END IF
END DO
IF (job <= 2) RETURN

!     .. order the entries according to column indices
!     bubble-sort is used

DO  i = 1, nrow
  DO  ipos = iao(i), indu(i)-1
    DO  j = indu(i)-1, ipos+1, -1
      k = j - 1
      IF (jao(k) > jao(j)) THEN
        ko = jao(k)
        jao(k) = jao(j)
        jao(j) = ko
        IF (value2 /= 0) THEN
          tmp = ao(k)
          ao(k) = ao(j)
          ao(j) = tmp
        END IF
      END IF
    END DO
  END DO
  DO  ipos = indu(i), iao(i+1)-1
    DO  j = iao(i+1)-1, ipos+1, -1
      k = j - 1
      IF (jao(k) > jao(j)) THEN
        ko = jao(k)
        jao(k) = jao(j)
        jao(j) = ko
        IF (value2 /= 0) THEN
          tmp = ao(k)
          ao(k) = ao(j)
          ao(j) = tmp
        END IF
      END IF
    END DO
  END DO
END DO

RETURN
!---- end of ssrcsr ----------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE ssrcsr
!-----------------------------------------------------------------------

SUBROUTINE xssrcsr (nrow,a,ja,ia,nzmax,ao,jao,iao,indu,ierr)

INTEGER, INTENT(IN)                      :: nrow
REAL(WP), INTENT(IN)                       :: a(*)
INTEGER, INTENT(IN)                      :: ja(*)
INTEGER, INTENT(IN)                      :: ia(nrow+1)
INTEGER, INTENT(IN)                      :: nzmax
REAL(WP), INTENT(OUT)                      :: ao(nzmax)
INTEGER, INTENT(OUT)                     :: jao(nzmax)
INTEGER, INTENT(OUT)                     :: iao(nrow+1)
INTEGER, INTENT(OUT)                     :: indu(nrow+1)
INTEGER, INTENT(OUT)                     :: ierr


!-----------------------------------------------------------------------
! Symmetric Sparse Row   to    (regular) Compressed Sparse Row
!-----------------------------------------------------------------------
! this subroutine converts  a symmetric  matrix in which only the lower
! part is  stored in compressed sparse row format, i.e.,
! a matrix stored in symmetric sparse format, into a fully stored matrix
! i.e., a matrix where both the lower and upper parts are stored in
! compressed sparse row format. the algorithm is in place (i.e. result
! may be overwritten onto the input matrix a, ja, ia ----- ).
! the output matrix delivered by ssrcsr is such that each row starts with
! the elements of the lower part followed by those of the upper part.
!-----------------------------------------------------------------------
! on entry:
!---------
!
! nrow  = row dimension of inout matrix
! a,
! ia,
! ja    = matrix in compressed sparse row format. This is assumed to be
!         a lower triangular matrix.

! nzmax = size of arrays ao and jao. ssrcsr will abort if the storage
!    provided in a, ja is not sufficient to store A. See ierr.
!
! on return:
!----------
! ao, iao,
!   jao = output matrix in compressed sparse row format. The resulting
!         matrix is symmetric and is equal to A+A**T - D, if
!         A is the original lower triangular matrix. ao, jao, iao,
!         can be the same as a, ja, ia in the calling sequence.

! indu  = integer array of length nrow+1. If the input matrix is such
!         that the last element in each row is its diagonal element then
!         on return, indu will contain the pointers to the diagonal
!         element in each row of the output matrix. Otherwise used as
!         work array.
! ierr  = integer. Serving as error message. If the length of the arrays
!         ao, jao exceeds nzmax, ierr returns the minimum value
!         needed for nzmax. otherwise ierr=0 (normal return).

!-----------------------------------------------------------------------
integer :: i, j, k, lenrow, nnz, kosav, klast, kfirst, ko, ipos

ierr = 0
DO  i=1,nrow+1
  indu(i) = 0
END DO

!     compute  number of elements in each row of strict upper part.
!     put result in indu(i+1)  for row i.

DO  i=1, nrow
  DO  k=ia(i),ia(i+1)-1
    j = ja(k)
    IF (j < i) indu(j+1) = indu(j+1)+1
  END DO
END DO
!-----------
!     find addresses of first elements of ouput matrix. result in indu
!-----------
indu(1) = 1
DO  i=1,nrow
  lenrow = ia(i+1)-ia(i)
  indu(i+1) = indu(i) + indu(i+1) + lenrow
END DO
!--------------------- enough storage in a, ja ? --------
nnz = indu(nrow+1)-1
IF (nnz > nzmax) THEN
  ierr = nnz
  RETURN
END IF

!     now copy lower part (backwards).

kosav = indu(nrow+1)
DO  i=nrow,1,-1
  klast = ia(i+1)-1
  kfirst = ia(i)
  iao(i+1) = kosav
  ko = indu(i)
  kosav = ko
  DO  k = kfirst, klast
    ao(ko) = a(k)
    jao(ko) = ja(k)
    ko = ko+1
  END DO
  indu(i) = ko
END DO
iao(1) = 1

!     now copy upper part. Go through the structure of ao, jao, iao
!     that has already been copied (lower part). indu(i) is the address
!     of the next free location in row i for ao, jao.

loop8:  DO  i=1,nrow
!     i-th row is now in ao, jao, iao structure -- lower half part
  DO  k=iao(i), iao(i+1)-1
    j = jao(k)
    IF (j >= i)  CYCLE loop8
    ipos = indu(j)
    ao(ipos) = ao(k)
    jao(ipos) = i
    indu(j) = indu(j) + 1
  END DO
END DO loop8
RETURN
!----- end of xssrcsr --------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE xssrcsr
!-----------------------------------------------------------------------

SUBROUTINE csrmsr (n,a,ja,ia,ao,jao,wk,iwk)

INTEGER, INTENT(IN)                      :: n
REAL(WP), INTENT(IN)                       :: a(*)
INTEGER, INTENT(IN)                      :: ja(*)
INTEGER, INTENT(IN)                      :: ia(n+1)
REAL(WP), INTENT(OUT)                      :: ao(*)
INTEGER, INTENT(OUT)                     :: jao(*)
REAL(WP), INTENT(OUT)                      :: wk(n)
INTEGER, INTENT(OUT)                     :: iwk(n+1)


!-----------------------------------------------------------------------
! Compressed Sparse Row   to      Modified - Sparse Row
!                                 Sparse row with separate main diagonal
!-----------------------------------------------------------------------
! converts a general sparse matrix a, ja, ia into
! a compressed matrix using a separated diagonal (referred to as
! the bell-labs format as it is used by bell labs semi conductor
! group. We refer to it here as the modified sparse row format.
! Note: this has been coded in such a way that one can overwrite
! the output matrix onto the input matrix if desired by a call of
! the form

!     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)

! In case ao, jao, are different from a, ja, then one can
! use ao, jao as the work arrays in the calling sequence:

!     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)

!-----------------------------------------------------------------------

! on entry :
!---------
! a, ja, ia = matrix in csr format. note that the
!      algorithm is in place: ao, jao can be the same
!            as a, ja, in which case it will be overwritten on it
!            upon return.
!
! on return :
!-----------

! ao, jao  = sparse matrix in modified sparse row storage format:
!    +  ao(1:n) contains the diagonal of the matrix.
!    +  ao(n+2:nnz) contains the nondiagonal elements of the
!             matrix, stored rowwise.
!    +  jao(n+2:nnz) : their column indices
!    +  jao(1:n+1) contains the pointer array for the nondiagonal
!             elements in ao(n+1:nnz) and jao(n+2:nnz).
!             i.e., for i .le. n+1 jao(i) points to beginning of row i
!       in arrays ao, jao.
!        here nnz = number of nonzero elements+1
! work arrays:
!------------
! wk = real work array of length n
! iwk   = integer work array of length n+1

! notes:
!-------
!        Algorithm is in place.  i.e. both:

!          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
!          (in which  ao, jao, are different from a, ja)
!           and
!          call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
!          (in which  wk, jwk, are different from a, ja)
!        are OK.
!--------
! coded by Y. Saad Sep. 1989. Rechecked Feb 27, 1990.
!-----------------------------------------------------------------------
integer :: i, j, k, icount, iptr, ii

icount = 0

! store away diagonal elements and count nonzero diagonal elements.

DO  i=1,n
  wk(i) = 0.0D0
  iwk(i+1) = ia(i+1)-ia(i)
  DO  k=ia(i),ia(i+1)-1
    IF (ja(k) == i) THEN
      wk(i) = a(k)
      icount = icount + 1
      iwk(i+1) = iwk(i+1)-1
    END IF
  END DO
END DO

! compute total length

iptr = n + ia(n+1) - icount

!     copy backwards (to avoid collisions)

DO  ii=n,1,-1
  DO  k=ia(ii+1)-1,ia(ii),-1
    j = ja(k)
    IF (j /= ii) THEN
      ao(iptr) = a(k)
      jao(iptr) = j
      iptr = iptr-1
    END IF
  END DO
END DO

! compute pointer values and copy wk(*)

jao(1) = n+2
DO  i=1,n
  ao(i) = wk(i)
  jao(i+1) = jao(i)+iwk(i+1)
END DO
RETURN
!------------ end of subroutine csrmsr ---------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrmsr
!-----------------------------------------------------------------------

SUBROUTINE msrcsr (n,a,ja,ao,jao,iao,wk,iwk)

INTEGER, INTENT(IN)                      :: n
REAL(WP), INTENT(IN)                       :: a(*)
INTEGER, INTENT(IN)                      :: ja(*)
REAL(WP), INTENT(OUT)                      :: ao(*)
INTEGER, INTENT(OUT)                     :: jao(*)
INTEGER, INTENT(OUT)                     :: iao(n+1)
REAL(WP), INTENT(OUT)                      :: wk(n)
INTEGER, INTENT(OUT)                     :: iwk(n+1)


!-----------------------------------------------------------------------
!       Modified - Sparse Row  to   Compressed Sparse Row

!-----------------------------------------------------------------------
! converts a compressed matrix using a separated diagonal
! (modified sparse row format) in the Compressed Sparse Row
! format.
! does not check for zero elements in the diagonal.


! on entry :
!---------
! n          = row dimension of matrix
! a, ja      = sparse matrix in msr sparse storage format
!              see routine csrmsr for details on data structure

! on return :
!-----------

! ao,jao,iao = output matrix in csr format.

! work arrays:
!------------
! wk       = real work array of length n
! iwk      = integer work array of length n+1

! notes:
!   The original version of this was NOT in place, but has
!   been modified by adding the vector iwk to be in place.
!   The original version had ja instead of iwk everywhere in
!   loop 500.  Modified  Sun 29 May 1994 by R. Bramley (Indiana).

!-----------------------------------------------------------------------
integer :: i, j, k, ii, iptr, idiag
LOGICAL :: added

DO  i=1,n
  wk(i) = a(i)
  iwk(i) = ja(i)
END DO
iwk(n+1) = ja(n+1)
iao(1) = 1
iptr = 1
!---------
DO  ii=1,n
  added = .false.
  idiag = iptr + (iwk(ii+1)-iwk(ii))
  DO  k=iwk(ii),iwk(ii+1)-1
    j = ja(k)
    IF (j < ii) THEN
      ao(iptr) = a(k)
      jao(iptr) = j
      iptr = iptr+1
    ELSE IF (added) THEN
      ao(iptr) = a(k)
      jao(iptr) = j
      iptr = iptr+1
    ELSE
! add diag element - only reserve a position for it.
      idiag = iptr
      iptr = iptr+1
      added = .true.
!     then other element
      ao(iptr) = a(k)
      jao(iptr) = j
      iptr = iptr+1
    END IF
  END DO
  ao(idiag) = wk(ii)
  jao(idiag) = ii
  IF (.NOT. added) iptr = iptr+1
  iao(ii+1) = iptr
END DO
RETURN
!------------ end of subroutine msrcsr ---------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE msrcsr
!-----------------------------------------------------------------------
subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
  integer,  intent(in)   :: n
  integer,  intent(in)   :: job
  integer,  intent(in)   :: ipos
  real(wp), intent(in)   :: a(:)
  integer,  intent(in)   :: ja(:)
  integer,  intent(in)   :: ia(n+1)
  real(wp), intent(out)  :: ao(:)
  integer,  intent(out)  :: jao(:)
  integer,  intent(out)  :: iao(n+1)

!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place.
!-----------------------------------------------------------------------
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n = dimension of A.
! job = integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2))
!   for any other normal usage, enter ipos=1.
! a = real array of length nnz (nnz=number of nonzero elements in input
!         matrix) containing the nonzero elements.
! ja = integer array of length nnz containing the column positions
!    of the corresponding elements in a.
! ia = integer of size n+1. ia(k) contains the position in a, ja of
!   the beginning of the k-th row.
!
! on return:
! ----------
! output arguments:
! ao = real array of size nzz containing the "a" part of the transpose
! jao = integer array of size nnz containing the column indices.
! iao = integer array of size n+1 containing the "ia" index array of
!   the transpose.
!-----------------------------------------------------------------------
  call csrcsc2 (n,n,job,ipos,a,ja,ia,ao,jao,iao)
end subroutine csrcsc
!-----------------------------------------------------------------------
subroutine csrcsc_r4 (n,job,ipos,a,ja,ia,ao,jao,iao)
  integer,  intent(in)   :: n
  integer,  intent(in)   :: job
  integer,  intent(in)   :: ipos
  real(sp), intent(in)   :: a(:)
  integer,  intent(in)   :: ja(:)
  integer,  intent(in)   :: ia(n+1)
  real(sp), intent(out)  :: ao(:)
  integer,  intent(out)  :: jao(:)
  integer,  intent(out)  :: iao(n+1)
  call csrcsc2_r4 (n,n,job,ipos,a,ja,ia,ao,jao,iao)
end subroutine csrcsc_r4
!-----------------------------------------------------------------------
subroutine csrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
  integer,  intent(in)   :: n
  integer,  intent(in)   :: n2
  integer,  intent(in)   :: job
  integer,  intent(in)   :: ipos
  real(wp), intent(in)   :: a(:)
  integer,  intent(in)   :: ja(:)
  integer,  intent(in)   :: ia(n+1)
  real(wp), intent(out)  :: ao(:)
  integer,  intent(out)  :: jao(:)
  integer,  intent(out)  :: iao(n2+1)

!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place.
!-----------------------------------------------------------------------
! Rectangular version.  n is number of rows of CSR matrix,
!                       n2 (input) is number of columns of CSC matrix.
!-----------------------------------------------------------------------
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n = number of rows of CSR matrix.
! n2    = number of columns of CSC matrix.
! job = integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2))
!   for any other normal usage, enter ipos=1.
! a = real array of length nnz (nnz=number of nonzero elements in input
!         matrix) containing the nonzero elements.
! ja = integer array of length nnz containing the column positions
!    of the corresponding elements in a.
! ia = integer of size n+1. ia(k) contains the position in a, ja of
!   the beginning of the k-th row.
!
! on return:
! ----------
! output arguments:
! ao = real array of size nzz containing the "a" part of the transpose
! jao = integer array of size nnz containing the column indices.
! iao = integer array of size n+1 containing the "ia" index array of
!   the transpose.
!-----------------------------------------------------------------------
  integer :: i, j, k, next
!----------------- compute lengths of rows of transp(A) ----------------
  iao(1:n2+1) = 0
  do  i=1, n
!NEC$ ivdep
     do  k=ia(i), ia(i+1)-1
        j = ja(k)+1
        iao(j) = iao(j)+1
     end do
  end do
!---------- compute pointers from lengths ------------------------------
  iao(1) = ipos
  do  i=1,n2
     iao(i+1) = iao(i) + iao(i+1)
  end do
!--------------- now do the actual copying -----------------------------
  do  i=1,n
!NEC$ ivdep
     do  k=ia(i),ia(i+1)-1
        j = ja(k)
        next = iao(j)
        if (job == 1)  ao(next) = a(k)
        jao(next) = i
        iao(j) = next+1
     end do
  end do
!-------------------------- reshift iao and leave ----------------------
  do  i=n2,1,-1
     iao(i+1) = iao(i)
  end do
  iao(1) = ipos
!--------------- end of csrcsc2 ----------------------------------------
!-----------------------------------------------------------------------
end subroutine csrcsc2
!-----------------------------------------------------------------------
subroutine csrcsc2_r4 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
  integer,  intent(in)   :: n
  integer,  intent(in)   :: n2
  integer,  intent(in)   :: job
  integer,  intent(in)   :: ipos
  real(sp), intent(in)   :: a(:)
  integer,  intent(in)   :: ja(:)
  integer,  intent(in)   :: ia(n+1)
  real(sp), intent(out)  :: ao(:)
  integer,  intent(out)  :: jao(:)
  integer,  intent(out)  :: iao(n2+1)
  !-----------------------------------------------------------------------
  integer :: i, j, k, next

  iao(1:n2+1) = 0
  do  i=1, n
!NEC$ ivdep
     do  k=ia(i), ia(i+1)-1
        j = ja(k)+1
        iao(j) = iao(j)+1
     end do
  end do
  iao(1) = ipos
  do  i=1,n2
     iao(i+1) = iao(i) + iao(i+1)
  end do
  do  i=1,n
!NEC$ ivdep
     do  k=ia(i),ia(i+1)-1
        j = ja(k)
        next = iao(j)
        if (job == 1)  ao(next) = a(k)
        jao(next) = i
        iao(j) = next+1
     end do
  end do
  do  i=n2,1,-1
     iao(i+1) = iao(i)
  end do
  iao(1) = ipos
end subroutine csrcsc2_r4
!-----------------------------------------------------------------------

SUBROUTINE csrbnd (n,a,ja,ia,job,abd,nabd,lowd,ml,mu,ierr)
  integer,  intent(in)        :: n
  real(wp), intent(in)        :: a(*)
  integer,  intent(in)        :: ja(*)
  integer,  intent(in)        :: ia(n+1)
  integer,  intent(in)        :: job
  integer,  intent(in)        :: nabd
  real(wp), intent(out)       :: abd(nabd,n)
  integer,  intent(inout)     :: lowd
  integer,  intent(inout)     :: ml
  integer,  intent(inout)     :: mu
  integer,  intent(out)       :: ierr

!-----------------------------------------------------------------------
!   Compressed Sparse Row  to  Banded (Linpack ) format.
!-----------------------------------------------------------------------
! this subroutine converts a general sparse matrix stored in
! compressed sparse row format into the banded format. for the
! banded format,the Linpack conventions are assumed (see below).
!-----------------------------------------------------------------------
! on entry:
!----------
! n = integer,the actual row dimension of the matrix.

! a,
! ja,
! ia    = input matrix stored in compressed sparse row format.

! job = integer. if job=1 then the values of the lower bandwith ml
!         and the upper bandwidth mu are determined internally.
!         otherwise it is assumed that the values of ml and mu
!         are the correct bandwidths on input. See ml and mu below.

! nabd  = integer. first dimension of array abd.

! lowd  = integer. this should be set to the row number in abd where
!         the lowest diagonal (leftmost) of A is located.
!         lowd should be  ( 1  .le.  lowd  .le. nabd).
!         if it is not known in advance what lowd should be
!         enter lowd = 0 and the default value lowd = ml+mu+1
!         will be chosen. Alternative: call routine getbwd from unary
!         first to detrermione ml and mu then define lowd accordingly.
!         (Note: the banded solvers in linpack use lowd=2*ml+mu+1. )

! ml = integer. equal to the bandwidth of the strict lower part of A
! mu = integer. equal to the bandwidth of the strict upper part of A
!         thus the total bandwidth of A is ml+mu+1.
!         if ml+mu+1 is found to be larger than lowd then an error
!         flag is raised (unless lowd = 0). see ierr.

! note:   ml and mu are assumed to have  the correct bandwidth values
!         as defined above if job is set to zero on entry.

! on return:
!-----------

! abd   = real array of dimension abd(nabd,n).
!         on return contains the values of the matrix stored in
!         banded form. The j-th column of abd contains the elements
!         of the j-th column of  the original matrix comprised in the
!         band ( i in (j-ml,j+mu) ) with the lowest diagonal at
!         the bottom row (row lowd). See details below for this format.

! ml = integer. equal to the bandwidth of the strict lower part of A
! mu = integer. equal to the bandwidth of the strict upper part of A
!         if job=1 on entry then these two values are internally computed.

! lowd  = integer. row number in abd where the lowest diagonal
!         (leftmost) of A is located on return. In case lowd = 0
!         on return, then it is defined to ml+mu+1 on return and the
!         lowd will contain this value on return. `

! ierr  = integer. used for error messages. On return:
!         ierr .eq. 0  :means normal return
!         ierr .eq. -1 : means invalid value for lowd. (either .lt. 0
!         or larger than nabd).
!         ierr .eq. -2 : means that lowd is not large enough and as
!         result the matrix cannot be stored in array abd.
!         lowd should be at least ml+mu+1, where ml and mu are as
!         provided on output.

!----------------------------------------------------------------------*
! Additional details on banded format.  (this closely follows the      *
! format used in linpack. may be useful for converting a matrix into   *
! this storage format in order to use the linpack  banded solvers).    *
!----------------------------------------------------------------------*
!             ---  band storage format  for matrix abd ---             *
! uses ml+mu+1 rows of abd(nabd,*) to store the diagonals of           *
! a in rows of abd starting from the lowest (sub)-diagonal  which  is  *
! stored in row number lowd of abd. the minimum number of rows needed  *
! in abd is ml+mu+1, i.e., the minimum value for lowd is ml+mu+1. the  *
! j-th  column  of  abd contains the elements of the j-th column of a, *
! from bottom to top: the element a(j+ml,j) is stored in  position     *
! abd(lowd,j), then a(j+ml-1,j) in position abd(lowd-1,j) and so on.   *
! Generally, the element a(j+k,j) of original matrix a is stored in    *
! position abd(lowd+k-ml,j), for k=ml,ml-1,..,0,-1, -mu.               *
! The first dimension nabd of abd must be .ge. lowd                    *
!                                                                      *
!     example [from linpack ]:   if the original matrix is             *
!                                                                      *
!              11 12 13  0  0  0                                       *
!              21 22 23 24  0  0                                       *
!               0 32 33 34 35  0     original banded matrix            *
!               0  0 43 44 45 46                                       *
!               0  0  0 54 55 56                                       *
!               0  0  0  0 65 66                                       *
!                                                                      *
! then  n = 6, ml = 1, mu = 2. lowd should be .ge. 4 (=ml+mu+1)  and   *
! if lowd = 5 for example, abd  should be:                             *
!                                                                      *
! untouched --> x  x  x  x  x  x                                       *
!               *  * 13 24 35 46                                       *
!               * 12 23 34 45 56    resulting abd matrix in banded     *
!              11 22 33 44 55 66    format                             *
!  row lowd--> 21 32 43 54 65  *                                       *
!                                                                      *
! * = not used                                                         *


!----------------------------------------------------------------------*
integer :: i, j, k, m, ii, mdiag

! first determine ml and mu.
!-----------------------------------------------------------------------
ierr = 0
!-----------
IF (job == 1) CALL getbwd(n,a,ja,ia,ml,mu)
m = ml+mu+1
IF (lowd == 0) lowd = m
IF (m > lowd)  ierr = -2
IF (lowd > nabd .OR. lowd < 0) ierr = -1
IF (ierr < 0) RETURN
!------------
DO   i=1,m
  ii = lowd -i+1
  DO  j=1,n
    abd(ii,j) = 0.0D0
  END DO
END DO
!---------------------------------------------------------------------
mdiag = lowd-ml
DO  i=1,n
  DO  k=ia(i),ia(i+1)-1
    j = ja(k)
    abd(i-j+mdiag,j) = a(k)
  END DO
END DO
RETURN
!------------- end of csrbnd -------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrbnd
!-----------------------------------------------------------------------

SUBROUTINE bndcsr (n,abd,nabd,lowd,ml,mu,a,ja,ia,LEN,ierr)

INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: nabd
REAL(WP), INTENT(IN)                     :: abd(nabd,*)
INTEGER, INTENT(IN)                      :: lowd
INTEGER, INTENT(IN)                      :: ml
INTEGER, INTENT(IN)                      :: mu
REAL(WP), INTENT(OUT)                    :: a(*)
INTEGER, INTENT(OUT)                     :: ja(*)
INTEGER, INTENT(OUT)                     :: ia(n+1)
INTEGER, INTENT(IN)                      :: LEN
INTEGER, INTENT(OUT)                     :: ierr
REAL(WP)  t

!-----------------------------------------------------------------------
! Banded (Linpack ) format   to    Compressed Sparse Row  format.
!-----------------------------------------------------------------------
! on entry:
!----------
! n = integer,the actual row dimension of the matrix.

! nabd  = first dimension of array abd.

! abd   = real array containing the values of the matrix stored in
!         banded form. The j-th column of abd contains the elements
!         of the j-th column of  the original matrix,comprised in the
!         band ( i in (j-ml,j+mu) ) with the lowest diagonal located
!         in row lowd (see below).

! lowd  = integer. this should be set to the row number in abd where
!         the lowest diagonal (leftmost) of A is located.
!         lowd should be s.t.  ( 1  .le.  lowd  .le. nabd).
!         The subroutines dgbco, ... of linpack use lowd=2*ml+mu+1.

! ml = integer. equal to the bandwidth of the strict lower part of A
! mu = integer. equal to the bandwidth of the strict upper part of A
!         thus the total bandwidth of A is ml+mu+1.
!         if ml+mu+1 is found to be larger than nabd then an error
!         message is set. see ierr.

! len   = integer. length of arrays a and ja. bndcsr will stop if the
!         length of the arrays a and ja is insufficient to store the
!         matrix. see ierr.

! on return:
!-----------
! a,
! ja,
! ia    = input matrix stored in compressed sparse row format.

! lowd  = if on entry lowd was zero then lowd is reset to the default
!         value ml+mu+l.

! ierr  = integer. used for error message output.
!         ierr .eq. 0 :means normal return
!         ierr .eq. -1 : means invalid value for lowd.
!   ierr .gt. 0 : means that there was not enough storage in a and ja
!         for storing the ourput matrix. The process ran out of space
!         (as indicated by len) while trying to fill row number ierr.
!         This should give an idea of much more storage might be required.
!         Moreover, the first irow-1 rows are correctly filled.

! notes:  the values in abd found to be equal to zero
! -----   (actual test: if (abd(...) .eq. 0.0d0) are removed.
!         The resulting may not be identical to a csr matrix
!         originally transformed to a bnd format.

!-----------------------------------------------------------------------
integer :: i, j, ko, irow

ierr = 0
!-----------
IF (lowd > nabd .OR. lowd <= 0) THEN
  ierr = -1
  RETURN
END IF
!-----------
ko = 1
ia(1) = 1
DO  irow=1,n
!-----------------------------------------------------------------------
  i = lowd
  DO   j=irow-ml,irow+mu
    IF (j <= 0 ) GO TO 19
    IF (j > n) EXIT
    t = abd(i,j)
    IF (t == 0.0D0) GO TO 19
    IF (ko > LEN) THEN
      ierr = irow
      RETURN
    END IF
    a(ko) = t
    ja(ko) = j
    ko = ko+1
19  i = i-1
  END DO
!     end for row irow
  ia(irow+1) = ko
END DO
RETURN
!------------- end of bndcsr -------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE bndcsr
!-----------------------------------------------------------------------
subroutine csrjad (nrow, a, ja, ia, idiag, iperm, ao, jao, iao)
  integer,  intent(in)   :: nrow
  real(wp), intent(in)   :: a(:)
  integer,  intent(in)   :: ja(:)
  integer,  intent(in)   :: ia(nrow+1)
  integer,  intent(out)  :: idiag
  integer,  intent(out)  :: iperm(nrow)
  real(wp), intent(out)  :: ao(:)
  integer,  intent(out)  :: jao(:)
  integer,  intent(out)  :: iao(nrow+1)
  !-----------------------------------------------------------------------
  !    Compressed Sparse Row  to   JAgged Diagonal storage.
  !-----------------------------------------------------------------------
  ! this subroutine converts a matrix stored in the compressed sparse
  ! row format to the jagged diagonal format. The data structure
  ! for the JAD (Jagged Diagonal storage) is as follows. The rows of
  ! the matrix are (implicitly) permuted so that their lengths are in
  ! decreasing order. The real entries ao(*) and their column indices
  ! jao(*) are stored in succession. The number of such diagonals is idiag.
  ! the lengths of each of these diagonals is stored in iao(*).
  ! For more details see [E. Anderson and Y. Saad,
  ! ``Solving sparse triangular systems on parallel computers'' in
  ! Inter. J. of High Speed Computing, Vol 1, pp. 73-96 (1989).]
  ! or  [Y. Saad, ``Krylov Subspace Methods on Supercomputers''
  ! SIAM J. on  Stat. Scient. Comput., volume 10, pp. 1200-1232 (1989).]
  !-----------------------------------------------------------------------
  ! on entry:
  !----------
  ! nrow    = row dimension of the matrix A.
  !
  ! a,
  ! ia,
  ! ja      = input matrix in compressed sparse row format.
  !
  ! on return:
  !----------
  ! idiag = integer. The number of jagged diagonals in the matrix.
  !
  ! iperm = integer array of length nrow containing the permutation
  !         of the rows that leads to a decreasing order of the
  !         number of nonzero elements.
  !
  ! ao    = real array containing the values of the matrix A in
  !         jagged diagonal storage. The j-diagonals are stored
  !         in ao in sequence.
  !
  ! jao   = integer array containing the column indices of the
  !         entries in ao.
  !
  ! iao   = integer array containing pointers to the beginning
  !         of each j-diagonal in ao, jao. iao is also used as
  !         a work array and it should be of length n+1 at least.
  !-----------------------------------------------------------------------
  !     ---- define initial iperm and get lengths of each row
  !     ---- jao is used as work vector to store these lengths

  integer :: i, j, jj, k, k0, k1, ilo

  do  j=1, nrow
     iperm(j) = j
  end do
  jao(1:nrow) = ia(2:nrow+1) - ia(1:nrow)
  ilo   = minval (jao(1:nrow))
  idiag = maxval (jao(1:nrow))

  ! Call sorter to get permutation. use iao as work array.

  call dcsort (jao, nrow, iao, iperm, ilo, idiag)

  ! Define output data structure. first lengths of j-diagonals

  iao(1:nrow+1) = 0
  do  k=1, nrow
     do  i=1, jao(iperm(k))
        iao(i) = iao(i)+1
     end do
  end do

  ! Get the output matrix itself

  k1 = 1
  k0 = k1
  do  jj=1, idiag               ! Loop over 'diagonals'
     do  k=1, iao(jj)           ! Loop over (permuted) rows
        i = ia(iperm(k))+jj-1
        ao (k1) = a (i)
        jao(k1) = ja(i)
        k1 = k1+1
     end do
     iao(jj) = k0
     k0 = k1
  end do
  iao(idiag+1) = k1
!----------end-of-csrjad------------------------------------------------
!-----------------------------------------------------------------------
end subroutine csrjad
!-----------------------------------------------------------------------
subroutine csrjad_r4 (nrow, a, ja, ia, idiag, iperm, ao, jao, iao)
  integer,  intent(in)   :: nrow
  real(sp), intent(in)   :: a(:)
  integer,  intent(in)   :: ja(:)
  integer,  intent(in)   :: ia(nrow+1)
  integer,  intent(out)  :: idiag
  integer,  intent(out)  :: iperm(nrow)
  real(sp), intent(out)  :: ao(:)
  integer,  intent(out)  :: jao(:)
  integer,  intent(out)  :: iao(nrow+1)

  integer :: i, j, jj, k, k0, k1, ilo

  do  j=1, nrow
     iperm(j) = j
  end do
  jao(1:nrow) = ia(2:nrow+1) - ia(1:nrow)
  ilo   = minval (jao(1:nrow))
  idiag = maxval (jao(1:nrow))

  call dcsort (jao, nrow, iao, iperm, ilo, idiag)

  iao(1:nrow+1) = 0
  do  k=1, nrow
     do  i=1, jao(iperm(k))
        iao(i) = iao(i)+1
     end do
  end do

  k1 = 1
  k0 = k1
  do  jj=1, idiag               ! Loop over 'diagonals'
     do  k=1, iao(jj)           ! Loop over (permuted) rows
        i = ia(iperm(k))+jj-1
        ao (k1) = a (i)
        jao(k1) = ja(i)
        k1 = k1+1
     end do
     iao(jj) = k0
     k0 = k1
  end do
  iao(idiag+1) = k1
end subroutine csrjad_r4
!-----------------------------------------------------------------------
subroutine jadcsr (nrow, idiag, a, ja, ia, iperm, ao, jao, iao)
  integer,  intent(in)   :: nrow
  integer,  intent(in)   :: idiag
  real(wp), intent(in)   :: a(:)
  integer,  intent(in)   :: ja(:)
  integer,  intent(in)   :: ia(idiag+1)
  integer,  intent(in)   :: iperm(nrow)
  real(wp), intent(out)  :: ao(:)
  integer,  intent(out)  :: jao(:)
  integer,  intent(out)  :: iao(nrow+1)
  !-----------------------------------------------------------------------
  !     Jagged Diagonal Storage   to     Compressed Sparse Row
  !-----------------------------------------------------------------------
  ! this subroutine converts a matrix stored in the jagged diagonal format
  ! to the compressed sparse row format.
  !-----------------------------------------------------------------------
  ! on entry:
  !----------
  ! nrow    = integer. the row dimension of the matrix A.
  !
  ! idiag   = integer. The  number of jagged diagonals in the data
  !           structure a, ja, ia.
  !
  ! a,
  ! ja,
  ! ia      = input matrix in jagged diagonal format.
  !
  ! iperm   = permutation of the rows used to obtain the JAD ordering.
  !
  ! on return:
  !----------
  ! ao, jao,
  ! iao     = matrix in CSR format.
  !-----------------------------------------------------------------------
  ! determine first the pointers for output matrix. Go through the
  ! structure once:

  integer :: i, j, jj, k, k1, kpos

  ! Compute the lengths of each row of output matrix -

  jao(1:nrow) = 0
!NEC$ novector
  do  i=1, idiag
!NEC$ ivdep
     do  k=1,ia(i+1)-ia(i)
        jao(iperm(k)) = jao(iperm(k))+1
     end do
  end do

  ! Remember to permute

  kpos = 1
  iao(1) = 1
  do  i=1, nrow
     kpos = kpos+jao(i)
     iao(i+1) = kpos
  end do

  ! Copy elements one at a time.

  do  jj = 1, idiag
     k1 = ia(jj)-1
!NEC$ ivdep
     do  k=1, ia(jj+1)-k1-1
        kpos = iao(iperm(k))
        ao (kpos) = a (k1+k)
        jao(kpos) = ja(k1+k)
        iao(iperm(k)) = kpos+1
     end do
  end do

  ! Rewind pointers

  do  j=nrow,1,-1
     iao(j+1) = iao(j)
  end do
  iao(1) = 1
!----------end-of-jadcsr------------------------------------------------
!-----------------------------------------------------------------------
end subroutine jadcsr
!-----------------------------------------------------------------------
subroutine jadcsr_r4 (nrow, idiag, a, ja, ia, iperm, ao, jao, iao)
  integer,  intent(in)   :: nrow
  integer,  intent(in)   :: idiag
  real(sp), intent(in)   :: a(:)
  integer,  intent(in)   :: ja(:)
  integer,  intent(in)   :: ia(idiag+1)
  integer,  intent(in)   :: iperm(nrow)
  real(sp), intent(out)  :: ao(:)
  integer,  intent(out)  :: jao(:)
  integer,  intent(out)  :: iao(nrow+1)

  integer :: i, j, jj, k, k1, kpos

  jao(1:nrow) = 0
!NEC$ novector
  do  i=1, idiag
!NEC$ ivdep
     do  k=1,ia(i+1)-ia(i)
        jao(iperm(k)) = jao(iperm(k))+1
     end do
  end do

  kpos = 1
  iao(1) = 1
  do  i=1, nrow
     kpos = kpos+jao(i)
     iao(i+1) = kpos
  end do

  do  jj = 1, idiag
     k1 = ia(jj)-1
!NEC$ ivdep
     do  k=1, ia(jj+1)-k1-1
        kpos = iao(iperm(k))
        ao (kpos) = a (k1+k)
        jao(kpos) = ja(k1+k)
        iao(iperm(k)) = kpos+1
     end do
  end do

  do  j=nrow,1,-1
     iao(j+1) = iao(j)
  end do
  iao(1) = 1
end subroutine jadcsr_r4
!-----------------------------------------------------------------------
subroutine dcsort (ival, n, icnt, index, ilo, ihi)
!-----------------------------------------------------------------------
!     Specifications for arguments:
!     ----------------------------
  integer, intent(in)  :: n
  integer, intent(in)  :: ival(n)
  integer, intent(in)  :: ilo
  integer, intent(in)  :: ihi
  integer, intent(out) :: icnt(ilo:ihi)
  integer, intent(out) :: index(n)
!-----------------------------------------------------------------------
!    This routine computes a permutation which, when applied to the
!    input vector ival, sorts the integers in ival in descending
!    order.  The permutation is represented by the vector index.  The
!    permuted ival can be interpreted as follows:
!      ival(index(i-1)) .ge. ival(index(i)) .ge. ival(index(i+1))
!
!    A specialized sort, the distribution counting sort, is used
!    which takes advantage of the knowledge that
!        1)  The values are in the (small) range [ ilo, ihi ]
!        2)  Values are likely to be repeated often
!
!    contributed to SPARSKIT by Mike Heroux. (Cray Research)
!    ---------------------------------------
!-----------------------------------------------------------------------
! Usage:
!------
!     call dcsort( ival, n, icnt, index, ilo, ihi )
!
! Arguments:
!-----------
!    ival  integer array (input)
!          On entry, ia is an n dimensional array that contains
!          the values to be sorted.  ival is unchanged on exit.
!
!    n     integer (input)
!          On entry, n is the number of elements in ival and index.
!
!    icnt  integer (work)
!          On entry, is an integer work vector of length
!          (ihi - ilo + 1).
!
!    index integer array (output)
!          On exit, index is an n-length integer vector containing
!          the permutation which sorts the vector ival.
!
!    ilo   integer (input)
!          On entry, ilo is .le. to the minimum value in ival.
!
!    ihi   integer (input)
!          On entry, ihi is .ge. to the maximum value in ival.
!
! Remarks:
!---------
!         The permutation is NOT applied to the vector ival.
!----------------------------------------------------------------
! Local variables:
!    Other integer values are temporary indices.
!
! Author:
!--------
!    Michael Heroux
!    Sandra Carney
!       Mathematical Software Research Group
!       Cray Research, Inc.
!
! References:
!    Knuth, Donald E., "The Art of Computer Programming, Volume 3:
!    Sorting and Searching," Addison-Wesley, Reading, Massachusetts,
!    1973, pp. 78-79.
!
! Revision history:
!    05/09/90: Original implementation.  A variation of the
!              Distribution Counting Sort recommended by
!              Sandra Carney. (Mike Heroux)
!-----------------------------------------------------------------
!     ----------------------------------
!     Specifications for local variables
!     ----------------------------------
  integer :: i, j, ivalj

!     --------------------------
!     First executable statement
!     --------------------------
  icnt(ilo:ihi) = 0
  do  i = 1, n
     icnt(ival(i)) = icnt(ival(i)) + 1
  end do

  do  i = ihi-1,ilo,-1
     icnt(i) = icnt(i) + icnt(i+1)
  end do

  do  j = n, 1, -1
     ivalj = ival(j)
     index(icnt(ivalj)) = j
     icnt(ivalj) = icnt(ivalj) - 1
  end do

end subroutine dcsort
!-------end-of-dcsort---------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE csorted(n, ia, ja, sorted)
  integer, intent(in)     :: n
  integer, intent(in)     :: ia(n+1)
  integer, intent(in)     :: ja(:)
  logical, intent(out)    :: sorted
!-----------------------------------------------------------------------
!     Checks if matrix in CSR format is sorted by columns.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     n       = number of rows in matrix
!     ia, ja  = sparsity structure of matrix in CSR format

!     On return:
!---------------
!     sorted  = indicates if matrix is sorted by columns

!-----------------------------------------------------------------------
!-----local variables
INTEGER :: i,j
!---------------------------------
DO i = 1, n
  DO j = ia(i)+1, ia(i+1)-1
    IF (ja(j-1) >= ja(j)) THEN
      sorted = .false.
      RETURN
    END IF
  END DO
END DO
sorted = .true.
RETURN
END SUBROUTINE csorted
!-----------------------------------------------------------------------
!------------------------end-of-csorted---------------------------------

SUBROUTINE dvperm (n, x, perm)
  INTEGER,  INTENT(IN)                  :: n
  REAL(WP), INTENT(INOUT)               :: x(n)
  INTEGER,  INTENT(INOUT)               :: perm(n)
!-----------------------------------------------------------------------
! this subroutine performs an in-place permutation of a real vector x
! according to the permutation array perm(*), i.e., on return,
! the vector x satisfies,

! x(perm(j)) :== x(j), j=1,2,.., n

!-----------------------------------------------------------------------
! on entry:
!---------
! n  = length of vector x.
! perm  = integer array of length n containing the permutation  array.
! x = input vector

! on return:
!----------
! x = vector x permuted according to x(perm(*)) :=  x(*)

!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables
  integer  :: ii, j, k, init, next
  REAL(WP) :: tmp, tmp1

  init      = 1
  tmp       = x(init)
  ii        = perm(init)
  perm(init)= -perm(init)
  k         = 0

! loop

  outer: do
     k = k+1

     ! save the chased element --

     tmp1   = x(ii)
     x(ii)  = tmp
     next   = perm(ii)
     IF (next < 0 ) GO TO 65

     ! test for end

     IF (k > n) exit outer
     tmp       = tmp1
     perm(ii)  = - perm(ii)
     ii        = next

     ! end loop

     cycle outer

     ! reinitialize cycle --

65   init = init+1
     IF (init > n) exit outer
     IF (perm(init) < 0) GO TO 65
     tmp = x(init)
     ii = perm(init)
     perm(init)=-perm(init)
  end do outer

  DO  j=1, n
     perm(j) = -perm(j)
  END DO
!-------------------end-of-dvperm---------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE dvperm
!-----------------------------------------------------------------------
SUBROUTINE getbwd(n,a,ja,ia,ml,mu)
!-----------------------------------------------------------------------
! gets the bandwidth of lower part and upper part of A.
! does not assume that A is sorted.
!-----------------------------------------------------------------------
! on entry:
!----------
! n = integer = the row dimension of the matrix
! a, ja,
!    ia = matrix in compressed sparse row format.

! on return:
!-----------
! ml = integer. The bandwidth of the strict lower part of A
! mu = integer. The bandwidth of the strict upper part of A

! Notes:
! ===== ml and mu are allowed to be negative or return. This may be
!       useful since it will tell us whether a band is confined
!       in the strict  upper/lower triangular part.
!       indeed the definitions of ml and mu are

!       ml = max ( (i-j)  s.t. a(i,j) .ne. 0  )
!       mu = max ( (j-i)  s.t. a(i,j) .ne. 0  )
!----------------------------------------------------------------------c
! Y. Saad, Sep. 21 1989                                                c
!----------------------------------------------------------------------c
  INTEGER,  INTENT(IN)     :: n
  REAL(WP), INTENT(IN)     :: a(*)
  INTEGER,  INTENT(IN)     :: ja(*)
  INTEGER,  INTENT(IN)     :: ia(n+1)
  INTEGER,  INTENT(OUT)    :: ml
  INTEGER,  INTENT(OUT)    :: mu

  INTEGER :: ldist,i,k

  ml = - n
  mu = - n
  DO  i=1,n
     DO  k=ia(i),ia(i+1)-1
        ldist = i-ja(k)
        ml = MAX(ml,ldist)
        mu = MAX(mu,-ldist)
     END DO
  END DO
!---------------end-of-getbwd ------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE getbwd
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE dump (i1,i2,values,a,ja,ia,iout)
  integer,  intent(in)                      :: i1
  integer,  intent(in)                      :: i2
  logical,  intent(in)                      :: values
  real(wp), intent(in)                      :: a(*)
  integer,  intent(in)                      :: ja(*)
  integer,  intent(in)                      :: ia(*)
  integer,  intent(in)                      :: iout
!-----------------------------------------------------------------------
! outputs rows i1 through i2 of a sparse matrix stored in CSR format
! (or columns i1 through i2 of a matrix stored in CSC format) in a file,
! one (column) row at a time in a nice readable format.
! This is a simple routine which is useful for debugging.
!-----------------------------------------------------------------------
! on entry:
!---------
! i1    = first row (column) to print out
! i2    = last row (column) to print out
! values= logical. indicates whether or not to print real values.
!         if value = .false. only the pattern will be output.
! a,
! ja,
! ia    =  matrix in CSR format (or CSC format)
! iout  = logical unit number for output.
!----------
! the output file iout will have written in it the rows or columns
! of the matrix in one of two possible formats (depending on the max
! number of elements per row. The values are output with only
! two digits of accuracy (D9.2). )
!-----------------------------------------------------------------------
!     local variables

  INTEGER :: maxr, i, k, k1, k2

! select mode horizontal or vertical

  maxr = 0
  DO  i=i1, i2
     maxr = MAX0(maxr,ia(i+1)-ia(i))
  END DO

  IF (maxr <= 8) THEN

! able to do one row acros line

     DO  i=i1, i2
        WRITE(iout,100) i
        k1=ia(i)
        k2 = ia(i+1)-1
        WRITE (iout,101) (ja(k),k=k1,k2)
        IF (values) WRITE (iout,102) (a(k),k=k1,k2)
     END DO
  ELSE

! unable to one row acros line. do three items at a time
! across a line
     DO  i=i1, i2
        IF (values) THEN
           WRITE(iout,200) i
        ELSE
           WRITE(iout,203) i
        END IF
        k1=ia(i)
        k2 = ia(i+1)-1
        IF (values) THEN
           WRITE (iout,201) (ja(k),a(k),k=k1,k2)
        ELSE
           WRITE (iout,202) (ja(k),k=k1,k2)
        END IF
     END DO
  END IF

! formats :

100 FORMAT (" ",34('-'),' row',i6,1X,34('-') )
101 FORMAT(' col:',8(i5,'     :'))
102 FORMAT(' val:',8(es9.2,' :') )
200 FORMAT (" ",30('-'),' row',i3,1X,30('-'),/ 3('  columns :    values  * '))
!-------------xiiiiiihhhhhhddddddddd-*-
201 FORMAT(3(" ",i6,'   :  ',es9.2,' * '))
202 FORMAT(6(" ",i5,'  *   '))
203 FORMAT (" ",30('-'),' row',i3,1X,30('-'),/ 3('  column  :  column   *') )
!----end-of-dump--------------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE dump
!-----------------------------------------------------------------------

  !============================================================================
  subroutine jad_stat (ia)
    integer, intent(in)  :: ia(:)       ! Indices of jagged diagonal heads
    !------------------------------------------------------------------
    ! Derive statistics from metadata of jagged diagonal storage format
    !------------------------------------------------------------------
    integer, allocatable :: jad_len(:)  ! Lengths of jagged diagonals
    integer, allocatable :: mult(:)     ! Multiplicities of lengths
    integer, allocatable :: hist(:,:)   ! Histogram multiplicity vs. length
    integer              :: jdiag       ! Number of jagged-diagonals
    integer :: i, k, n, n1, n2
    integer :: lmin, lmax, mmin, mmax, mlarge, mdiff
    real    :: l_av, l_rms, m_av, m_rms

    jdiag = size (ia) - 1
    allocate (jad_len(jdiag))
    jad_len(1:jdiag) = ia(2:jdiag+1) - ia(1:jdiag)
    lmin  = minval (jad_len)
    lmax  = maxval (jad_len)
    l_av  = real (ia(jdiag+1) - ia(1)) / jdiag
    l_rms = sqrt (sum (real (jad_len)**2) / jdiag)
    allocate (mult(0:lmax))
    mult = 0
    do i = 1, jdiag
       k = jad_len(i)
       mult(k) = mult(k) + 1
    end do
    mmax = maxval (mult)
    mmin = minval (mult, mult > 0)
    mlarge = 0
    do i = lmax, 1, -1
       if (mult(i) > 0) then
          mlarge = i
          exit
       end if
    end do
    mdiff = count (mult > 0)
    m_av  = 0
    m_rms = 0
    do i = 1, mlarge
       m_av  = m_av  + mult(i)
       m_rms = m_rms + mult(i)**2
    end do
    m_av  =       m_av  / mdiff
    m_rms = sqrt (m_rms / mdiff)

    n2 = int (log (2*lmax - 0.99)/log (2.0))    ! 2^(n2-1) < lmax <= 2^n2
    allocate (hist(0:n2,5))
    hist = 0
    do i = 1, lmax
       if (mult(i) > 0) then
          k = int (log (2*i - 0.99)/log (2.0))
          n = min (mult(i), 5)
          hist(k,n) = hist(k,n) + 1
       end if
    end do

    print *
    print '(A)',       "  Sparse Matrix metadata statistics (JAD):"
    print '(A)',       "============================================"
    print '(A,i10)',   "Number of Jagged Diagonals       :", jdiag
    print '(A,i10)',   "Length of longest  diagonal      :", lmax
    print '(A,i10)',   "Length of shortest diagonal      :", lmin
    print '(A,f10.3)', "Mean length of diagonals         :", l_av
    print '(A,f10.3)', "RMS  length of diagonals         :", l_rms
    print '(A,i10)',   "Maximum multiplicity of lengths  :", mmax
    print '(A,i10)',   "Minimum multiplicity of lengths  :", mmin
    print '(A,i10)',   "Multiplicity of longest diagonal :", mult(mlarge)
    print '(A,i10)',   "Diagonals with multiplicity unity:", sum (hist(:,1))
    print '(A,i10)',   "Count of different multiplicities:", mdiff
    print '(A,f10.3)', "Mean multiplicity of lengths     :", m_av
    print '(A,f10.3)', "RMS  multiplicity of lengths     :", m_rms
    print *
    print '(A)', "Diagonal length vs. multiplicity"
    print *
    print '(4x,A)',"Length range:    Multiplicity:   1     2     3     4   >=5"
    n1 = 0
    do n = 0, n2
       if (sum (hist(n,:)) == 0) cycle
       n1 = n
       exit
    end do
    do n = n1, n2
       print '(i7,A,i7,15x,5i6)', 2**(n-1)+1, " ..", min (2**n,lmax), hist(n,:)
    end do
    print *
  end subroutine jad_stat
  !============================================================================

subroutine test_convert ()
  !-----------------------------------------------------------------------
  ! This routine is derived from chkfmt1.f and rmatvec.f
  !-----------------------------------------------------------------------
  integer :: ndns, nxmax, nmx, nnzmax
  parameter (nxmax = 10, nmx = nxmax*nxmax, nnzmax=10*nmx)
  integer :: ia(nmx+1),ja(nnzmax),ia1(nnzmax),ja1(nnzmax),  &
             iwk(nmx*2+1),iperm(nmx),ia2(nnzmax),ja2(nnzmax)
  real(wp) :: a(nnzmax), a1(nnzmax), a2(nnzmax), dns(20,20)
  real(wp) :: tmp
  integer  :: k,nx,ny,n,iout,nnz,i,ierr, ndiag,k1,k2
  real(wp), allocatable :: adense(:,:), x(:), y(:), z(:), yy(:)
  !-----------------------------------------------------------------------
  data ndns/20/
  !-----------------------------------------------------------------------
  print *, "Test of selected conversion routines of SPARSKIT"
  print *
  !---- Generate test matrix in dense representation
  nx = 5
  ny = 5
  allocate (adense (nx,ny))
  adense = 0
  do i = 1, nx
     if (i <= ny) adense(i,i) = i + 0.1*i
  end do
  adense(1,2) = 1.2
  adense(3,2) = 3.2
  adense(4,1) = 4.1
  adense(1,5) = 1.5
  adense(5,5) = 0       ! Last row all zeros
  ! These are useful for testing the grouped version:
  adense(4,3) = 4.3
  adense(2,3) = 2.3
  print *, "Dense test matrix A:"
  print *, "===================="
  do i = 1, nx
     print '(99f4.1)', adense(i,:)
  end do
  ! Test vectors for matrix-vector multiplication
  allocate (x(ny), y(nx), z(nx), yy(nx))
  call random_number (x)
  z = matmul (adense, x)
  print *
  print *, "Test vector x:"
  print *, real (x)
  print *
  print *, "   z = A*x:", real (z)
  print *
  n = nx
  ia = -1
  ja = -1
  print *, "Conversion Dense -> CSR"
  call dnscsr (nx,ny,nnzmax,adense,nx,a,ja,ia,ierr)
  print *
  print '(A,99i4)',   "ia: ", ia(1:n+1)
  print '(A,99i4)',   "ja: ", ja(1:ia(n+1)-1)
  print '(A,99f4.1)', " a: ",  a(1:ia(n+1)-1)
  print *
  !call dump(1,nx,.false.,a,ja,ia,6)
  call dump(1,nx,.true. ,a,ja,ia,6)

  print *
  print *, "Conversion CSR -> Dense:"
  call csrdns (n,n,a,ja,ia,dns,ndns,ierr)
  if (any (dns(1:n,1:n) /= adense)) then
     print *, "FAIL!"
     stop
  end if
  print *, "OK."
  print *
  print *, "Testing amux:"
  print *
  call amux (n, x, y, a, ja, ia)
  print *, "   y = A*x:", REAL (y)
  tmp = sum (ABS (y-z)**2)
  print *, "||y-z||^2 =", tmp
  print *
  if (tmp > 1e-12) then
     print *, "===== FAIL! ====="
     stop
  end if

  print *
  print *, "Conversion CSR -> CSC:"
  call csrcsc (n,1,1,a,ja,ia,a1,ja1,ia1)
  print *, "Conversion CSC -> CSR:"
  call csrcsc (n,1,1,a1,ja1,ia1,a2,ja2,ia2)
  nnz = ia(n+1)-1
  if ( any (ia2(1:n+1) /= ia(1:n+1)) .or. &
       any (ja2(1:nnz) /= ja(1:nnz)) .or. &
       any ( a2(1:nnz) /=  a(1:nnz))) then
     call dump (1,n,.true.,a2,ja2,ia2,6)
     print *, "FAIL!"
     stop
  end if
  print *, "OK."
  print *
  print *, "Checking atmux"
  call atmux  (n, x, y, a1, ja1, ia1)
  print *
  print *, "y = (A~)^T*x:", real (y)
  tmp = sum (abs (y-z)**2)
  print *, "||y-z||^2 =", tmp
  print *
  if (tmp > 1e-12) then
     print *, "===== FAIL! ====="
     stop
  end if

  print *
  print *, "Conversion CSR -> JAD format:"
  print *
  call csrjad (n, a, ja, ia, ndiag, iwk, a1, ja1, ia1)

  iperm(1:n) = iwk(1:n)
  print *, "ndiag =", ndiag
  print '(a,99i4)',   "ia: ", ia1(1:n+1)
  print '(a,99i4)',   "ja: ", ja1(1:ia(n+1)-1)
  print '(a,99f4.1)', " a: ",  a1(1:ia(n+1)-1)
  print *
  iout = 6
  write (iout,*) '* PERMUTATION ARRAY *'
  write (iout,'(17i4)') (iperm(k),k=1,n)
  print *
  do i=1,ndiag
     write (iout,*) ' J-diagonal number: ', i
     k1 = ia1(i)
     k2 = ia1(i+1)-1
     write (iout,'("JCOEF = ",99i4)')  (ja1(k),k=k1,k2)
     write (iout,'("COEF  = ",99f4.1)') (a1(k),k=k1,k2)
  end do

  print *
  call jad_stat (ia1(1:ndiag+1))
  print *
  call amuxj (n, x, y, ndiag, a1, ja1, ia1)
  yy(iperm(1:n)) = y(1:n)
  call dvperm (n, y, iperm)
  print *, "   y = A*x:", real (y)
  tmp = sum (abs (y-z)**2)
  print *, "||y-z||^2 =", tmp
  print *
  if (tmp > 1e-12) then
     print *, "===== FAIL! ====="
     stop
  end if
  if (maxval (abs (y-yy)) > 0) then
     print *, "===== Permutation FAILED! ====="
     stop
  end if
  print *, "Permutation (dvperm) OK."
  print *
  print *, "Conversion JAD format -> CSR:"
  call jadcsr (n, ndiag, a1, ja1, ia1, iperm, a2, ja2, ia2)
  nnz = ia(n+1)-1
  if ( any (ia2(1:n+1) /= ia(1:n+1)) .or. &
       any (ja2(1:nnz) /= ja(1:nnz)) .or. &
       any ( a2(1:nnz) /=  a(1:nnz))) then
     write (*,*) '-------------------------------------------'
     write (*,*) '+++ matrix after conversion from jadcsr +++'
     write (*,*) '-------------------------------------------'
     call dump (1,n,.true.,a2,ja2,ia2,6)
     print *, "FAIL!"
     stop
  end if
  print *, "OK."
  print *
end subroutine test_convert

  !-----------------------------------------------------------------------

end module mo_sparskit

!program chkfmt
!  use mo_sparskit
!  implicit none
!  call test_convert ()
!end program chkfmt
