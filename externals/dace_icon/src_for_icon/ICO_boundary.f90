!
!+ GNSS Radio occultation operator: process boundary points of ICO grid
!
MODULE ICO_boundary
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Processing boundary points of ICO grid
!   for producing enlarged arrays.
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
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Reference:
!   Michael E. Gorbunov and Luis Kornblueh
!   Principles of variational assimilation of GNSS radio occultation data.
!   Max-Planck-Institut fuer Meteorologie, Hamburg, Report No. 350 (2003)
!
! Author:
! Michael E. Gorbunov  2004  original code
! Changes:
! Andreas Rhodin             adapted to DWD 3D-VAR
!==============================================================================


!
! Module ICO_boundary
!
! Processing boundary points of ICO grid
! for producing enlarged arrays.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 29 Nov 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, Pi, dtr, rtd
!----------------------------------------------------------
Implicit None
Private
Public :: Set_Boundaries
!----------------------------------------------------------
! Public Parameters:
!
Integer, Parameter :: &
   nd = 10              ! Number of diamonds
!
! --- Numbers of diamond neighbours in different directions
!
Integer, Parameter :: & ! Poleward westward neighbours
   mpw(nd) = (/ 5, 1, 2, 3, 4,10, 6, 7, 8, 9 /)
Integer, Parameter :: & ! Poleward eastward neighbours
   mpe(nd) = (/ 2, 3, 4, 5, 1, 7, 8, 9,10, 6 /)
Integer, Parameter :: & ! Antipoleward westward neighbours
   maw(nd) = (/10, 6, 7, 8, 9, 1, 2, 3, 4, 5 /)
Integer, Parameter :: & ! Antipoleward eastward neighbours
   mae(nd) = (/ 6, 7, 8, 9,10, 2, 3, 4, 5, 1 /)
!
! Public Scalars:
!
Integer :: &
    np_1,          & ! number of points to process - 1 boundary line
    np_2,          & ! number of points to process - 2 boundary lines
    np_max           ! max. number of points to process
!
! Public Arrays:
!
! --- Permutation indexes for enlarging arrays at boundaries
! idx[1/2/3,m] - j1/j2/jd indices of m-th element in buffer
!
Integer, Allocatable ::  &
    idx_recv(:,:), & ! Indices of received points
    idx_send(:,:)    ! Indices of points to send
!----------------------------------------------------------
!
Contains


!==========================================================
Subroutine Setup_Boundary &
  (kgg1s, kgg1e,  & ! <-- = kg1s, kg1e (for compatibility)
   kgg2s, kgg2e,  & ! <-- = kg2s, kg2e (for compatibility)
   knd,           & ! <-- Number of diamonds (must be = nd)
   kg1s, kg1e,    & ! <-- Core interval of 1st index 0..ni
   kg2s, kg2e,    & ! <-- Core interval of 2nd index 1..ni+1
   kni)             ! <-- Triangles / diamond side
!
! Initialization of index arrays for enlarging
! grid boundaries by 1 or 2 rows.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 29 Nov 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   kgg1s, kgg1e     ! = kg1s, kg1e (for compatibility)
!
Integer, Intent(In) :: &
   kgg2s, kgg2e     ! = kg2s, kg2e (for compatibility)
!
Integer, Intent(In) :: &
   knd              ! Number of diamonds (must be = nd)
!
Integer, Intent(In) :: &
   kg1s, kg1e       ! Core interval of 1st index 0..ni
!
Integer, Intent(In) :: &
   kg2s, kg2e       ! Core interval of 2nd index 1..ni+1
!
Integer, Intent(In) :: &
   kni              ! Triangles / diamond side
!----------------------------------------------------------
! Local Scalars:
!
Integer :: jd, j1, j2, j  ! Gridpoint indices
Integer :: jr             ! Number of exchange lines
Integer :: mbsize         ! Exchange buffer size
Integer :: mb             ! Exchange buffer counter
Integer :: jb             ! 1-diamond exchange buffer index
Integer :: mi1sc, mi1ec   ! First index inner limits
Integer :: mi2sc, mi2ec   ! Second index inner limits
!
! Local Arrays:
!
! --- Core + extended indices:
! marr(1,.,.,.) processor number (=0, for compatibility)
! marr(2/3/4,.,.,.) j1/j2/jd indices
!
Integer :: marr(4,kgg1s-2:kgg1e+2,kgg2s-2:kgg2e+2,1:knd)
!
! --- Boundary indices for 1 diamond
! idx_bound(1/2,.) == j1/j2-index of boundary point
! idx_bound(3,.)   == 1/2 if point needed in 1/2 line exchange
! idx_bound(3,.)   == 3 used in 2 line exchange for compatibility
!
Integer :: idx_bound(3,(kg1e-kg1s+5)*(kg2e-kg2s+5) )
!----------------------------------------------------------
! Global variables used:
!
!   idx_recv    ! Indices of received points
!   idx_send    ! Indices of points to send
!----------------------------------------------------------


!----------------------------------------------------------
! 1. SETTING MARR ARRAY
!----------------------------------------------------------


!--- 1.0 Initialization

marr(:,:,:,:) = -1


!--- 1.1. Diamond cores

Do jd = 1, knd
   Do j2 = 1, kni
      Do j1 = 1, kni
         marr(1,j1,j2,jd) = 0
         marr(2,j1,j2,jd) = j1
         marr(3,j1,j2,jd) = j2
         marr(4,j1,j2,jd) = jd
      End Do
   End Do
End Do

!--- 1.2. Pole

Do jd =  1, knd
   marr(1,0,1,jd) = 0 ! Pole is always owned by proc 0
   marr(2,0,1,jd) = 0
   marr(3,0,1,jd) = 1
   If(jd <= 5) then
      marr(4,0,1,jd) = 1 ! Use always diamond 1 for north pole
   Else
      marr(4,0,1,jd) = 6 ! Use always diamond 6 for south pole
   End If
End Do

!--- 1.3. Boundaries

Do jd = 1, knd

   !--- 1.3.1. Poleward west

   Do j= 1, kni
      marr(:,j  , 0,jd) = marr(:,1,j,mpw(jd))
      marr(:,j+1,-1,jd) = marr(:,2,j,mpw(jd))
   End Do

   !--- 1.3.2. Poleward east

   Do j= 1, kni
      marr(:, 0,j+1,jd) = marr(:,j,1,mpe(jd))
      marr(:,-1,j+2,jd) = marr(:,j,2,mpe(jd))
      marr(:,-2,j+3,jd) = marr(:,j,3,mpe(jd))
   End Do

   !--- 1.3.3. Anti-poleward west

   Do j= 1, kni
      marr(:,kni+1,kni-j+1,jd) = marr(:,j,kni  ,maw(jd))
      marr(:,kni+2,kni-j+1,jd) = marr(:,j,kni-1,maw(jd))
   End Do

   !--- 1.3.4. Anti-poleward east

   Do j= 1, kni
      marr(:,kni-j+1,kni+1,jd) = marr(:,kni  ,j,mae(jd))
      marr(:,kni-j+1,kni+2,jd) = marr(:,kni-1,j,mae(jd))
      marr(:,kni-j+1,kni+3,jd) = marr(:,kni-2,j,mae(jd))
   End Do

End Do


!--- 1.4. Special points

Do jd = 1, knd

   !--- 1.4.1. Pole

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

   !--- 1.4.2. West Corner

   marr(:,kni+1, 0,jd) = marr(:,kni,0,jd) ! Mirror

   marr(:,kni+2,-1,jd) = marr(:,kni+1,-1,jd) ! Undefined
   marr(:,kni+2, 0,jd) = marr(:,kni+1,-1,jd) ! Mirror

   !--- 1.4.3. Antipole

   marr(:,kni+1,kni+1,jd) = marr(:,kni  ,kni+2,jd) ! Mirror
   marr(:,kni+1,kni+2,jd) = marr(:,kni  ,kni+2,jd) ! Undefined

   marr(:,kni+2,kni+1,jd) = marr(:,kni  ,kni+3,jd) ! Mirror
   marr(:,kni+2,kni+2,jd) = marr(:,kni  ,kni+3,jd) ! Undefined
   marr(:,kni+2,kni+3,jd) = marr(:,kni  ,kni+3,jd) ! Undefined
   marr(:,kni+1,kni+3,jd) = marr(:,kni  ,kni+3,jd) ! Undefined

   !--- 1.4.4. East Corner

   marr(:, 0,kni+2,jd) = marr(:,-1,kni+2,jd) ! Copy
   marr(:,-1,kni+2,jd) = marr(:,-1,kni+1,jd) ! Mirror

   marr(:, 0,kni+3,jd) = marr(:,-2,kni+3,jd) ! Copy
   marr(:,-2,kni+3,jd) = marr(:,-2,kni+2,jd) ! Undefined
   marr(:,-1,kni+3,jd) = marr(:,-2,kni+2,jd) ! Mirror

End Do


!----------------------------------------------------------
! 2. SETTING BOUNDARY INDEX ARRAY
!----------------------------------------------------------


!--- 2.1. Set computational boundaries

mi1sc = MAX(kg1s,1)
mi1ec = kg1e
mi2sc = kg2s
mi2ec = MIN(kg2e,kni)


!--- 2.2. Processing boundary points

mb = 0

Do j1 = kg1s-2, kg1e+2
   Do j2 = kg2s-2, kg2e+2

      If (j1 < mi1sc .OR. j1 > mi1ec .OR. &
          j2 < mi2sc .OR. j2 > mi2ec) then

         !--- 2.2.1. Adding boundary point indices

         mb = mb+1
         idx_bound(1,mb) = j1
         idx_bound(2,mb) = j2

         !--- 2.2.2. Line exchange mode determination

         If (j1 < mi1sc-2 .OR. j1 > mi1ec+2 .OR. &
             j2 < mi2sc-2 .OR. j2 > mi2ec+2) then
            idx_bound(3,mb) = 3
         Else If (j1 < mi1sc-1 .OR. j1 > mi1ec+1 .OR.   &
                  j2 < mi2sc-1 .OR. j2 > mi2ec+1) then
            idx_bound(3,mb) = 2
         Else
            idx_bound(3,mb) = 1
         End If

      End If

   End Do
End Do

mbsize = mb


!--- 2.3. Line exchange mode correction

Do j = 1, mbsize

   If (idx_bound(1,j) == -1 .AND. idx_bound(2,j) <= 2) then
      idx_bound(3,j) = 1 ! was previously set to 2
   End If

   If (idx_bound(1,j) == -2 .AND. idx_bound(2,j) <= 3) then
      idx_bound(3,j) = 2 ! was previously set to 3
   End If

End Do


!----------------------------------------------------------
! 2. SETTING PERMUTATION INDEX ARRAYS
!----------------------------------------------------------


!--- 2.0. Memory allocation

Allocate(idx_recv(3,knd*mbsize))
Allocate(idx_send(3,knd*mbsize))


!--- 2.1. Setting permutation index arrays

mb = 0

LineExchange: Do jr = 1,2
   Diamonds: Do jd = 1,knd

      Do jb = 1, mbsize

         j1 = idx_bound(1,jb)
         j2 = idx_bound(2,jb)

         If ((jr==1 .AND. idx_bound(3,jb)==1) .OR.   &
             (jr==2 .AND. idx_bound(3,jb)/=1)) then

            mb = mb+1

            idx_recv(1,mb) = j1
            idx_recv(2,mb) = j2
            idx_recv(3,mb) = jd

            idx_send(1,mb) = marr(2,j1,j2,jd)
            idx_send(2,mb) = marr(3,j1,j2,jd)
            idx_send(3,mb) = marr(4,j1,j2,jd)

         End If

      End Do

   End Do Diamonds

   If (jr==1) then
      np_1 = mb
   Else
      np_2 = mb
   End If

End Do LineExchange

np_max = np_2


!--- 2.2. Internal check

If (np_max /= knd*mbsize) then
   Write (0,*) 'set_boundaries: Internal error!'
   Stop
End If


End Subroutine Setup_Boundary



!==========================================================
Subroutine Set_Boundaries &
  (p,               & ! <-> Array to be enlarged on boundaries
   kip1s, kip1e,    & ! <-- First grid index limits
   kip2s, kip2e,    & ! <-- Second grid index limits
   kip3s, kip3e,    & ! <-- Component inex limits
   kr)                ! <-- Number of additional rows
!
! Enlarging field by adding rows on boundaries
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 29 Nov 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   kip1s, kip1e         ! First grid index limits
!
Integer, Intent(In) :: &
   kip2s, kip2e         ! Second grid index limits
!
Integer, Intent(In) :: &
   kip3s, kip3e         ! Component inex limits
                        ! Each component is enlarged independetly.
!
! Inout arguments:
!
Real(Double), Intent(InOut) :: &
   p(kip1s:kip1e, &      ! Array to be enlarged on boundaries
     kip2s:kip2e, &
     kip3s:kip3e, &
     nd)
!
Integer, Intent(In) :: &
   kr                   ! Number of additional rows
                        ! Can be 1 or 2.
!----------------------------------------------------------
! Local Scalars:
!
Integer :: j1, j2, jd   ! Grid indices
Integer :: j3           ! Component index
Integer :: j            ! Exchange buffer index
Integer :: mpr          ! Exchange buffer size
!
! Local Arrays:
!
Real(Double) :: &
   zbuff(np_max)        ! Exchange buffer
!----------------------------------------------------------
! Global variables used:
!
!    np_1           ! number of points to process - 1 boundary line
!    np_2           ! number of points to process - 2 boundary lines
!    np_max         ! max. number of points to process
!    idx_recv(:,:)  ! Indices of received points
!    idx_send(:,:)  ! Indices of points to send
!----------------------------------------------------------


!----------------------------------------------------------
! 1. BUFFERED EXCHANGE
!----------------------------------------------------------


!--- 1.1. Determination of buffer size

If (kr == 1) then
   mpr = np_1
Else
   mpr = np_2
End If


!--- 1.2. Exchange for each component

Do j3 = kip3s, kip3e

   Do j = 1, mpr
      j1 = idx_send(1,j)
      j2 = idx_send(2,j)
      jd = idx_send(3,j)
      zbuff(j) = p(j1,j2,j3,jd)
   End Do

   Do j = 1, mpr
      j1 = idx_recv(1,j)
      j2 = idx_recv(2,j)
      jd = idx_recv(3,j)
      p(j1,j2,j3,jd) = zbuff(j)
   End Do

End Do


End Subroutine Set_Boundaries



End Module ICO_boundary

