!
!+ GNSS Radio occultation observation operator: Error processing
!
MODULE Errors
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Error processing.
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
! Module Errors
!
! Error processing.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 20 Mar 1999 | Original code.
!----------------------------------------------------------
Implicit None
Private
Public :: Error_Status, Enter_Callee, Exit_Callee, Error, &
          Display_Status, Clear_Status
!----------------------------------------------------------
! Public Type Definitions:
!
Type Error_Status
   Character(Len=255)          :: Routine
   Integer                     :: User_Code
   Integer                     :: System_Code
   Character(Len=255)          :: Description
   Type(Error_Status), Pointer :: Caller_Status
   Type(Error_Status), Pointer :: Callee_Status
End Type Error_Status
!----------------------------------------------------------
!
Contains

!==========================================================
Subroutine Enter_Callee &
  (Routine,   & ! <-- User routine
   Stat)        ! <-> Pointer to callee status
!
! Adding callee error status to linked list.
!----------------------------------------------------------
! Method:
!   Linked lists.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 20 Mar 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Character(Len=*), Intent(In) :: &
   Routine     ! User routine name
!
! Inout arguments:
!
Type(Error_Status), Pointer :: &
   Stat        ! Input: Caller status (may be not associated)
               ! Ouput: Callee status linked to caller status
!----------------------------------------------------------
! Local Scalars:
!
Type(Error_Status), Pointer :: Callee_Status ! New callee status
!----------------------------------------------------------


!----------------------------------------------------------
! 1. MEMORY ALLOCATION AND FIELD INITIALIZATION
!----------------------------------------------------------

Allocate(Callee_Status)

Callee_Status%Routine        = Routine
Callee_Status%User_Code      = 0
Callee_Status%System_Code    = 0
Callee_Status%Description    = ''

Nullify(Callee_Status%Callee_Status)


!----------------------------------------------------------
! 2. LINKING
!----------------------------------------------------------

If (Associated(Stat)) then
   Callee_Status%Caller_Status => Stat
   Stat%Callee_Status          => Callee_Status
Else
   Nullify(Callee_Status%Caller_Status)
End If

Stat => Callee_Status


End Subroutine Enter_Callee



!==========================================================
Subroutine Exit_Callee &
  (Stat,         & ! <-> Pointer to callee status
   User_Code,    & ! <~~ User error code
   System_Code,  & ! <~~ System error code
   Description)    ! <~~ Error description
!
! Adding callee error status to linked list.
!----------------------------------------------------------
! Method:
!   Linked lists.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 20 Mar 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Inout arguments:
!
Type(Error_Status), Pointer :: &
   Stat           ! Input: Callee status linked to caller status.
                  ! Ouput: Caller status.
                  !        Callee status is nullified if no error
                  !        occurred.
!
! Input arguments:
!
Integer, Optional, Intent(In) :: &
   User_Code      ! User error code
!
Integer, Optional, Intent(In) :: &
   System_Code    ! System error code
!
Character(Len=*), Optional, Intent(In) :: &
   Description    ! Error description
!----------------------------------------------------------

If (.not. Associated(Stat)) then
   Return
End If

If ((.not. Present(User_Code))   .and.  &
    (.not. Present(System_Code)) .and.  &
    (.not. Associated(Stat%Callee_Status))) then
   Stat => Stat%Caller_Status
   Deallocate(Stat%Callee_Status)
Else
   If (Present(User_Code)) then
      Stat%User_Code   = User_Code
   End If
   If (Present(System_Code)) then
      Stat%System_Code = System_Code
   End If
   If (Present(Description)) then
      Stat%Description = Description
   End If
   Stat => Stat%Caller_Status
End If


End Subroutine Exit_Callee



!==========================================================
Function Error &
  (Stat)     ! <-- Pointer to status
!
! Checking callee error status.
!----------------------------------------------------------
! Method:
!   Linked lists.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 20 Mar 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Error_Status), Pointer :: &
   Stat           ! Caller status.
!
! Function result:
!
Logical :: &
   Error    ! .True. if callee error occurred.
!----------------------------------------------------------


Error = Associated(Stat%Callee_Status)


End Function Error



!==========================================================
Subroutine Clear_Status &
  (Stat)    ! --> Status to clear
!
! Clearing error codes and deallocation of linked chain
! of callee statuses.
!----------------------------------------------------------
! Method:
!   Deallocation of linked lists.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 04 Jan 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Inout arguments:
!
Type(Error_Status), Pointer :: &
   Stat     ! Status to clear
!----------------------------------------------------------
! Local Scalars:
!
Type(Error_Status), Pointer :: &
   Current,   &  ! Current status in chain
   Next          ! Next status in chain
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CHAIN TRACING
!----------------------------------------------------------

If (.not. Associated(Stat%Callee_Status)) then
   Return
End If

Current => Stat%Callee_Status


Clear: Do
   If (.not. Associated(Current%Callee_Status)) then
      Exit Clear
   Else
      Next => Current%Callee_Status
   End If
   Deallocate(Current)
   Current => Next
End Do Clear


End Subroutine Clear_Status



!==========================================================
Subroutine Display_Status &
  (Stat)    ! --> Status to display
!
! Displaying callee statuses.
!----------------------------------------------------------
! Method:
!   Linked lists.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 04 Jan 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Inout arguments:
!
Type(Error_Status), Pointer :: &
   Stat     ! Status to print
!----------------------------------------------------------
! Local Scalars:
!
Type(Error_Status), Pointer :: &
   Current,   &  ! Current status in chain
   Next          ! Next status in chain
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CHAIN TRACING
!----------------------------------------------------------

If (.not. Associated(Stat%Callee_Status)) then
   Return
End If

Current => Stat%Callee_Status

Display: Do
   Write (0,'(A,": ",A/,"  User: ",I5,"  System: ",I5)')    &
      Trim(Current%Routine), Trim(Current%Description),     &
      Current%User_Code, Current%System_Code
   If (.not. Associated(Current%Callee_Status)) then
      Exit Display
   Else
      Next => Current%Callee_Status
   End If
   Current => Next
End Do Display


End Subroutine Display_Status



End Module Errors


