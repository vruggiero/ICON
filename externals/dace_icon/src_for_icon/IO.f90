!
!+ GNSS Radio occultation observation operator: I/O routines
!
MODULE IO
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   I/O routines.
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
! Module IO
!
! I/O routines.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Dec 1997 | Original code.
!   1.1   | 10 Sep 1998 | Modified FileOpen
!   2.0   | 01 Oct 1998 | Split_Pathname
!   2.1   | 01 Oct 1998 | FileOpen parameters same as Open
!   3.0   ! 15 Oct 1998 | PutXY and DeleteFile
!   4.0   | 18 Jan 1999 | Single and Double versions of PutXY.
!----------------------------------------------------------
Implicit None
Private
Public :: Fileopen
!----------------------------------------------------------
! Interfaces:
!
Interface PutXY
   Module Procedure PutXY_Double
   Module Procedure PutXY_Single
End Interface
!----------------------------------------------------------
Contains

!==========================================================
Subroutine FileOpen(  &
   Name,     & ! <-- File pathname
   Status,   & ! <~~ Open status
   Access,   & ! <~~ Access mode
   Form,     & ! <~~ File form
   RecL,     & ! <~~ Record length
   Position, & ! <~~ Open position
   Action,   & ! <~~ Read/write permissions
   Unit,     & ! --> Unit number
   IOStat)     ! --> Error code
!
! Opening a file with given options.
!----------------------------------------------------------
! Method:
!   Search for a free unit, opening the file,
!   checking for possible errors (see also
!   Fortran-90 Reference Manual).
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Dec 1997 | Original code.
!   2.0   | 10 Sep 1998 | Subroutine, more options.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
! (Default values are marked with *)
!
Character(Len=*), Intent(In)           :: &
   Name     ! The pathname of the file to open

Character(Len=*), Intent(In), Optional :: &
   Status   ! File status at opening
            ! (OLD, NEW, REPLACE, SCRATCH, UNKNOWN*)

Character(Len=*), Intent(In), Optional :: &
   Access   ! File access mode
            ! (SEQUENTIAL*, DIRECT)

Character(Len=*), Intent(In), Optional :: &
   Form     ! File form
            ! (FORMATTED   default if Access==SEQUENTIAL,
            !  UNFORMATTED default if Access==DIRECT)

Integer,          Intent(In), Optional :: &
   RecL     ! File record length

Character(Len=*), Intent(In), Optional :: &
   Position ! Initial position
            ! (ASIS*, REWIND, APPEND)

Character(Len=*), Intent(In), Optional :: &
   Action   ! File action
            ! (READ, WRITE, READWRITE*)
!
! Output arguments:
!
Integer, Intent(Out)                   :: &
   Unit     ! File unit number:
            !   >0: A correct unit number
            !    0: No free units
            !   -1: Error opening file (see IOStat)

Integer, Intent(Out)                   :: &
   IOStat   ! Error code:
            !    0: No IO error (but maybe no free units)
            !   >0: System-dependent error code
!----------------------------------------------------------
! Local Scalars:
!
Integer           :: i         ! Unit number tested
Logical           :: Busy      ! Unit status
Character(Len=20) :: FStatus   ! Variable for Status
Character(Len=20) :: FAccess   ! Variable for Access
Character(Len=20) :: FForm     ! Variable for Form
Character(Len=20) :: FPosition ! Variable for Position
Character(Len=20) :: FAction   ! Variable for Action
!----------------------------------------------------------


!----------------------------------------------------------
! 1. SEARCHING FOR FREE UNIT
!----------------------------------------------------------

Do i=10,99
   Unit = i
   Inquire(Unit = Unit, Opened = Busy, IOStat = IOStat)
   If ((IOStat == 0) .and. (.not. Busy)) then
      Exit
   End If
End Do

If (Busy .or. (IOStat /= 0)) then
  Unit = 0
  Return
End If


!----------------------------------------------------------
! 2. SETTINGS OPTIONS
!----------------------------------------------------------

If (Present(Status)) then
   FStatus = Status
Else
   FStatus = 'UNKNOWN'
End If

If (Present(Access)) then
   FAccess = Access
Else
   FAccess = 'SEQUENTIAL'
End If

If (Present(Form)) then
   FForm = Form
Else
   If (FAccess == 'SEQUENTIAL') then
      FForm = 'FORMATTED'
   Else
      FForm = 'UNFORMATTED'
   End If
End If

If (Present(Position)) then
   FPosition = Position
Else
   FPosition = 'ASIS'
End If

If (Present(Action)) then
   FAction = Action
Else
   FAction = 'READWRITE'
End If

!----------------------------------------------------------
! 3. OPENING THE FILE
!----------------------------------------------------------

If (Present(RecL)) then
   Open(Unit     = Unit,      &
        File     = Name,      &
        Status   = FStatus,   &
        Access   = FAccess,   &
        Form     = FForm,     &
        RecL     = RecL,      &
        Position = FPosition, &
        Action   = FAction,   &
        IOStat   = IOStat)
Else
   Open(Unit     = Unit,      &
        File     = Name,      &
        Status   = FStatus,   &
        Access   = FAccess,   &
        Form     = FForm,     &
        Position = FPosition, &
        Action   = FAction,   &
        IOStat   = IOStat)
End If

Return


End Subroutine FileOpen



!==========================================================
Subroutine Split_Pathname &
  (Pathname,   & ! <-- Pathname to split
   Path,       & ! --> Path
   Filename)     ! --> File name
!
! Splitting of a pathname into path and filename.
!----------------------------------------------------------
! Method:
!   Search for last position of path delimiter.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 01 Oct 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Character(Len=*), Intent(In) :: &
   Pathname  ! Pathname to split
!
! Inout arguments:
!
Character(Len=*), Intent(InOut) :: &
   Path      ! Path
!
! Inout arguments:
!
Character(Len=*), Intent(InOut) :: &
   Filename  ! File name
!----------------------------------------------------------
! Local Scalars:
!
Integer :: Pos  ! Last delimiter position in the pathname
!----------------------------------------------------------


Pos = Scan(Pathname, '/', .True.)

If (Pos /= 0) then
   Path     = Pathname(1:Pos)
   Filename = Pathname(Pos+1:)
Else
   Path     = './'
   Filename = Pathname
End If


End Subroutine Split_Pathname



!==========================================================
Subroutine PutXY_Double &
  (Name,  & ! <-- Pathname of output file
   X,     & ! <-- Array of X-cooridnates
   Y1,    & ! <~~ Array of Y1(X)-function
   Y2,    & ! <~~ Array of Y2(X)-function
   Y3,    & ! <~~ Array of Y3(X)-function
   Y4,    & ! <~~ Array of Y4(X)-function
   Y5,    & ! <~~ Array of Y5(X)-function
   Y6,    & ! <~~ Array of Y5(X)-function
   XFmt,  & ! <~~ Output format of X data
   YFmt,  & ! <~~ Output format of Y data
   Stat)    ! --> Error code
!
! Writing a data file in X-Y1-YN format.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 15 Oct 1998 | Original code.
!   2.0   | 25 Oct 1998 | All Yi optional, Y6.
!   3.0   | 18 Jan 1999 | Single and Double versions.
!   3.1   | 02 Oct 1999 | Delimiters between output fields.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Character(Len=*), Intent(In)       :: &
   Name     ! Output file pathname
!
Real(Double), Intent(In)           :: &
   X(1:)    ! Array of X-cooridnates
!
Real(Double), Optional, Intent(In) :: &
   Y1(1:)    ! Array of Y1(X)-function
!
Real(Double), Optional, Intent(In) :: &
   Y2(1:)    ! Array of Y2(X)-function
!
Real(Double), Optional, Intent(In) :: &
   Y3(1:)    ! Array of Y3(X)-function
!
Real(Double), Optional, Intent(In) :: &
   Y4(1:)    ! Array of Y4(X)-function
!
Real(Double), Optional, Intent(In) :: &
   Y5(1:)    ! Array of Y5(X)-function
!
Real(Double), Optional, Intent(In) :: &
   Y6(1:)    ! Array of Y6(X)-function
!
Character(Len=*), Optional, Intent(In)       :: &
   XFmt      ! Output format of X data
!
Character(Len=*), Optional, Intent(In)       :: &
   YFmt      ! Output format of Y data
!
! Output arguments:
!
Integer, Intent(Out) :: &
   Stat      ! Error Code
             !    0 - no error
             !   >0 - error code from Open, Write, or Close
             !   <0 - array size mismatch
!----------------------------------------------------------
! Local Scalars:
!
Integer           :: N     ! Number of X-grid points
Integer           :: i     ! Y-array number
Integer           :: k     ! Grid index
Integer           :: IU    ! File unit number
Character(Len=80) :: XAFmt ! Actual output format of X data
Character(Len=80) :: YAFmt ! Actual output format of Y data
!
! Local Arrays:
!
Target :: &
   Y1, Y2, Y3, Y4, Y5, Y6
!
Type ArrayPointer
   Real(Double), Pointer :: Y(:)  ! Y-array pointer
End Type ArrayPointer
!
Type(ArrayPointer) :: &
   P(6)        ! Y-Array pointers
!
Logical :: &
   YPresent(6) ! Indicator of presense
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION, ARRAY SIZE CHECK
!----------------------------------------------------------

N = Size(X)

YPresent = (/ Present(Y1), Present(Y2), &
              Present(Y3), Present(Y4), &
              Present(Y5), Present(Y6)  /)

Do i=1,6
   If (YPresent(i)) then
      Select Case(i)
         Case(1)
            P(i)%Y => Y1
         Case(2)
            P(i)%Y => Y2
         Case(3)
            P(i)%Y => Y3
         Case(4)
            P(i)%Y => Y4
         Case(5)
            P(i)%Y => Y5
         Case(6)
            P(i)%Y => Y6
      End Select
      If (Size(P(i)%Y) /= N) then
         Stat = -i
         Return
      End If
   End If
End Do

If (Present(XFmt)) then
   XAFmt = XFmt
Else
   XAFmt = '(ES23.14E3)'
End If

If (Present(YFmt)) then
   YAFmt = YFmt
Else
   YAFmt = '(ES23.14E3)'
End If

!----------------------------------------------------------
! 2. OPENING FILE
!----------------------------------------------------------

Call FileOpen(  &
   Name,                    & ! <-- File pathname
   Access   = 'SEQUENTIAL', & ! <~~ Access mode
   Form     = 'FORMATTED',  & ! <~~ File form
   RecL     = 256,          & ! <~~ Record length
   Position = 'REWIND',     & ! <~~ Open position
   Action   = 'WRITE',      & ! <~~ Read/write permissions
   Unit     = IU,           & ! --> Unit number
   IOStat   = Stat)           ! --> Error code

If (Stat /= 0) then
   Return
End If


!----------------------------------------------------------
! 3. WRITING OUTPUT DATA
!----------------------------------------------------------


Do k=1,N
   Write (IU, XAFmt, Advance='NO',IOStat=Stat) X(k)
   If (Stat /= 0) then
      Close (Unit = IU)
      Return
   End If
   Do i=1,6
      If (YPresent(i)) then
         Write (IU, '(1X)', Advance='NO')
         Write (IU, YAFmt, Advance='NO',IOStat=Stat) P(i)%Y(k)
         If (Stat /= 0) then
            Close (Unit = IU)
            Return
         End If
      End If
   End Do
   Write(IU, '()')
End Do

Close (Unit = IU)


End Subroutine PutXY_Double



!==========================================================
Subroutine PutXY_Single &
  (Name,  & ! <-- Pathname of output file
   X,     & ! <-- Array of X-cooridnates
   Y1,    & ! <~~ Array of Y1(X)-function
   Y2,    & ! <~~ Array of Y2(X)-function
   Y3,    & ! <~~ Array of Y3(X)-function
   Y4,    & ! <~~ Array of Y4(X)-function
   Y5,    & ! <~~ Array of Y5(X)-function
   Y6,    & ! <~~ Array of Y5(X)-function
   XFmt,  & ! <~~ Output format of X data
   YFmt,  & ! <~~ Output format of Y data
   Stat)    ! --> Error code
!
! Writing a data file in X-Y1-YN format.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 15 Oct 1998 | Original code.
!   2.0   | 25 Oct 1998 | All Yi optional, Y6.
!   3.0   | 18 Jan 1999 | Single and Double versions.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Single
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Character(Len=*), Intent(In)       :: &
   Name     ! Output file pathname
!
Real(Single), Intent(In)           :: &
   X(1:)    ! Array of X-cooridnates
!
Real(Single), Optional, Intent(In) :: &
   Y1(1:)    ! Array of Y1(X)-function
!
Real(Single), Optional, Intent(In) :: &
   Y2(1:)    ! Array of Y2(X)-function
!
Real(Single), Optional, Intent(In) :: &
   Y3(1:)    ! Array of Y3(X)-function
!
Real(Single), Optional, Intent(In) :: &
   Y4(1:)    ! Array of Y4(X)-function
!
Real(Single), Optional, Intent(In) :: &
   Y5(1:)    ! Array of Y5(X)-function
!
Real(Single), Optional, Intent(In) :: &
   Y6(1:)    ! Array of Y6(X)-function
!
Character(Len=*), Optional, Intent(In)       :: &
   XFmt      ! Output format of X data
!
Character(Len=*), Optional, Intent(In)       :: &
   YFmt      ! Output format of Y data
!
! Output arguments:
!
Integer, Intent(Out) :: &
   Stat      ! Error Code
             !    0 - no error
             !   >0 - error code from Open, Write, or Close
             !   <0 - array size mismatch
!----------------------------------------------------------
! Local Scalars:
!
Integer           :: N     ! Number of X-grid points
Integer           :: i     ! Y-array number
Integer           :: k     ! Grid index
Integer           :: IU    ! File unit number
Character(Len=80) :: XAFmt ! Actual output format of X data
Character(Len=80) :: YAFmt ! Actual output format of Y data
!
! Local Arrays:
!
Target :: &
   Y1, Y2, Y3, Y4, Y5, Y6
!
Type ArrayPointer
   Real(Single), Pointer :: Y(:)  ! Y-array pointer
End Type ArrayPointer
!
Type(ArrayPointer) :: &
   P(6)        ! Y-Array pointers
!
Logical :: &
   YPresent(6) ! Indicator of presense
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION, ARRAY SIZE CHECK
!----------------------------------------------------------

N = Size(X)

YPresent = (/ Present(Y1), Present(Y2), &
              Present(Y3), Present(Y4), &
              Present(Y5), Present(Y6)  /)

Do i=1,6
   If (YPresent(i)) then
      Select Case(i)
         Case(1)
            P(i)%Y => Y1
         Case(2)
            P(i)%Y => Y2
         Case(3)
            P(i)%Y => Y3
         Case(4)
            P(i)%Y => Y4
         Case(5)
            P(i)%Y => Y5
         Case(6)
            P(i)%Y => Y6
      End Select
      If (Size(P(i)%Y) /= N) then
         Stat = -i
         Return
      End If
   End If
End Do

If (Present(XFmt)) then
   XAFmt = XFmt
Else
   XAFmt = '(ES23.14E3)'
End If

If (Present(YFmt)) then
   YAFmt = YFmt
Else
   YAFmt = '(ES23.14E3)'
End If

!----------------------------------------------------------
! 2. OPENING FILE
!----------------------------------------------------------

Call FileOpen(  &
   Name,                    & ! <-- File pathname
   Access   = 'SEQUENTIAL', & ! <~~ Access mode
   Form     = 'FORMATTED',  & ! <~~ File form
   RecL     = 256,          & ! <~~ Record length
   Position = 'REWIND',     & ! <~~ Open position
   Action   = 'WRITE',      & ! <~~ Read/write permissions
   Unit     = IU,           & ! --> Unit number
   IOStat   = Stat)           ! --> Error code

If (Stat /= 0) then
   Return
End If


!----------------------------------------------------------
! 3. WRITING OUTPUT DATA
!----------------------------------------------------------


Do k=1,N
   Write (IU, XAFmt, Advance='NO',IOStat=Stat) X(k)
   If (Stat /= 0) then
      Close (Unit = IU)
      Return
   End If
   Do i=1,6
      If (YPresent(i)) then
         Write (IU, '(1X)', Advance='NO')
         Write (IU, YAFmt, Advance='NO',IOStat=Stat) P(i)%Y(k)
         If (Stat /= 0) then
            Close (Unit = IU)
            Return
         End If
      End If
   End Do
   Write(IU, '()')
End Do

Close (Unit = IU)


End Subroutine PutXY_Single



!==========================================================
Subroutine DeleteFile &
  (Name,     & ! <-- File name
   IOStat)     ! --> Error code
!
! Deleting a file.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 15 Oct 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Character(Len=*), Intent(In) :: &
   Name     ! File name
!
! Output arguments:
!
Integer, Intent(Out) :: &
   IOStat   ! Error code from FileOpen or Close
!----------------------------------------------------------
! Local Scalars:
!
Integer :: IU  ! File unit number
!----------------------------------------------------------


Call FileOpen(Name, Unit=IU, IOStat=IOStat)
If (IOStat /= 0) then
   Return
End If

Close (Unit = IU, Status='DELETE', IOStat=IOStat)


End Subroutine DeleteFile



End Module IO

