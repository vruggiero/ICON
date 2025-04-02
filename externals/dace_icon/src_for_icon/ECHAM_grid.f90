!
!+ GNSS Radio occultation observation operator: interpolation coefficients
!
MODULE ECHAM_grid
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Calculation of ECHAM vertical and horizontal grids
!   and interpolation coefficients.
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
! V1_6         2009/06/10 Harald Anlauf
!  Fix precision of FP constants
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_22        2013-02-13 Harald Anlauf
!  add vertical coordinate type, geopotential height
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
! Module ECHAM_grid
!
! Calculation of ECHAM vertical and horizontal grids
! and interpolation coefficients.
!----------------------------------------------------------
! (C) Copyright 1998-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 16 Dec 1998 | Original code.
!   2.0   | 20 Mar 1999 | Error check.
!   3.0   | 07 Apr 2002 | Latitude_Interpolation,
!         |             | Longitude_Interpolation, Weight4.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Single, Double, WorkPr, &
    Pi, rtd
!
Use Errors, only: &
! Imported Type Definitions:
    Error_Status,  &
! Imported Routines:
    Enter_Callee,  &
    Exit_Callee
!----------------------------------------------------------
Implicit None
Private
Public :: longitude_interpolation, latitude_interpolation
!----------------------------------------------------------
! Public Parameters:
!
Integer, Parameter :: &
   err_AB_Size  = 5101, &
   err_No_AB    = 5102, &
   err_Bad_Type = 5103
!----------------------------------------------------------
!
Contains

!==========================================================
Subroutine Vertical_Coordinates &
  (A,     & ! --> A(k+1/2) coefficients
   B,     & ! --> B(k+1/2) coefficients
   ErrStat) ! --> Error status
!
! Calculation of vertical coordinates A(i) and B(i):
! P(i) = A(i) + B(i)*Psur, P and Psur [Pa].
!----------------------------------------------------------
! Method:
!   Ready vertical cooridate set for given
!   number of vertical levels.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 05 Oct 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Output arguments:
!
Real(Double), Intent(Out) :: &
   A(0:)    ! A(k+1/2) coefficients
!
Real(Double), Intent(Out) :: &
   B(0:)    ! B(k+1/2) coefficients
!
Type(Error_Status), Pointer :: &
   ErrStat  ! Error status
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   A19(0:19) = &  ! A(k+1/2) for 19 levels
     (/   0.000000_Double,      2000.000000_Double,  &
       4000.000000_Double,      6046.110595_Double,  &
       8267.927560_Double,     10609.513232_Double,  &
      12851.100169_Double,     14698.498086_Double,  &
      15861.125180_Double,     16116.236610_Double,  &
      15356.924115_Double,     13621.460403_Double,  &
      11101.561987_Double,      8127.144155_Double,  &
       5125.141747_Double,      2549.969411_Double,  &
        783.195032_Double,         0.000000_Double,  &
          0.000000_Double,         0.000000_Double /)
!
Real(Double), Parameter :: &
   B19(0:19) = &  ! B(k+1/2) for 19 levels
     (/   0.0000000000_Double,     0.0000000000_Double,  &
          0.0000000000_Double,     0.0003389933_Double,  &
          0.0033571866_Double,     0.0130700434_Double,  &
          0.0340771467_Double,     0.0706498323_Double,  &
          0.1259166826_Double,     0.2011954093_Double,  &
          0.2955196487_Double,     0.4054091989_Double,  &
          0.5249322235_Double,     0.6461079479_Double,  &
          0.7596983769_Double,     0.8564375573_Double,  &
          0.9287469142_Double,     0.9729851852_Double,  &
          0.9922814815_Double,     1.0000000000_Double /)
!
!
Real(Double), Parameter :: &
   A31(0:31) = &  ! A(k+1/2) for 31 levels
     (/   0.000000000000_Double,      2000.000000000000_Double,  &
       4000.000000000000_Double,      6000.000000000000_Double,  &
       8000.000000000000_Double,      9976.136718750000_Double,  &
      11820.539062500000_Double,     13431.394531250000_Double,  &
      14736.355468750000_Double,     15689.207031250000_Double,  &
      16266.609375000000_Double,     16465.003906250000_Double,  &
      16297.621093750000_Double,     15791.597656250000_Double,  &
      14985.269531250000_Double,     13925.519531250000_Double,  &
      12665.292968750000_Double,     11261.230468750000_Double,  &
       9771.406250000000_Double,      8253.210937500000_Double,  &
       6761.339843750000_Double,      5345.914062500000_Double,  &
       4050.717773437500_Double,      2911.569335937500_Double,  &
       1954.805175781250_Double,      1195.889892578125_Double,  &
        638.148925781250_Double,       271.626464843750_Double,  &
         72.063583374023_Double,         0.000000000000_Double,  &
          0.000000000000_Double,         0.000000000000_Double  /)
!
Real(Double), Parameter :: &
   B31(0:31) = &  ! B(k+1/2) for 31 levels
     (/   0.000000000000_Double,         0.000000000000_Double,  &
          0.000000000000_Double,         0.000000000000_Double,  &
          0.000000000000_Double,         0.000390858157_Double,  &
          0.002919700695_Double,         0.009194131941_Double,  &
          0.020319156349_Double,         0.036974858493_Double,  &
          0.059487640858_Double,         0.087894976139_Double,  &
          0.122003614902_Double,         0.161441504955_Double,  &
          0.205703258514_Double,         0.254188597202_Double,  &
          0.306235373020_Double,         0.361145019531_Double,  &
          0.418202280998_Double,         0.476688146591_Double,  &
          0.535886585712_Double,         0.595084249973_Double,  &
          0.653564572334_Double,         0.710594415665_Double,  &
          0.765405237675_Double,         0.817166984081_Double,  &
          0.864955842495_Double,         0.907715857029_Double,  &
          0.944213211536_Double,         0.972985208035_Double,  &
          0.992281496525_Double,         1.000000000000_Double  /)
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------

Call Enter_Callee &
  ('Vertical_Coordinates',  & ! <-- User routine
   ErrStat)                   ! <-> Pointer to callee status


!----------------------------------------------------------
! 1. ARRAY SIZE CHECK
!----------------------------------------------------------

If (Size(A) /= Size(B)) then
   Call Exit_Callee &
     (ErrStat,           & ! <-> Pointer to callee status
      err_AB_Size,       & ! <~~ User error code
      0,                 & ! <~~ System error code
      'A/B size mismatch') ! <~~ Error description
   Return
End If


!----------------------------------------------------------
! 2. CHOICE OF COEFFICIENT SET
!----------------------------------------------------------

Select Case (UBound(A,1))
   Case(19)
      A(:) = A19(:)
      B(:) = B19(:)
   Case(31)
      A(:) = A31(:)
      B(:) = B31(:)
   Case Default
      Call Exit_Callee &
        (ErrStat,     & ! <-> Pointer to callee status
         err_No_AB,   & ! <~~ User error code
         0,           & ! <~~ System error code
         'No A/B')      ! <~~ Error description
      Return
End Select


Call Exit_Callee(ErrStat)


End Subroutine Vertical_Coordinates



!==========================================================
Subroutine Rectangular_Grid &
  (Grid_Type, & ! <-- Data representation type
   XLon,      & ! --> Longitude grid [deg]
   XLat,      & ! --> Latitude grid [deg]
   ErrStat)     ! --> Error status
!
! Generation of surface grid.
!----------------------------------------------------------
! Method:
!   Calculation of homogeneous lat/lon or gaussian grid.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 16 Dec 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   Grid_Type ! Data representation type
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   XLon(1:)  ! Longitude grid [deg]
!
Real(Double), Intent(Out) :: &
   XLat(1:)  ! Latitude grid [deg]
!
Type(Error_Status), Pointer :: &
   ErrStat   ! Error status
!----------------------------------------------------------
! Local Parameters:
!
Integer, Parameter :: &
   grd_LatLon   = 0,  &
   grd_Gaussian = 4
!
! Local Scalars:
!
Integer :: NLon  ! Longitude dimension
Integer :: NLat  ! Latitude dimension
Integer :: i     ! Array index
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------

Call Enter_Callee &
  ('Vertical_Coordinates',  & ! <-- User routine
   ErrStat)                   ! <-> Pointer to callee status


!----------------------------------------------------------
! 1. CALCULATION OF GRID DIMENSIONS
!----------------------------------------------------------

NLon = Size(XLon)
NLat = Size(XLat)

!----------------------------------------------------------
! 2. CALCULATION OF GRID
!----------------------------------------------------------

Select Case (Grid_Type)

   Case (grd_LatLon)
      Do i=1,NLon
         XLon(i) = 360.0_Double*(i-1)/NLon
      End Do
      Do i=1,NLat
         XLat(i) = 90.0_Double - 180.0_Double*(i-1)/(NLat-1)
      End Do

   Case (grd_Gaussian)
      Do i=1,NLon
         XLon(i) = 360.0_Double*(i-1)/NLon
      End Do
      Call GaussGrid(XLat)

   Case Default
      Call Exit_Callee &
        (ErrStat,            & ! <-> Pointer to callee status
         err_Bad_Type,       & ! <~~ User error code
         0,                  & ! <~~ System error code
         'Unknown data type')  ! <~~ Error description
      Return

End Select


Call Exit_Callee(ErrStat)


End Subroutine Rectangular_Grid



!==========================================================
Subroutine GaussGrid &
  (Gauss)  ! --> Gaussian grid of latitudes [deg]
!
!  Compute gaussian grid of latitudes.
!----------------------------------------------------------
! Method:
!   Search for zeroes of Legendre polynomial.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!         |    Sep 1996 | Fortran-77: F. VAN DEN BERGHE
!   1.0   | 19 Dec 1997 | Fortran-90: M. E. Gorbunov
!   1.1   | 26 Dec 1997 | Simplifications of code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Output arguments:
!
Real(Double), Intent(Out) :: &
   Gauss(1:) ! Gaussian grid of latitudes [deg]
!----------------------------------------------------------
! Local Parameters:
Real(Double), Parameter :: &
   Eps = 1.0e-12_Double ! Accuracy
!
! Local Scalars:
!
Integer      :: K2      ! Full grid dimension
Integer      :: KHalf   ! Halg grid dimension
Integer      :: K       ! Grid index
Real(Double) :: Theta   ! Polar latitude (colatitude) [rad]
Real(Double) :: Step    ! Step for seek of zeros
Real(Double) :: P1, P2  ! Values of Legendre polynomial
!----------------------------------------------------------


KHalf = Size(Gauss)/2
K2    = 2*KHalf
Theta = 0

Do K=1,KHalf

   Step = Pi/(8*KHalf)

   Seek_Zero: Do

      P2 = Legendre(K2, Theta)

      Seek_Sign: Do
         P1    = P2
         Theta = Theta + Step
         P2 = Legendre(K2, Theta)
         If (Sign(1.0_Double, P1) /= Sign(1.0_Double, P2)) then
            Exit Seek_Sign
         End If
      End Do Seek_Sign

      If (Step < Eps) then
         Exit Seek_Zero
      End If

      Theta = Theta - Step
      Step  = Step/2.71828_Double

   End Do Seek_Zero

   Gauss(K)      = 90.0_Double - Theta*rtd
   Gauss(K2-K+1) = - Gauss(K)

End Do


End Subroutine GaussGrid



!==========================================================
Function Legendre  &
  (N,      & ! <-- Degree of Legendre polynomial
   Theta)  & ! <-- Colatitude [rad]
Result(P)    ! --> Value of Legendre polynomial
!
! Calculation of the Legendre polynomial.
!----------------------------------------------------------
! Method:
!   Recurrent formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!         |    Sep 1996 | Fortran-77. F. VAN DEN BERGHE
!   1.0   | 19 Dec 1997 | Fortran-90. M. E. Gorbunov
!   1.1   | 26 Dec 1997 | Simplifications of code.
!   2.0   | 16 Dec 1998 | Function
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In)       :: &
   N        ! Degree of Legendre polynomial
!
Real(Double), Intent(In)  :: &
   Theta    ! Colatitude [rad]
!
! Function result:
!
Real(Double)              :: &
   P        ! Value of Legendre polynomial
!----------------------------------------------------------
! Local Scalars:
Real(Double) :: X, Y1, Y2, Y3, G
Integer      :: i
!----------------------------------------------------------


X  = Cos(Theta)
Y1 = 1.0_Double
Y2 = X

Do i=2,N
   G  = X*Y2
   Y3 = G - Y1 + G - (G - Y1)/Real(i, Double)
   Y1 = Y2
   Y2 = Y3
End Do

P = Y3


End Function Legendre



!==========================================================
Subroutine Longitude_Interpolation &
  (XLon,      & ! <-- Longitude grid [deg]
   PLon,      & ! <-- Longiude of point [deg]
   ILon,      & ! --> Interpolation subgrid
   WLon,      & ! --> Weights of subgrid points
   WLon1,     & ! ~~> Weight derivatives
   WLon2)       ! ~~> Weight 2nd derivatives
!
! Longitude interpolation: calculation of interpolation
! weights of grid points and, optionally, the 1st and 2nd
! derivatives of weighting function.
!----------------------------------------------------------
! Method:
!   For N-point scheme, the central cell containing
!   interpolation point is found. Coefficients of
!   polynomial of N-1 power are found so as for
!   derivatives up to (N-2)th to be continuous.
!----------------------------------------------------------
! (C) Copyright 1998-99, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 13 Nov 1998 | Original code.
!   2.0   ! 20 Feb 1999 ! Faster index search.
!----------------------------------------------------------
! Modules used:
!
Use Interpolation, only: &
! Imported Routines:
    Seek_Index
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Workpr), Intent(In) :: &
   XLon(1:)   ! Longitude grid [deg]
!
Real(Double), Intent(In) :: &
   PLon       ! Longiude of point [deg]
!
! Output arguments:
!
Integer, Intent(Out) :: &
   ILon(1:)   ! Indices of grid points used
              ! Size(ILon) must be even
!
Real(Double), Intent(Out) :: &
   WLon(1:)   ! Weights of grid points used
              ! Size(WLon) must be = Size(ILon)
!
Real(Double), Intent(Out), Optional :: &
   WLon1(1:)  ! Weight derivatives
              ! Size(WLon1) must be = Size(ILon)
!
Real(Double), Intent(Out), Optional :: &
   WLon2(1:)  ! Weight 2nd derivatives
              ! Size(WLon2) must be = Size(ILon)
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: NLon  ! Longitude grid dimension
Integer      :: DLon  ! Longitude grid direction
Integer      :: NI    ! Interpolation scheme dimension
Integer      :: i     ! Grid index
Integer      :: IC    ! Lower index of central cell of subgrid
Real(Double) :: TLon  ! Longitude in 0-360 deg interval
!
! Local Arrays:
!
Real(Double) :: &
   SLon(Size(ILon))    ! Interpolation subgrid
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF GRID SIZE AND DIRECTION
!----------------------------------------------------------


!--- 1.1. Main grid

NLon = Size(XLon)
DLon = Nint(Sign(1.0_Workpr, XLon(NLon)-XLon(1)))


!--- 1.2. Interpolation subgrid

NI   = Size(ILon)


!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATION CELL
!----------------------------------------------------------


!--- 2.1. Transform of longitude to 0-360 interval

TLon = XLon(1) + Modulo(PLon-XLon(1), 360.0_Double)


!--- 2.2. Calculation of central cell of interpolation subgrid

IC       = NI/2

!ILon(IC) = Sum(MaxLoc(DLon*Xlon(:), &
!                      DLon*XLon(:) <= DLon*TLon))
If (DLon*TLon >= DLon*XLon(NLon)) then
   ILon(IC) = NLon
Else
   ILon(IC) = Seek_Index (XLon, Real(TLon, Workpr))
End If
ILon(IC) = Max(1, ILon(IC))
SLon(IC) = XLon(ILon(IC))


!--- 2.3. Calculation of interpolation subgrid

Do i=1,NI
   ILon(i) = 1  + Modulo(ILon(IC)+i-IC-1, NLon)
   SLon(i) = SLon(IC)-180.0_Double + &
      Modulo(XLon(ILon(i))-SLon(IC)+180.0_Double, 360.0_Double)
End Do


!----------------------------------------------------------
! 3. CALCULATION OF INTERPOLATION WEIGHTS
!----------------------------------------------------------


Select Case (NI)


!--- 3.1. 2-point linear interpolation

Case(2)

   WLon(1) = (SLon(2) - TLon)/(SLon(2)-SLon(1))
   WLon(2) = (TLon - SLon(1))/(SLon(2)-SLon(1))

   If (Present(WLon1)) then
      WLon1(:) = 0
   End If

   If (Present(WLon2)) then
      WLon2(:) = 0
   End If


!--- 3.2. 4-point linear interpolation

Case(4)

   Call Weight4 &
     (TLon,      & ! <-- Interpolation point
      SLon,      & ! <-- Interpolation subgrid
      0,         & ! <-- Border indicator
      WLon,      & ! --> Interpolation weights
      WLon1)       ! ~~> Interpolation weight derivatives

   If (Present(WLon2)) then
      WLon2(:) = 0
   End If


!--- 3.3. Default

Case Default

   WLon(:) = 0

   If (Present(WLon1)) then
      WLon1(:) = 0
   End If

   If (Present(WLon2)) then
      WLon2(:) = 0
   End If


End Select


End Subroutine Longitude_Interpolation



!==========================================================
Subroutine Latitude_Interpolation &
  (XLat,      & ! <-- Latitude grid [deg]
   PLat,      & ! <-- Latitude of point [deg]
   ILat,      & ! --> Interpolation subgrid
   WLat,      & ! --> Weights of subgrid points
   WLat1,     & ! ~~> Weight derivatives
   WLat2)       ! ~~> Weight 2nd derivatives
!
! Latitude interpolation: calculation of interpolation
! weights of grid points and, optionally, the 1st and 2nd
! derivatives of weighting function.
!----------------------------------------------------------
! Method:
!   For N-point scheme, the central cell containing
!   interpolation point is found. Coefficients of
!   polynomial of N-1 power are found so as for
!   derivatives up to (N-2)th to be continuous.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 13 Nov 1998 | Original code.
!   2.0   ! 20 Feb 1999 ! Faster index search.
!----------------------------------------------------------
! Modules used:
!
Use Interpolation, only: &
! Imported Routines:
    Seek_Index
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Workpr), Intent(In) :: &
   XLat(1:)   ! Latitude grid [deg]
!
Real(Double), Intent(In) :: &
   PLat       ! Latitude of point [deg]
!
! Output arguments:
!
Integer, Intent(Out) :: &
   ILat(1:)   ! Indices of grid points used
              ! Size(ILat) must be even
!
Real(Double), Intent(Out) :: &
   WLat(1:)   ! Weights of grid points used
              ! Size(WLat) must be = Size(ILat)
!
Real(Double), Intent(Out), Optional :: &
   WLat1(1:)  ! Weight derivatives
              ! Size(WLat1) must be = Size(ILat)
!
Real(Double), Intent(Out), Optional :: &
   WLat2(1:)  ! Weight 2nd derivatives
              ! Size(WLat2) must be = Size(ILat)
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: NLat  ! Latitude grid dimension
Integer      :: DLat  ! Latitude grid direction
Integer      :: NI    ! Interpolation scheme dimension
Integer      :: i     ! Grid index
Integer      :: IC    ! Lower index of central cell of subgrid
Integer      :: B     ! Border indicator
!
! Local Arrays:
!
Real(Double) :: &
   SLat(Size(ILat))   ! Interpolation subgrid
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF GRID SIZE AND DIRECTION
!----------------------------------------------------------


!--- 1.1. Main grid

NLat = Size(XLat)
DLat = Nint(Sign(1.0_Workpr, XLat(NLat)-XLat(1)))


!--- 1.2. Interpolation subgrid

NI   = Size(ILat)


!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATION CELL
!----------------------------------------------------------


!--- 2.1. Calculation of central cell of interpolation subgrid

IC       = NI/2

!ILat(IC) = Sum(MaxLoc(DLat*XLat(:), &
!                      DLat*XLat(:) <= DLat*PLat))
If (DLat*PLat >= DLat*XLat(NLat)) then
   ILat(IC) = NLat - 1
Else If (DLat*PLat <= DLat*XLat(1)) then
   ILat(IC) = 1
Else
   ILat(IC) = Seek_Index (XLat, Real(PLat, Workpr))
End If
ILat(IC) = Min(Max(1, ILat(IC)), NLat-1)
SLat(IC) = XLat(ILat(IC))


!--- 2.2. Calculation of interpolation subgrid

Do i=1,NI
   ILat(i) = Min(Max(1, ILat(IC)+i-IC), NLat)
   SLat(i) = XLat(ILat(i))
End Do


!--- 2.3. Determination of border position

If (ILat(IC) == 1) then
   B = -1
Else If (ILat(IC) == NLat-1) then
   B = 1
Else
   B = 0
End If


!----------------------------------------------------------
! 3. CALCULATION OF INTERPOLATION WEIGHTS
!----------------------------------------------------------


Select Case (NI)


!--- 3.1. 2-point linear interpolation

Case(2)

   WLat(1) = (SLat(2) - PLat)/(SLat(2)-SLat(1))
   WLat(2) = (PLat - SLat(1))/(SLat(2)-SLat(1))

   If (Present(WLat1)) then
      WLat1(:) = 0
   End If

   If (Present(WLat2)) then
      WLat2(:) = 0
   End If


Case(4)

   Call Weight4 &
     (PLat,    & ! <-- Interpolation point
      SLat,    & ! <-- Interpolation subgrid
      B,       & ! <-- Border indicator
      WLat,    & ! --> Interpolation weights
      WLat1)     ! ~~> Interpolation weight derivatives

   If (Present(WLat2)) then
      WLat2(:) = 0
   End If


!--- 3.x. Default

Case Default

   WLat(:) = 0

   If (Present(WLat1)) then
      WLat1(:) = 0
   End If

   If (Present(WLat2)) then
      WLat2(:) = 0
   End If


End Select


End Subroutine Latitude_Interpolation



!==========================================================
Subroutine Weight4 &
  (XP,        & ! <-- Interpolation point
   X,         & ! <-- Interpolation subgrid
   B,         & ! <-- Border indicator
   W,         & ! --> Interpolation weights
   DW)          ! ~~> Interpolation weight derivatives
!
! Calculation of weights for 4-point interpolation.
!----------------------------------------------------------
! Method:
!   Cubic interpolation between 2nd and 3rd points using
!   2 derivatives from finite differences and 2 values
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov, A. K. Steiner.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 09 Mar 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   XP       ! Interpolation point
!
Real(Double), Intent(In) :: &
   X(:)     ! Interpolation subgrid
!
Integer, Intent(In) :: &
   B        ! Border indicator
            !   -1 - left border
            !    0 - middle
            !    1 - right border
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   W(:)     ! Interpolation weights
!
Real(Double), Optional, Intent(Out) :: &
   DW(:)    ! Interpolation weights
!----------------------------------------------------------



Select Case(B)

   Case (0)
      W(1) = 0.5_Double*(XP-X(2))*(XP-X(3))**2/((X(1)-X(2))*(X(2)-X(3))**2)

      W(2) = ((XP-X(3))/(X(2)-X(3))**3)*     &
             ((XP-X(3))*(3*X(2)-X(3)-2*XP) + &
             0.5_Double*(XP-X(2))*((XP-X(3))*(X(3)-2*X(2)+X(1))/(X(1)-X(2)) + &
                XP-X(2)))

      W(3) = ((XP-X(2))/(X(3)-X(2))**3)*     &
             ((XP-X(2))*(3*X(3)-X(2)-2*XP) + &
             0.5_Double*(XP-X(3))*((XP-X(2))*(X(4)-2*X(3)+X(2))/(X(4)-X(3)) + &
                XP-X(3)))

      W(4) = 0.5_Double*(XP-X(3))*(XP-X(2))**2/((X(4)-X(3))*(X(3)-X(2))**2)

      If (Present(DW)) then

         DW(1) = ((2*X(2) + X(3) - 3*XP)*(X(3) - XP))/   &
                 (2*(X(1) - X(2))*(X(2) - X(3))**2)

         DW(2) = (-X(2)**3 + X(2)**2*(6*X(3) - 4*XP) +   &
                 X(2)*XP*(-4*X(3) + 3*XP) +              &
                 X(1)*(X(2)**2 - 8*X(2)*X(3) + X(3)**2 + &
                 6*X(2)*XP + 6*X(3)*XP - 6*XP**2) +      &
                 X(3)*(X(3)**2 - 4*X(3)*XP + 3*XP**2))/  &
                 (2*(X(1) - X(2))*(X(2) - X(3))**3)

         DW(3) = (X(2)**3 - X(3)**3 + X(2)**2*(X(4) - 4*XP) + &
                 X(3)**2*(X(4) - 4*XP) - 6*X(4)*XP**2 +       &
                 3*X(3)*XP*(2*X(4) + XP) + X(2)*(6*X(3)**2 -  &
                 4*X(3)*(2*X(4) + XP) + &
                 3*XP*(2*X(4) + XP)))/(2*(X(2) - X(3))**3*(X(3) - X(4)))

         DW(4) = -(((X(2) + 2*X(3) - 3*XP)*(X(2) - XP))/ &
                  (2*(X(2) - X(3))**2*(X(3) - X(4))))

      End If

   Case(-1)

      W(1) = 0

      W(2) = 0.5_Double*(XP-X(3))*(X(2)-2*X(3)+XP)/(X(3)-X(2))**2

      W(3) = ((XP-X(2))/(X(3)-X(2))**2)*                &
             (2*X(3)-X(2)-XP + 0.5_Double*(2*X(3)-X(4)-X(2))*  &
             (XP-X(3))/(X(3)-X(4)))

      W(4) = 0.5_Double*(XP-X(3))*(XP-X(2))/((X(4)-X(3))*(X(3)-X(2)))

      If (Present(DW)) then

         DW(1) = 0

         DW(2) = (X(2) - 3*X(3) + 2*XP)/(2*(X(2) - X(3))**2)

         DW(3) = (X(2)**2 + 2*X(3)**2 - 3*X(3)*X(4) + X(2)* &
           (-X(3) + X(4) - 2*XP) + 2*X(4)*XP)/              &
           (2*(X(2) - X(3))**2*(X(3) - X(4)))

         DW(4) = -((X(2) + X(3) - 2*XP)/ &
                  (2*(X(2) - X(3))*(X(3) - X(4))))

      End If

   Case(1)

      W(1) = 0.5_Double*(XP-X(2))*(XP-X(3))/((X(1)-X(2))*(X(2)-X(3)))

      W(2) = ((XP-X(3))/(X(2)-X(3))**2)*                &
             (2*X(2)-X(3)-XP + 0.5_Double*(2*X(2)-X(3)-X(1))*  &
             (XP-X(2))/(X(2)-X(1)))

      W(3) = 0.5_Double*(XP-X(2))*(X(3)-2*X(2)+XP)/(X(3)-X(2))**2

      W(4) = 0

      If (Present(DW)) then

         DW(1) = -((X(2) + X(3) - 2*XP)/(2*(X(1) - X(2))*(X(2) - X(3))))

         DW(2) = (-2*X(2)**2 + X(2)*X(3) + X(1)*                &
                 (3*X(2) - X(3) - 2*XP) - X(3)*(X(3) - 2*XP))/  &
                 (2*(X(1) - X(2))*(X(2) - X(3))**2)

         DW(3) = (-3*X(2) + X(3) + 2*XP)/(2*(X(2) - X(3))**2)

         DW(4) = 0

      End If

End Select

!DEBUG - NO HORIZONTAL GRADIENTS
!DW(:) = 0


End Subroutine Weight4



End Module ECHAM_grid
