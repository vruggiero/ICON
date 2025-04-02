!
!+ GNSS Radio occultation observation operator: Adjoint version of ray-tracer
!
MODULE ECHAM_Rays_adj
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Adjoint version of ray-tracer for global fields.
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
! V1_9         2010/04/20 Harald Anlauf
!  eliminate some NULLIFYs
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
! Module ECHAM_Rays_adj
!
! Adjoint version of ray-tracer for ECHAM global fields.
!----------------------------------------------------------
! (C) Copyright 2000-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 03 May 2000 | Original code.
!   2.0   | 09 Apr 2002 | Icosahedral grid included.
!   3.0   | 17 Apr 2002 | RAStat.
!   3.1   | 20 Apr 2002 | Allocate_Profile_rays_adj.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double
!
Use ECHAM_fields, only: &
! Imported Type Definitions:
    Matrix,          &
    Profile,         &
! Imported Parameters:
    NG1,             &
    fst_Null,        &
    fst_Allocated,   &
    fst_Initialized, &
! Imported Array Variables:
    gf,              &
! Imported Array Variables:
    PStat
!
Use ECHAM_fields_adj, only: &
! Imported Array Variables:
    PAStat
!
Use mo_wmo_tables, only: &
! Imported Parameters:
    WMO6_LATLON,                   &
    WMO6_GAUSSIAN,                 &
!   WMO6_HARMONIC,                 &
    DWD6_ICOSAHEDRON
!
Use mo_exception, only: &
    finish
!----------------------------------------------------------
Implicit None
Private
Public :: Allocate_Profile_rays_adj, ECHAM_rays_adj_init, &
          ECHAM_Rays_adj_cleanup, ECHAM_Refraction
!----------------------------------------------------------
! Private Type Definitions:
!
Type, Private :: Node
   Type(Node),   Pointer :: Previous  ! Backward pointer
   Integer               :: NS        ! Step number
   Integer,      Pointer :: &
      IDX(:,:,:)           ! Subgrid indices [mu, point, index]
   Real(Double), Pointer :: &
      NP_T(:,:,:),       & ! --> d(NP_mu)/d(T(IGP,k))
      NP_Q(:,:,:),       & ! --> d(NP_mu)/d(Q(IGP,k))
      NP_P(:,:),         & ! --> d(NP_mu)/d(Psur(IGP))
      NGN_T(:,:,:,:),    & ! --> d(NGN_mu(i))/d(T(IGP,k))
      NGN_Q(:,:,:,:),    & ! --> d(NGN_mu(i))/d(Q(IGP,k))
      NGN_P(:,:,:)         ! --> d(NGN_mu(i))/d(Psur(IGP))
   Real(Double)          :: &
      B(6,6),            & ! B_{n} matrix
      CMA(4,6,3)           ! C^{mu}_{n}*A matrices
End Type Node
!
! Public Arrays:
!
Type(Matrix), Pointer, Public  :: &
     ZN_T(:,:,:) => NULL(),       & ! d(ZN(i))/d(T(j))
     ZN_Q(:,:,:) => NULL()          ! d(ZN(i))/d(Q(j))
Type(Profile), Pointer, Public :: &
     ZN_P(:,:,:) => NULL()          ! d(ZN(i))/d(Psur(j))
!
Integer, Allocatable,   Public :: &
     RAStat(:,:,:)                  ! Dynamical field status of ZN_{T,Q,P}
!
! Private Scalars:
!
Logical, Private :: &
   Undefined_Pointers = .TRUE. ! Pointer status indicator.
!----------------------------------------------------------
!
Contains

!==========================================================
Subroutine ECHAM_rays_adj_cleanup
!----------------------------------------------------------
! 2. CLEARING DATA
!----------------------------------------------------------

! moved from ECHAM_rays_adj_init to ECHAM_rays_adj_cleanup

Integer :: ErrCode

If (Undefined_Pointers) then
   Nullify(ZN_T)
   Nullify(ZN_Q)
   Nullify(ZN_P)
   Undefined_Pointers = .FALSE.
Else
   Call Clear_Matrices(ZN_T)
   Call Clear_Matrices(ZN_Q)
   Call Clear_Profiles(ZN_P)
End If

Deallocate(RAStat, Stat = ErrCode)

End Subroutine ECHAM_rays_adj_cleanup
!==========================================================
Subroutine Clear_Profiles(F)
   Type(Profile), Pointer :: &
      F(:,:,:)  ! Field of profiles to clear

   Integer :: I1      ! 1st grid index
   Integer :: I2      ! 2nd grid index
   Integer :: ID      ! Diamond index
   Integer :: ErrCode ! System error status

   If (Associated(F)) then
      Do I1 = LBound(F,1),UBound(F,1)
         Do I2 = LBound(F,2),UBound(F,2)
            Do ID = LBound(F,3),UBound(F,3)
!              if (associated (F(I1, I2, ID)%P)) then
               if (RAStat(I1, I2, ID) /= FST_NULL) then
                    Deallocate(F(I1, I2, ID)%P, Stat = ErrCode)
               end if
            End Do
         End Do
      End Do
   End If

   Deallocate(F, Stat = ErrCode)
   Nullify(F)

End Subroutine Clear_Profiles
!==========================================================
Subroutine Clear_Matrices(F)
   Type(Matrix), Pointer :: &
      F(:,:,:)  ! Field of profiles to clear

   Integer :: I1      ! 1st grid index
   Integer :: I2      ! 2nd grid index
   Integer :: ID      ! Diamond index
   Integer :: ErrCode ! System error status

   If (Associated(F)) then
      Do I1 = LBound(F,1),UBound(F,1)
         Do I2 = LBound(F,2),UBound(F,2)
            Do ID = LBound(F,3),UBound(F,3)
!              if (associated (F(I1, I2, ID)%M)) then
               if (RAStat(I1, I2, ID) /= FST_NULL) then
                    Deallocate(F(I1, I2, ID)%M, Stat = ErrCode)
               end if
            End Do
         End Do
      End Do
   End If

   Deallocate(F, Stat = ErrCode)
   Nullify(F)

End Subroutine Clear_Matrices
!==========================================================
Subroutine ECHAM_rays_adj_init
!
! Initialization of ECHAM_rays_adj module.
!----------------------------------------------------------
! Method:
!   Clearing module global data.
!----------------------------------------------------------
! (C) Copyright 2000-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 26 Apr 2000 | Original code.
!   2.0   | 09 Apr 2002 | Icosahedral grid included.
!   2.1   | 17 Apr 2002 | RAStat.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Local Scalars:
!
Integer   :: N1s, N1e ! 1st grid index limits
Integer   :: N2       ! 2nd grid dimension
Integer   :: ND       ! Number of diamonds
!Integer   :: I1, I2   ! Grid indices
!Integer   :: ID       ! Diamond index
Integer   :: ErrCode  ! System error status
!----------------------------------------------------------
! Global variables used:
!
!  ZN_T(:,:)       ! d(ZN(i))/d(T(j))
!----------------------------------------------------------


!----------------------------------------------------------
! 1. DETERMINATION OF FIELD DIMENSIONS
!----------------------------------------------------------

N1S  = gf% lbg(1)
N1E  = gf% ubg(1)
N2   = gf% ubg(2) - gf% lbg(2) + 1
ND   = gf% ubg(3) - gf% lbg(3) + 1


!----------------------------------------------------------
! 2. CLEARING DATA
!----------------------------------------------------------

! moved from ECHAM_rays_adj_init to ECHAM_rays_adj_cleanup

call ECHAM_rays_adj_cleanup

!----------------------------------------------------------
! 2. MEMORY ALLOCATION FOR SURFACE ARRAYS
!----------------------------------------------------------


!--- 2.1. Memory allocation

Allocate(RAStat(N1S:N1E, N2, ND))
RAStat(:,:,:) = fst_Null

Allocate(ZN_T(N1S:N1E, N2, ND),   &
         ZN_Q(N1S:N1E, N2, ND),   &
         ZN_P(N1S:N1E, N2, ND),   &
         Stat = ErrCode)


!--- 2.2. Nullifying vertical profiles and matrices

!Do I1=N1S,N1E
!   Do I2=1,N2
!      Do ID=1,ND
!         Nullify(ZN_T(I1, I2, ID)%M)
!         Nullify(ZN_Q(I1, I2, ID)%M)
!         Nullify(ZN_P(I1, I2, ID)%P)
!      End Do
!   End Do
!End Do


End Subroutine ECHAM_rays_adj_init


!==========================================================
Subroutine Allocate_Profile_rays_adj &
  (I1,       & ! <-- 1st grid index
   I2,       & ! <-- 2nd grid index
   ID)         ! <-- Diamond index

!
! Allocation of profiles for a grid point.
!----------------------------------------------------------
! (C) Copyright 2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 20 Apr 2002 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   I1            ! 1st grid index
!
Integer, Intent(In) :: &
   I2            ! 2nd grid index
!
Integer, Intent(In) :: &
   ID            ! Diamond index
!----------------------------------------------------------
! Local Scalars:
!
Integer        :: NLev    ! Number of levels
!----------------------------------------------------------


!----------------------------------------------------------
! 1. MEMORY ALLOCATION
!----------------------------------------------------------

NLev = gf% nz

Allocate(ZN_T(I1, I2, ID)%M(6,NLev))
Allocate(ZN_Q(I1, I2, ID)%M(6,NLev))
Allocate(ZN_P(I1, I2, ID)%P(6))

RAStat(I1, I2, ID) = fst_Allocated


End Subroutine Allocate_Profile_rays_adj




!==========================================================
Subroutine Ray_Trace_adj &
  (XT,        & ! <-- Transmitter position
   UT,        & ! <-- Transmitter ray direction
   XR,        & ! <-- Receiver position
   S,         & ! <-- Integration step
   XN,        & ! --> Ray point nearest to receiver
   UN,        & ! --> Ray direction at XN
   ZN_Z0,     & ! --> d(XN,UN)/d(XT,UT)
   Stat)        ! --> Error status
!
! Tracing ray with given initial conditions: Adjoint
! version.
!----------------------------------------------------------
! Method:
!   Numerical intergration of the differential equation
!   of rays in cartesian coordinates with given intial
!   conditions, stopping at the point nearest to the
!   receiver.
!----------------------------------------------------------
! (C) Copyright 2000-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   2.0   | 22 Feb 1999 | Basic non-adjoint version.
!   1.0   | 03 May 2000 | Adjoint version.
!   2.0   | 09 Apr 2002 | Icosahedral grid included.
!   2.1   | 17 Apr 2002 | RAStat.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Routines:
    Vector_Normed, &
    Vector_Norm,   &
! Imported Operators:
    Operator(*),   &
    Operator(+),   &
    Operator(-),   &
    Operator(/),   &
    Operator(.x.), &
    Assignment(=)
!
Use Matrix, only: &
! Imported Routines:
    Tensor_Product
!
Use Earth, only: &
! Imported Parameters:
    R_Earth,   &
    H_atm
!
Use ECHAM_fields_adj, only: &
! Imported Routines:
    ECHAM_NGradN_adj
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   XT         ! Transmitter position
              ! Transmitter must be outside atmosphere:
              ! |XT| > R_Earth + H_atm
!
Type(Cartesian), Intent(In) :: &
   UT         ! Transmitter ray direction
!
Type(Cartesian), Intent(In) :: &
   XR         ! Receiver position
              ! Receiver must be outside atmosphere:
              ! |XR| > R_Earth + H_atm
!
Real(Double), Intent(In) :: &
   S          ! Ray integration step
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   XN         ! Ray point nearest to receiver
!
Type(Cartesian), Intent(Out) :: &
   UN         ! Ray direction at XN
!
Real(Double), Intent(Out)    :: &
   ZN_Z0(6,6) ! d(XN,UN)/d(XT,UT)
!
Integer, Intent(Out) :: &
   Stat       ! Error status:
              !    0 - no error
              !    1 - ray collided with Earth
              !    2 - transmitter inside atmosphere
              !    3 - receiver inside atmosphere
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   R_atm = R_Earth + H_atm    ! Radius of atmosphere
!
Real(Double), Parameter :: &
   E3(3,3) = Reshape((/1, 0, 0,  &
                       0, 1, 0,  &
                       0, 0, 1 /),  &
                     (/3,3/)),  & ! Unit 3x3 matrix
   E6(6,6) = Reshape((/1, 0, 0, 0, 0, 0,  &
                       0, 1, 0, 0, 0, 0,  &
                       0, 0, 1, 0, 0, 0,  &
                       0, 0, 0, 1, 0, 0,  &
                       0, 0, 0, 0, 1, 0,  &
                       0, 0, 0, 0, 0, 1 /),  &
                     (/6,6/))     ! Unit 6x6 matrix
!
! Local Scalars:
!
Type(Cartesian) :: X      ! Current ray point coordinate
Type(Cartesian) :: U      ! Current ray direction
Type(Cartesian) :: U1     ! Unnormed vector U
Real(Double)    :: p      ! Starting ray impact parameter
Real(Double)    :: Sa     ! Distance from transmitter to atmosphere
Real(Double)    :: NP     ! Interpolated N
Real(Double)    :: Sr     ! Distance to reciever
Integer         :: NLev   ! Number of model leves
Type(Node), Pointer :: &
           SD, SDnew      ! Derivatives at current step
Integer         :: Mu     ! Substep index
Integer         :: NGP    ! Subgrid dimension
Integer         :: IGP    ! Subgrid index
Integer         :: J1, J2 ! Grid indices
Integer         :: JD     ! Diamond index
!
! Local Arrays:
!
Type(Cartesian) :: &
    DX(4),    & ! Runge-Kutta intermediate dX/dt
    DU(4)       ! Runge-Kutta intermediate dU/dt
Real(Double)    :: &
    NHN(4,3,3)  ! Interpolated Grad x (1+N)Grad(N)
!
! --- Matrices for adjoint calculations
!
Real(Double) :: &
   R(6,6)                          ! Renorming matrix
Real(Double) :: &
   Sa_X(3),                      & ! d(Sa)/d(XT)
   Sa_U(3),                      & ! d(Sa)/d(U)
   Sr_X(3),                      & ! d(Sr)/d(X)
   Sr_U(3)                         ! d(Sr)/d(U)
Real(Double) :: &
   BM(4,6,6),                    & ! B^{mu} matrices
   B21(6,6), B32(6,6), B43(6,6), & ! Combinations of 2 BM
   B321(6,6), B432(6,6),         & ! Combinations of 3 BM
   B4321(6,6)                      ! Combinations of 4 BM
Real(Double) :: &
   CM(4,6,6),                    & ! C^{mu} matrices
   A(6,3)                          ! A matrix
Real(Double) :: &
   BN(6,6),                      & ! Product of B_{n} matrices
   CB(4,6,3)                       ! Product of C^{mu}_{n} and B_{n}
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------

!--- 0.1. Parameter check

If (Vector_Norm(XT) < R_atm) then
   Stat = 2
   Return
End If

If (Vector_Norm(XR) < R_atm) then
   Stat = 3
   Return
End If

Stat = 0


!--- 0.2 Determination of model level number

NLev = gf% nz


!--- 0.3. Calculation of A matrix

A(1:3,1:3) = 0
A(4:6,1:3) = E3(:,:)


!--- 0.4. Subgrid dimension

Select Case (gf% Grid_Type)
   Case (WMO6_GAUSSIAN, WMO6_LATLON)
      NGP = NG1**2
   Case (DWD6_ICOSAHEDRON)
      NGP = 3
   case default
      call finish('Ray_Trace_adj','invalid data type')
End Select


!----------------------------------------------------------
! 1. FIRST STEP FROM TRANSMITTER TO ATMOSPHERE
!----------------------------------------------------------


!--- 1.1. First step

U  = Vector_Normed(UT)
p  = Sqrt((XT*XT) - (XT*U)**2)
Sa = -(XT*U) - Sqrt(Dim(R_atm**2, p**2))
X  = XT + Sa*U


!--- 1.2. Initialization step list

Allocate(SD)
Nullify(SD%Previous)
SD%NS = 1
Nullify  &
  (SD%IDX,      & ! Subgrid indices [IGP, index]
   SD%NP_T,     & ! d(NP)/d(T(IGP,k))
   SD%NP_Q,     & ! d(NP)/d(Q(IGP,k))
   SD%NP_P,     & ! d(NP)/d(Psur(IGP))
   SD%NGN_T,    & ! d(NGN(i))/d(T(IGP,k))
   SD%NGN_Q,    & ! d(NGN(i))/d(Q(IGP,k))
   SD%NGN_P)      ! d(NGN(i))/d(Psur(IGP))


!--- 1.3. Calculation of B_{1} matrix

R(:,:) = 0
R(1:3,1:3) = E3(:,:)
R(4:6,4:6) = &
   (E3(:,:) - Tensor_Product(UT%X(:),UT%X(:))/Sqrt(UT*UT))/ &
   Sqrt(UT*UT)

Sa_X(:) = -U%X(:) - &
          (-XT%X(:) + (XT*U)*U%X(:))/Sqrt(Dim(R_atm**2, p**2))
Sa_U(:) = -XT%X(:) - &
          ((XT*U)*XT%X(:))/Sqrt(Dim(R_atm**2, p**2))

SD%B(1:3,1:3) = E3(:,:) + Tensor_Product(U%X(:),Sa_X(:))
SD%B(1:3,4:6) = Sa*E3(:,:) + Tensor_Product(U%X(:),Sa_U(:))
SD%B(4:6,4:6) = E3(:,:)
SD%B(4:6,1:3) = 0
SD%CMA(:,:,:) = 0

SD%B(:,:) = MatMul(SD%B(:,:),R(:,:))


!----------------------------------------------------------
! 2. RAY INTEGRATION
!----------------------------------------------------------


Ray_Integrate: Do


   !--- 2.1. Checking for ray outgoing from atmosphere

   If (((X*U) > 0) .and. (Vector_Norm(X) > R_atm)) then
      Exit Ray_Integrate
   End If


   !--- 2.2. Adding node to step list

   Allocate(SDnew)
   SDnew%Previous => SD
   SDnew%NS       =  SD%NS + 1
   SD             => SDnew
   Allocate  &
     (SD%IDX(4,NGP,3),         & ! Subgrid indices [IGP, index]
      SD%NP_T(4,NGP,NLev),     & ! d(NP)/d(T(IGP,k))
      SD%NP_Q(4,NGP,NLev),     & ! d(NP)/d(Q(IGP,k))
      SD%NP_P(4,NGP),          & ! d(NP)/d(Psur(IGP))
      SD%NGN_T(4,3,NGP,NLev),  & ! d(NGN(i))/d(T(IGP,k))
      SD%NGN_Q(4,3,NGP,NLev),  & ! d(NGN(i))/d(Q(IGP,k))
      SD%NGN_P(4,3,NGP))         ! d(NGN(i))/d(Psur(IGP))


    !--- 2.3. Making a step of Runge-Kutta ray integration

   DX(1) = U
   Call ECHAM_NGradN_adj &
     (X,                    & ! <-- Cartesian coordinates of point
      DU(1),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(1,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat,                 & ! --> Error status
      SD%IDX(1,:,:),        & ! --> Subgrid indices [IGP, index]
      SD%NP_T(1,:,:),       & ! --> d(NP)/d(T(IGP,k))
      SD%NP_Q(1,:,:),       & ! --> d(NP)/d(Q(IGP,k))
      SD%NP_P(1,:),         & ! --> d(NP)/d(Psur(IGP))
      SD%NGN_T(1,:,:,:),    & ! --> d(NGN(i))/d(T(IGP,k))
      SD%NGN_Q(1,:,:,:),    & ! --> d(NGN(i))/d(Q(IGP,k))
      SD%NGN_P(1,:,:))        ! --> d(NGN(i))/d(Psur(IGP))
   DX(2) = U + DU(1)*(S/2)
   Call ECHAM_NGradN_adj &
     (X + DX(1)*(S/2),      & ! <-- Cartesian coordinates of point
      DU(2),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(2,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat,                 & ! --> Error status
      SD%IDX(2,:,:),        & ! --> Subgrid indices [IGP, index]
      SD%NP_T(2,:,:),       & ! --> d(NP)/d(T(IGP,k))
      SD%NP_Q(2,:,:),       & ! --> d(NP)/d(Q(IGP,k))
      SD%NP_P(2,:),         & ! --> d(NP)/d(Psur(IGP))
      SD%NGN_T(2,:,:,:),    & ! --> d(NGN(i))/d(T(IGP,k))
      SD%NGN_Q(2,:,:,:),    & ! --> d(NGN(i))/d(Q(IGP,k))
      SD%NGN_P(2,:,:))        ! --> d(NGN(i))/d(Psur(IGP))
   DX(3) = U + DU(2)*(S/2)
   Call ECHAM_NGradN_adj &
     (X + DX(2)*(S/2),      & ! <-- Cartesian coordinates of point
      DU(3),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(3,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat,                 & ! --> Error status
      SD%IDX(3,:,:),        & ! --> Subgrid indices [IGP, index]
      SD%NP_T(3,:,:),       & ! --> d(NP)/d(T(IGP,k))
      SD%NP_Q(3,:,:),       & ! --> d(NP)/d(Q(IGP,k))
      SD%NP_P(3,:),         & ! --> d(NP)/d(Psur(IGP))
      SD%NGN_T(3,:,:,:),    & ! --> d(NGN(i))/d(T(IGP,k))
      SD%NGN_Q(3,:,:,:),    & ! --> d(NGN(i))/d(Q(IGP,k))
      SD%NGN_P(3,:,:))        ! --> d(NGN(i))/d(Psur(IGP))
   DX(4) = U + DU(3)*S
   Call ECHAM_NGradN_adj &
     (X + DX(3)*S,          & ! <-- Cartesian coordinates of point
      DU(4),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(4,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat,                 & ! --> Error status
      SD%IDX(4,:,:),        & ! --> Subgrid indices [IGP, index]
      SD%NP_T(4,:,:),       & ! --> d(NP)/d(T(IGP,k))
      SD%NP_Q(4,:,:),       & ! --> d(NP)/d(Q(IGP,k))
      SD%NP_P(4,:),         & ! --> d(NP)/d(Psur(IGP))
      SD%NGN_T(4,:,:,:),    & ! --> d(NGN(i))/d(T(IGP,k))
      SD%NGN_Q(4,:,:,:),    & ! --> d(NGN(i))/d(Q(IGP,k))
      SD%NGN_P(4,:,:))        ! --> d(NGN(i))/d(Psur(IGP))

   X = X + (DX(1) + 2*DX(2) + 2*DX(3) + DX(4))*(S/6)
   U = U + (DU(1) + 2*DU(2) + 2*DU(3) + DU(4))*(S/6)


   !--- 2.4. Renorming U

   U1 = U
   U  = Vector_Normed(U)*(1+NP)


   !--- 2.5. Checking for ray collision with Earth

   If (Stat /= 0) then
      UN = U
      XN = X
      Clear: Do
         If(.not. Associated(SD%Previous)) then
            Deallocate(SD)
            Exit Clear
         End If
         SDnew => SD%Previous
         Deallocate  &
           (SD%IDX,     & ! Subgrid indices [IGP, index]
            SD%NP_T,    & ! d(NP)/d(T(KLon,KLat,k))
            SD%NP_Q,    & ! d(NP)/d(Q(KLon,KLat,k))
            SD%NP_P,    & ! d(NP)/d(Psur(KLon,KLat))
            SD%NGN_T,   & ! d(NGN(i))/d(T(KLon,KLat,k))
            SD%NGN_Q,   & ! d(NGN(i))/d(Q(KLon,KLat,k))
            SD%NGN_P)     ! d(NGN(i))/d(Psur(KLon,KLat))
         Deallocate(SD)
         SD => SDnew
      End Do Clear
      Return
   End If


   !--- 2.6. Calculation of BM matrices

   BM(:,:,:) = 0
   Do Mu=1,4
      BM(Mu,1:3,4:6) = E3(:,:)
      BM(Mu,4:6,1:3) = NHN(Mu,:,:)
   End Do

   B21(:,:)   = MatMul(BM(2,:,:), BM(1,:,:))
   B32(:,:)   = MatMul(BM(3,:,:), BM(2,:,:))
   B43(:,:)   = MatMul(BM(4,:,:), BM(3,:,:))

   B321(:,:)  = MatMul(BM(3,:,:), B21(:,:))
   B432(:,:)  = MatMul(BM(4,:,:), B32(:,:))

   B4321(:,:) = MatMul(BM(4,:,:), B321(:,:))


   !--- 2.7. Calculation of B matrix

   SD%B(:,:) = &
      E6(:,:) + &
      (S/6)*(BM(1,:,:) + 2*BM(2,:,:) + 2*BM(3,:,:) + BM(4,:,:)) + &
      (S**2/6)*(B21(:,:) + B32(:,:) + B43(:,:)) + &
      (S**3/12)*(B321(:,:) + B432(:,:)) + &
      (S**4/24)*B4321(:,:)


   !--- 2.8. Calculation of C^{mu} matrices

   CM(1,:,:) = &
      (S/6)*E6(:,:)      + &
      (S**2/6)*BM(2,:,:) + &
      (S**3/12)*B32(:,:) + &
      (S**4/24)*B432(:,:)
   CM(2,:,:) = &
      (S/3)*E6(:,:)      + &
      (S**2/6)*BM(3,:,:) + &
      (S**3/12)*B43(:,:)
   CM(3,:,:) = &
      (S/3)*E6(:,:)      + &
      (S**2/6)*BM(4,:,:)
   CM(4,:,:) = &
      (S/6)*E6(:,:)


   !--- 2.9. Multiplication of C^{mu} with A

   Do Mu=1,4
      SD%CMA(Mu,:,:) = MatMul(CM(Mu,:,:),A(:,:))
   End Do


   !--- 2.10. Calculation of renorming matrix

   R(1:3,1:3) = E3(:,:)
   R(1:3,4:6) = 0
   R(4:6,1:3) = &
      Tensor_Product(U1%X(:),DU(4)%X(:))/ &
      ((1 + NP)*Sqrt(U1*U1))
   R(4:6,4:6) = &
      (1 + NP)*(E3(:,:) - Tensor_Product(U1%X(:),U1%X(:))/Sqrt(U1*U1))/ &
      Sqrt(U1*U1)


   !--- 2.11. Multiplication of B and C^{mu} with R

   SD%B(:,:) = MatMul(R(:,:),SD%B(:,:))
   Do Mu=1,4
      SD%CMA(Mu,:,:) = MatMul(R(:,:),SD%CMA(Mu,:,:))
   End Do


End Do Ray_Integrate


!----------------------------------------------------------
! 3. LAST STEP FROM ATMOSPHERE TO RECEIVER
!----------------------------------------------------------


!--- 3.1. Last step

UN = Vector_Normed(U)
XN = X

Sr = (XR - X)*UN
XN = XN + Sr*UN


!--- 3.2. Calculation of B_{N} matrix

R(:,:) = 0
R(1:3,1:3) = E3(:,:)
R(4:6,4:6) = &
   (E3(:,:) - Tensor_Product(U%X(:),U%X(:))/Sqrt(U*U))/ &
   Sqrt(U*U)

Sr_X(:) = -UN%X(:)
Sr_U(:) = XR%X(:) - X%X(:)

BN(1:3,1:3) = E3(:,:) + Tensor_Product(U%X(:),Sr_X(:))
BN(1:3,4:6) = Sr*E3(:,:) + Tensor_Product(U%X(:),Sr_U(:))
BN(4:6,4:6) = E3(:,:)
BN(4:6,1:3) = 0

BN(:,:) = MatMul(BN(:,:),R(:,:))


!----------------------------------------------------------
! 4. BACKWARD ADJOINT INTEGRATION
!----------------------------------------------------------


Adj_Integrate: Do


   !--- 4.1. Calculation of C matrices

   Do Mu=1,4
      CB(Mu,:,:) = MatMul(BN(:,:), SD%CMA(Mu,:,:))
   End Do


   !--- 4.2. Calculation of B matrix

   BN(:,:) = MatMul(BN(:,:), SD%B(:,:))


   !--- 4.3. Exit if the first step reached

   If(.not. Associated(SD%Previous)) then
      Deallocate(SD)
      Exit Adj_Integrate
   End If


   !--- 4.4. Calculation of d(ZN)/d(T(i,j,k),Q(i,j,k),Psur(i,j))

   Do Mu=1,4
      Do IGP=1,NGP
         J1 = SD%IDX(Mu,IGP,1)
         J2 = SD%IDX(Mu,IGP,2)
         JD = SD%IDX(Mu,IGP,3)
         If (RAStat(J1, J2, JD) == fst_Null) then
            Allocate(ZN_T(J1,J2,JD)%M(6,NLev))
            Allocate(ZN_Q(J1,J2,JD)%M(6,NLev))
            Allocate(ZN_P(J1,J2,JD)%P(6))
            RAStat(J1, J2, JD)    = fst_Allocated
         End If
         If (RAStat(J1, J2, JD) == fst_Allocated) then
            ZN_T(J1,J2,JD)%M(:,:) = 0
            ZN_Q(J1,J2,JD)%M(:,:) = 0
            ZN_P(J1,J2,JD)%P(:)   = 0
            RAStat(J1, J2, JD)    = fst_Initialized
         End If
         ZN_T(J1,J2,JD)%M(:,:) = &
            ZN_T(J1,J2,JD)%M(:,:) + &
            MatMul(CB(Mu,:,:),SD%NGN_T(Mu,:,IGP,:))
         ZN_Q(J1,J2,JD)%M(:,:) = &
            ZN_Q(J1,J2,JD)%M(:,:) + &
            MatMul(CB(Mu,:,:),SD%NGN_Q(Mu,:,IGP,:))
         ZN_P(J1,J2,JD)%P(:)   = &
            ZN_P(J1,J2,JD)%P(:) + &
            MatMul(CB(Mu,:,:),SD%NGN_P(Mu,:,IGP))
      End Do
   End Do


   !--- 4.5. Deallocation of current step data and
   !---      taking next.

   SDnew => SD%Previous

   Deallocate  &
     (SD%IDX,     & ! Subgrid indices [IGP, index]
      SD%NP_T,    & ! d(NP)/d(T(KLon,KLat,k))
      SD%NP_Q,    & ! d(NP)/d(Q(KLon,KLat,k))
      SD%NP_P,    & ! d(NP)/d(Psur(KLon,KLat))
      SD%NGN_T,   & ! d(NGN(i))/d(T(KLon,KLat,k))
      SD%NGN_Q,   & ! d(NGN(i))/d(Q(KLon,KLat,k))
      SD%NGN_P)     ! d(NGN(i))/d(Psur(KLon,KLat))
   Deallocate(SD)

   SD => SDnew


End Do Adj_Integrate


ZN_Z0(:,:) = BN(:,:)


End Subroutine Ray_Trace_adj



!==========================================================
Subroutine Ray_Trace_ZNZ0 &
  (XT,        & ! <-- Transmitter position
   UT,        & ! <-- Transmitter ray direction
   XR,        & ! <-- Receiver position
   S,         & ! <-- Integration step
   XN,        & ! --> Ray point nearest to receiver
   UN,        & ! --> Ray direction at XN
   ZN_Z0,     & ! --> d(XN,UN)/d(XT,UT)
   Stat)        ! --> Error status
!
! Tracing ray with given initial conditions and calculation
! of derivative d(ZN)/d(Z0).
!----------------------------------------------------------
! Method:
!   Numerical intergration of the differential equation
!   of rays in cartesian coordinates with given intial
!   conditions, stopping at the point nearest to the
!   receiver.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   2.0   | 22 Feb 1999 | Basic non-adjoint version.
!   1.0   | 03 May 2000 | Reduced adjoint version.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Routines:
    Vector_Normed, &
    Vector_Norm,   &
! Imported Operators:
    Operator(*),   &
    Operator(+),   &
    Operator(-),   &
    Operator(/),   &
    Operator(.x.), &
    Assignment(=)
!
Use Matrix, only: &
! Imported Routines:
    Tensor_Product
!
Use Earth, only: &
! Imported Parameters:
    R_Earth,   &
    H_atm
!
Use ECHAM_fields, only: &
! Imported Routines:
    ECHAM_NGHN
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   XT         ! Transmitter position
              ! Transmitter must be outside atmosphere:
              ! |XT| > R_Earth + H_atm
!
Type(Cartesian), Intent(In) :: &
   UT         ! Transmitter ray direction
!
Type(Cartesian), Intent(In) :: &
   XR         ! Receiver position
              ! Receiver must be outside atmosphere:
              ! |XR| > R_Earth + H_atm
!
Real(Double), Intent(In) :: &
   S          ! Ray integration step
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   XN         ! Ray point nearest to receiver
!
Type(Cartesian), Intent(Out) :: &
   UN         ! Ray direction at XN
!
Real(Double), Intent(Out)    :: &
   ZN_Z0(6,6) ! d(XN,UN)/d(XT,UT)
!
Integer, Intent(Out) :: &
   Stat       ! Error status:
              !    0 - no error
              !    1 - ray collided with Earth
              !    2 - transmitter inside atmosphere
              !    3 - receiver inside atmosphere
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   R_atm = R_Earth + H_atm    ! Radius of atmosphere
!
Real(Double), Parameter :: &
   E3(3,3) = Reshape((/1, 0, 0,  &
                       0, 1, 0,  &
                       0, 0, 1 /),  &
                     (/3,3/)),  & ! Unit 3x3 matrix
   E6(6,6) = Reshape((/1, 0, 0, 0, 0, 0,  &
                       0, 1, 0, 0, 0, 0,  &
                       0, 0, 1, 0, 0, 0,  &
                       0, 0, 0, 1, 0, 0,  &
                       0, 0, 0, 0, 1, 0,  &
                       0, 0, 0, 0, 0, 1 /),  &
                     (/6,6/))     ! Unit 6x6 matrix
!
! Local Scalars:
!
Type(Cartesian) :: X      ! Current ray point coordinate
Type(Cartesian) :: U      ! Current ray direction
Type(Cartesian) :: U1     ! Unnormed vector U
Real(Double)    :: p      ! Starting ray impact parameter
Real(Double)    :: Sa     ! Distance from transmitter to atmosphere
Real(Double)    :: NP     ! Interpolated N
Real(Double)    :: Sr     ! Distance to reciever
Integer         :: Mu     ! Integration substep index
!
! Local Arrays:
!
Type(Cartesian) :: &
    DX(4),    & ! Runge-Kutta intermediate dX/dt
    DU(4)       ! Runge-Kutta intermediate dU/dt
Real(Double)    :: &
    NHN(4,3,3)  ! Interpolated Grad x (1+N)Grad(N)
!
! --- Matrices for adjoint calculations
!
Real(Double) :: &
   R(6,6)                          ! Renorming matrix
Real(Double) :: &
   Sa_X(3),                      & ! d(Sa)/d(XT)
   Sa_U(3),                      & ! d(Sa)/d(U)
   Sr_X(3),                      & ! d(Sr)/d(X)
   Sr_U(3)                         ! d(Sr)/d(U)
Real(Double) :: &
   BM(4,6,6),                    & ! B^{mu} matrices
   B21(6,6), B32(6,6), B43(6,6), & ! Combinations of 2 BM
   B321(6,6), B432(6,6),         & ! Combinations of 3 BM
   B4321(6,6)                      ! Combinations of 4 BM
Real(Double) :: &
   B(6,6),                       & ! B_{n} matrix
   BN(6,6)                         ! Product of B_{n} matrices
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------

!--- 0.1. Parameter check

If (Vector_Norm(XT) < R_atm) then
   Stat = 2
   Return
End If

If (Vector_Norm(XR) < R_atm) then
   Stat = 3
   Return
End If

Stat = 0


!----------------------------------------------------------
! 1. FIRST STEP FROM TRANSMITTER TO ATMOSPHERE
!----------------------------------------------------------


!--- 1.1. First step

U  = Vector_Normed(UT)
p  = Sqrt((XT*XT) - (XT*U)**2)
Sa = -(XT*U) - Sqrt(Dim(R_atm**2, p**2))
X  = XT + Sa*U


!--- 1.2. Calculation of B_{1} matrix

R(:,:) = 0
R(1:3,1:3) = E3(:,:)
R(4:6,4:6) = &
   (E3(:,:) - Tensor_Product(UT%X(:),UT%X(:))/Sqrt(UT*UT))/ &
   Sqrt(UT*UT)

Sa_X(:) = -U%X(:) - &
          (-XT%X(:) + (XT*U)*U%X(:))/Sqrt(Dim(R_atm**2, p**2))
Sa_U(:) = -XT%X(:) - &
          ((XT*U)*XT%X(:))/Sqrt(Dim(R_atm**2, p**2))

B(1:3,1:3) = E3(:,:) + Tensor_Product(U%X(:),Sa_X(:))
B(1:3,4:6) = Sa*E3(:,:) + Tensor_Product(U%X(:),Sa_U(:))
B(4:6,4:6) = E3(:,:)
B(4:6,1:3) = 0

B(:,:) = MatMul(B(:,:),R(:,:))


!--- 1.3. Initialization of BN matrix

BN(:,:) = B(:,:)


!----------------------------------------------------------
! 2. RAY INTEGRATION
!----------------------------------------------------------


Ray_Integrate: Do


   !--- 2.1. Checking for ray outgoing from atmosphere

   If (((X*U) > 0) .and. (Vector_Norm(X) > R_atm)) then
      Exit Ray_Integrate
   End If


    !--- 2.3. Making a step of Runge-Kutta ray integration

   DX(1) = U
   Call ECHAM_NGHN &
     (X,                    & ! <-- Cartesian coordinates of point
      DU(1),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(1,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat)                   ! --> Error status
   DX(2) = U + DU(1)*(S/2)
   Call ECHAM_NGHN &
     (X + DX(1)*(S/2),      & ! <-- Cartesian coordinates of point
      DU(2),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(2,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat)                   ! --> Error status
   DX(3) = U + DU(2)*(S/2)
   Call ECHAM_NGHN &
     (X + DX(2)*(S/2),      & ! <-- Cartesian coordinates of point
      DU(3),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(3,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat)                   ! --> Error status
   DX(4) = U + DU(3)*S
   Call ECHAM_NGHN &
     (X + DX(3)*S,          & ! <-- Cartesian coordinates of point
      DU(4),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(4,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat)                   ! --> Error status

   X = X + (DX(1) + 2*DX(2) + 2*DX(3) + DX(4))*(S/6)
   U = U + (DU(1) + 2*DU(2) + 2*DU(3) + DU(4))*(S/6)


   !--- 2.4. Renorming U

   U1 = U
   U  = Vector_Normed(U)*(1+NP)


   !--- 2.5. Checking for ray collision with Earth

   If (Stat /= 0) then
      UN = U
      XN = X
      Return
   End If


   !--- 2.6. Calculation of BM matrices

   BM(:,:,:) = 0
   Do Mu=1,4
      BM(Mu,1:3,4:6) = E3(:,:)
      BM(Mu,4:6,1:3) = NHN(Mu,:,:)
   End Do

   B21(:,:)   = MatMul(BM(2,:,:), BM(1,:,:))
   B32(:,:)   = MatMul(BM(3,:,:), BM(2,:,:))
   B43(:,:)   = MatMul(BM(4,:,:), BM(3,:,:))

   B321(:,:)  = MatMul(BM(3,:,:), B21(:,:))
   B432(:,:)  = MatMul(BM(4,:,:), B32(:,:))

   B4321(:,:) = MatMul(BM(4,:,:), B321(:,:))


   !--- 2.7. Calculation of B matrix

   B(:,:) = &
      E6(:,:) + &
      (S/6)*(BM(1,:,:) + 2*BM(2,:,:) + 2*BM(3,:,:) + BM(4,:,:)) + &
      (S**2/6)*(B21(:,:) + B32(:,:) + B43(:,:)) + &
      (S**3/12)*(B321(:,:) + B432(:,:)) + &
      (S**4/24)*B4321(:,:)


   !--- 2.8. Calculation of renorming matrix

   R(1:3,1:3) = E3(:,:)
   R(1:3,4:6) = 0
   R(4:6,1:3) = &
      Tensor_Product(U1%X(:),DU(4)%X(:))/ &
      ((1 + NP)*Sqrt(U1*U1))
   R(4:6,4:6) = &
      (1 + NP)*(E3(:,:) - Tensor_Product(U1%X(:),U1%X(:))/Sqrt(U1*U1))/ &
      Sqrt(U1*U1)


   !--- 2.9. Multiplication of B and C^{mu} with R

   B(:,:) = MatMul(R(:,:),B(:,:))


   !--- 2.10. Calculation of BN

   BN(:,:) = MatMul(B(:,:),BN(:,:))


End Do Ray_Integrate


!----------------------------------------------------------
! 3. LAST STEP FROM ATMOSPHERE TO RECEIVER
!----------------------------------------------------------


!--- 3.1. Last step

UN = Vector_Normed(U)
XN = X

Sr = (XR - X)*UN
XN = XN + Sr*UN


!--- 3.2. Calculation of B_{N} matrix

R(:,:) = 0
R(1:3,1:3) = E3(:,:)
R(4:6,4:6) = &
   (E3(:,:) - Tensor_Product(U%X(:),U%X(:))/Sqrt(U*U))/ &
   Sqrt(U*U)

Sr_X(:) = -UN%X(:)
Sr_U(:) = XR%X(:) - X%X(:)

B(1:3,1:3) = E3(:,:) + Tensor_Product(U%X(:),Sr_X(:))
B(1:3,4:6) = Sr*E3(:,:) + Tensor_Product(U%X(:),Sr_U(:))
B(4:6,4:6) = E3(:,:)
B(4:6,1:3) = 0

B(:,:) = MatMul(B(:,:),R(:,:))


!--- 3.3. Calculation of BN

BN(:,:) = MatMul(B(:,:),BN(:,:))


!--- 3.4. Setting ZN_Z0

ZN_Z0(:,:) = BN(:,:)


End Subroutine Ray_Trace_ZNZ0



!==========================================================
Subroutine Find_Ray_P &
  (XT,     & ! <-- Transmitter position
   VT,     & ! <-- Transmitter velocity
   XR,     & ! <-- Receiver position
   VR,     & ! <-- Receiver velocity
   S,      & ! <-- Integration step
   XLC,    & ! <-- Local curvature center
   P,      & ! <-- Impact parameter [km]
   UT,     & ! --> Transmitter ray direction
   XN,     & ! --> Ray point nearest to receiver
   UN,     & ! --> Ray direction at XN
   Stat)     ! --> Error status
!
! Finding ray with prescribed impact parameter.
!----------------------------------------------------------
! Method:
!   Iterative Newton solution.
!----------------------------------------------------------
! (C) Copyright 2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 04 May 2000 | Original code.
!   1.1   | 20 Dec 2001 | Accurate combination of
!         |             | Newton method and binary search.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Routines:
    Rotate,        &
    Vector_Normed, &
    Vector_Norm,   &
    Vector_Angle,  &
! Imported Operators:
    Operator(*),   &
    Operator(+),   &
    Operator(-),   &
    Operator(/),   &
    Operator(.x.), &
    Assignment(=)
!
Use Occ_Refraction_adj, only: &
! Imported Routines:
    Ray_Refraction_adj,   &
    Ray_Derivatives
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In)  :: &
   XT        ! Transmitter position
Type(Cartesian), Intent(In)  :: &
   VT        ! Transmitter velocity
Type(Cartesian), Intent(In)  :: &
   XR        ! Receiver position
Type(Cartesian), Intent(In)  :: &
   VR        ! Receiver velocity
Real(Double), Intent(In)     :: &
   S         ! Integration step
Type(Cartesian), Intent(In)  :: &
   XLC       ! Local curvature center
Real(Double), Intent(In)     :: &
   P         ! Impact parameter [km]
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   UT        ! Transmitter ray direction
Type(Cartesian), Intent(Out) :: &
   XN        ! Ray point nearest to receiver
Type(Cartesian), Intent(Out) :: &
   UN        ! Ray direction at XN
Integer, Intent(Out)         :: &
   Stat      ! Error status
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   Pacc = 1d-6         ! Accuracy of impact parameter
Integer, Parameter :: &
   Nit  = 20           ! Maximum number of iterations
!
! Local Scalars:
!
Real(Double) :: Alpha        ! Leveling angle
Real(Double) :: ER           ! Ray refraction angle
Real(Double) :: PR           ! Ray impact parameter
Real(Double) :: DeltaP       ! Impact parameter correction
Integer      :: i            ! Iteration index
Real(Double) :: D            ! Estimate of observation distance
Real(Double) :: AlphaMin     ! Lower estimate of Alpha
Real(Double) :: AlphaMax     ! Upper estimate of Alpha
Integer      :: RStat        ! Ray-tracing status
!
! Local Arrays:
!
! --- Geometric parameters
!
Real(Double) :: Z0(6)        ! Initial ray point (X0,U0)
Real(Double) :: ZN(6)        ! Final ray point (XN,UN)
!
! --- Ray derivatives
!
Real(Double) :: ZN_Z0(6,6)   ! d(XN,UN)/d(XT,UT)
Real(Double) :: E_Z0(6)      ! d(Eps)/d(Z0)
Real(Double) :: E_ZN(6)      ! d(Eps)/d(ZN)
Real(Double) :: P_Z0(6)      ! d(P)/d(Z0)
Real(Double) :: P_ZN(6)      ! d(P)/d(ZN)
Real(Double) :: FZ0_P(6)     ! D(Z0)/D(P)
Real(Double) :: FE_ZN(6)     ! D(E)/D(ZN)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIAL APPROXIMATION
!----------------------------------------------------------

Alpha = ASin(P/Vector_Norm(XT-XLC))

UT = Vector_Normed(Rotate(-(XT-XLC), (XR-XLC).x.(XT-XLC), Alpha))

D = Sqrt((XT-XLC)*(XT-XLC) - P**2)

AlphaMin = Alpha - 20.0_Double/D
AlphaMax = Alpha + 20.0_Double/D


!----------------------------------------------------------
! 2. ITERATIVE SOLUTION
!----------------------------------------------------------


Stat = 1

Newton: Do i=1,Nit

   Alpha = Vector_Angle(UT, -(XT-XLC))

   If ((Alpha < AlphaMin) .or. (Alpha > AlphaMax)) then
      UT = Vector_Normed(Rotate(-(XT-XLC),     &
                         (XR-XLC).x.(XT-XLC),  &
                         (AlphaMin + AlphaMax)/2))
      Alpha = Vector_Angle(UT, -(XT-XLC))
   End If

   Call Ray_Trace_ZNZ0 &
     (XT,        & ! <-- Transmitter position
      UT,        & ! <-- Transmitter ray direction
      XR,        & ! <-- Receiver position
      S,         & ! <-- Integration step
      XN,        & ! --> Ray point nearest to receiver
      UN,        & ! --> Ray direction at XN
      ZN_Z0,     & ! --> d(XN,UN)/d(XT,UT)
      RStat)       ! --> Error status

   If (RStat /= 0) then
      AlphaMin = Max(AlphaMin,Alpha)
      UT = Vector_Normed(Rotate(-(XT-XLC),     &
                         (XR-XLC).x.(XT-XLC),  &
                         (AlphaMin + AlphaMax)/2))
      Cycle Newton
   End If

   Z0(1:3) = XT - XLC
   Z0(4:6) = UT
   ZN(1:3) = XN - XLC
   ZN(4:6) = UN

   Call Ray_Refraction_adj &
     (Z0,        & ! <-- Initial ray point (X0,U0)
      ZN,        & ! <-- Final ray point (XN,UN)
      VT,        & ! <-- Transmitter velocity [km/s]
      VR,        & ! <-- Receiver velocity [km/s]
      ER,        & ! --> Refraction angle [rad]
      PR,        & ! --> Impact parameter [km]
      E_Z0,      & ! --> d(Eps)/d(Z0)
      E_ZN,      & ! --> d(Eps)/d(ZN)
      P_Z0,      & ! --> d(P)/d(Z0)
      P_ZN)        ! --> d(P)/d(ZN)

   DeltaP = P - PR

   If (PR < P) then
      AlphaMin = Max(AlphaMin,Alpha)
   Else
      AlphaMax = Min(AlphaMax,Alpha)
   End If

   If (Abs(DeltaP) < Pacc) then
      Stat = 0
      Exit Newton
   Else
      Stat = 1
   End If

   Call Ray_Derivatives &
     (Z0,        & ! <-- Initial ray point (X0,U0)
      ZN,        & ! <-- Final ray point (XN,UN)
      ZN_Z0,     & ! <-- d(ZN)/d(Z0)
      P_Z0,      & ! <-- d(P)/d(Z0)
      P_ZN,      & ! <-- d(P)/d(ZN)
      E_Z0,      & ! <-- d(E)/d(Z0)
      E_ZN,      & ! <-- d(D)/d(ZN)
      FZ0_P,     & ! --> D(Z0)/D(P)
      FE_ZN)       ! --> D(E)/D(ZN)

   UT = UT + Cartesian(DeltaP*FZ0_P(4:6))


End Do Newton



End Subroutine Find_Ray_P



!==========================================================
Subroutine ECHAM_Refraction &
  (PRO,       & ! <-- RO impact parameter [km]
   ERLC,      & ! <-- Curvature center in Earth frame
   ERLEO,     & ! <-- LEO coordinates in Earth frame
   EVLEO,     & ! <-- LEO velocities in Earth frame
   ERGPS,     & ! <-- GPS coordinates in Earth frame
   EVGPS,     & ! <-- GPS velocities in Earth frame
   EM,        & ! --> Model refraction angle
   FE_ZN,     & ! --> D(E)/D(ZN)
   Stat)        ! --> Error status
!
! The observational operator for the radiooccultation
! measurements and its adjoint.
!----------------------------------------------------------
! Method:
!   Finding ray with given impact parameter and
!   calculation of ray derivatives.
!----------------------------------------------------------
! (C) Copyright 2000-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 05 May 2000 | Original code.
!   2.0   | 09 Apr 2002 | Icosahedral grid included.
!   2.1   | 17 Apr 2002 | RAStat.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Operators:
    Operator(*),   &
    Operator(+),   &
    Operator(-),   &
    Operator(/),   &
    Operator(.x.), &
    Assignment(=)
!
Use Occ_Refraction_adj, only: &
! Imported Routines:
    Ray_Refraction_adj,   &
    Ray_Derivatives
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   PRO        ! Impact parameters [km]
!
Type(Cartesian), Intent(In) :: &
   ERLC       ! Curvature center in Earth frame
!
Type(Cartesian), Intent(In) :: &
   ERLEO      ! LEO coordinates in Earth fram
!
Type(Cartesian), Intent(In) :: &
   EVLEO      ! LEO velocities in Earth frame
!
Type(Cartesian), Intent(In) :: &
   ERGPS      ! GPS coordinates in Earth frame
!
Type(Cartesian), Intent(In) :: &
   EVGPS      ! GPS velocities in Earth frame
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   EM         ! Model refraction angle
!
Real(Double), Intent(Out) :: &
   FE_ZN(6)   ! D(E)/D(ZN)
!
Integer, Intent(Out)      :: &
   Stat       ! Error status
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   S = 15.0_Double         ! Ray integration step
!!!   S = 10.0_Double         ! Ray integration step
!
! Local Scalars:
!
! --- Geometric parameters
!
Type(Cartesian) :: XT      ! Transmitter position
Type(Cartesian) :: XR      ! Receiver position
Type(Cartesian) :: XN      ! Final ray point
Type(Cartesian) :: UT      ! Ray direction at transmitter
Type(Cartesian) :: UN      ! Ray direction at XN
Type(Cartesian) :: VT      ! Transmitter velocity [km/s]
Type(Cartesian) :: VR      ! Receiver velocity [km/s]
Real(Double)    :: PM      ! Model impact parameter
!
! Local Arrays:
!
! --- Geometric parameters
!
Real(Double) :: Z0(6)      ! Initial ray point (X0,U0)
Real(Double) :: ZN(6)      ! Final ray point (XN,UN)
!
! --- Ray derivatives
!
Real(Double) :: ZN_Z0(6,6) ! d(XN,UN)/d(XT,UT)
Real(Double) :: E_Z0(6)    ! d(Eps)/d(Z0)
Real(Double) :: E_ZN(6)    ! d(Eps)/d(ZN)
Real(Double) :: P_Z0(6)    ! d(P)/d(Z0)
Real(Double) :: P_ZN(6)    ! d(P)/d(ZN)
Real(Double) :: FZ0_P(6)   ! D(Z0)/D(P)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. FINDING RAY DIRECTION FOR GIVEN IMPACT PARAMETER
!----------------------------------------------------------


!--- 1.1. Setting transmitter and receiver
!---      coordinates/velocities

XT = ERGPS
VT = EVGPS
XR = ERLEO
VR = EVLEO


!--- 1.2. Finding ray

Call Find_Ray_P &
  (XT,      & ! <-- Transmitter position
   VT,      & ! <-- Transmitter velocity
   XR,      & ! <-- Receiver position
   VR,      & ! <-- Receiver velocity
   S,       & ! <-- Integration step
   ERLC,    & ! <-- Local curvature center
   PRO,     & ! <-- Impact parameter [km]
   UT,      & ! --> Transmitter ray direction
   XN,      & ! --> Ray point nearest to receiver
   UN,      & ! --> Ray direction at XN
   Stat)      ! --> Error status

If (Stat /= 0) then
   Return
End If


!----------------------------------------------------------
! 2. REINITIALIZATION OF ECHAM
!----------------------------------------------------------

Where(PStat(:,:,:) == fst_Initialized)
   PStat (:,:,:) = fst_Allocated
End Where
Where(PAStat(:,:,:) == fst_Initialized)
   PAStat(:,:,:) = fst_Allocated
End Where
Where(RAStat(:,:,:) == fst_Initialized)
   RAStat(:,:,:) = fst_Allocated
End Where


!----------------------------------------------------------
! 3. RUNNIG ADJOINT RAY-TRACER
!----------------------------------------------------------


!--- 3.1. Running adjoint ray-tracer

Call Ray_Trace_adj &
  (XT,        & ! <-- Transmitter position
   UT,        & ! <-- Transmitter ray direction
   XR,        & ! <-- Receiver position
   S,         & ! <-- Integration step
   XN,        & ! --> Ray point nearest to receiver
   UN,        & ! --> Ray direction at XN
   ZN_Z0,     & ! --> d(XN,UN)/d(XT,UT)
   Stat)        ! --> Error status


!--- 3.2. Calculation of ray derivatives

Z0(1:3) = XT - ERLC
Z0(4:6) = UT
ZN(1:3) = XN - ERLC
ZN(4:6) = UN

Call Ray_Refraction_adj &
  (Z0,        & ! <-- Initial ray point (X0,U0)
   ZN,        & ! <-- Final ray point (XN,UN)
   VT,        & ! <-- Transmitter velocity [km/s]
   VR,        & ! <-- Receiver velocity [km/s]
   EM,        & ! --> Refraction angle [rad]
   PM,        & ! --> Impact parameter [km]
   E_Z0,      & ! --> d(Eps)/d(Z0)
   E_ZN,      & ! --> d(Eps)/d(ZN)
   P_Z0,      & ! --> d(P)/d(Z0)
   P_ZN)        ! --> d(P)/d(ZN)

Call Ray_Derivatives &
  (Z0,        & ! <-- Initial ray point (X0,U0)
   ZN,        & ! <-- Final ray point (XN,UN)
   ZN_Z0,     & ! <-- d(ZN)/d(Z0)
   P_Z0,      & ! <-- d(P)/d(Z0)
   P_ZN,      & ! <-- d(P)/d(ZN)
   E_Z0,      & ! <-- d(E)/d(Z0)
   E_ZN,      & ! <-- d(D)/d(ZN)
   FZ0_P,     & ! --> D(Z0)/D(P)
   FE_ZN)       ! --> D(E)/D(ZN)


End Subroutine ECHAM_Refraction



!!!TEST
Subroutine Print_Matrix(A,Header)
   Implicit None
   Real(Double), Intent(In) :: &
      A(:,:)
   Character(Len=*)         :: &
      Header
   Integer :: i
   Write(*,'(3A)') '---- ', Header, ' ----'
   Do i=1,Size(A,1)
      Write(*,'(10(:ES13.6:1X))') A(i,:)
   End Do
End Subroutine Print_Matrix



End Module ECHAM_Rays_adj
