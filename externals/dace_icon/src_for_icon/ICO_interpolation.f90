!
!+ GNSS Radio occultation operator: Interpolation on icosahedral-hexagonal grid
!
MODULE ICO_interpolation
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Interpolation on a icosahedral-hexagonal grid.
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
!  Merge some optimizations from mo_grid_intpol
!  Setup_Interpolation: increase threshold to avoid excessive warnings
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  turn error message into warning (inaccurate triangle search)
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
! Module ICO_interpolation
!
! Interpolation on a icosahedral-hexagonal grid.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 05 Dec 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, Pi, dtr, rtd
!
Use mo_exception, only: &
    finish
!
Use ICO_grid, only: gf, &
! Imported Scalar Variables:
!   nd,             &
!   ni, ni2, ni3,   &
!   ng1s,   ng1e,   &
!   ng1sm1, ng1ep1, &
!   ng1sm2, ng1ep2, &
!   ng2s,   ng2e,   &
!   ng2sm1, ng2ep1, &
!   ng2sm2, ng2ep2, &
! Imported Array Variables:
    nspoke
!   xnglob,    &
!   xn,        &
!   rlon,      &
!   rlat
use mo_exception, only: message
!----------------------------------------------------------
Implicit None
Private
Public :: setup_interpolation
!----------------------------------------------------------
!
Contains

!==========================================================
Subroutine Setup_Interpolation &
  (NP,            & ! <-- Number of interpolation points
   slon,          & ! <-- Point longitudes [deg]
   slat,          & ! <-- Point latitudes [deg]
   jg,            & ! --> Triangle mesh indices
   SWS)             ! --> Symplectic weights for triangle vertices
!
! Find triangle grid meshes containing points
! and corresponding symplectic weights.
!----------------------------------------------------------
! Method:
!   Binary search.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 06 Dec 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   NP           ! Number of interpolation points
!
Real(Double), Intent(In) :: &
   slon(NP)     ! Point longitude [deg]
!
Real(Double), Intent(In) :: &
   slat(NP)     ! Point latitude [deg]
!
! Output arguments:
!
Integer, Intent(Out) :: &
   jg(4,NP)     ! Triangle mesh indices
                ! jg(1:3) - (j1,j2,jd) - top vertex indices
                ! jg(4)   - mt - triangle number
!
Real(Double), Intent(Out) :: &
   SWS(3,NP)    ! Symplectic weights for triangle vertices
!----------------------------------------------------------
! Local Parameters:
!
!
! Local Scalars:
!
Integer :: ji         ! Binary search index
Integer :: IP         ! Interpolation points index
Integer :: ipa        ! Pole-/antipoleward index
Integer :: iv         ! Vertex index
!
! Local Arrays:
!
Integer :: &
   JICO(2,3,2)        ! Icosahedral triangle vertices
                      ! [pole-/antipoleward,vertex,index]
Integer :: jd(NP)     ! Diamond index
Integer :: kpa(NP)    ! Pole-/antipoleward triangle discriminator
!
Real(Double) :: &
   zr(3,NP)           ! Cartesian coordinates of points
Integer      :: &
   jb(3,2,NP),      & ! Vertex indices of big triangle
   js(3,2,NP),      & ! Vertex indices of sub-triangle
                      !    (Top/Left/Right,j1/j2)
!  jtv(NP)            ! Top vertex index
   jtv                ! Top vertex index
Real(Double) :: &
!  rtv(gf% nd,2,3,3),& ! Icosahedral triangle vertices
!                      ! [diamond, pole-/antipoleward,vertex,component]
   rtv_(gf% nd,2,3,3),&! Icosahedral triangle vertices
                       ! [diamond, pole-/antipoleward,component,vertex]
   SWB(gf% nd,2,3,NP),&! Symplectic weights of vertices
                       ! [diamond, pole-/antipoleward,vertex,point]
   SWMB(gf% nd,2,NP)   ! Minimum symplectic weight of triangle
                       ! [diamond,pole-/antipoleward,point]
Real(Double) :: &
   SWM(NP)            ! Minimum symplectic weights
                      ! of triangles containing points
Integer      :: &
!  ib(2,NP)           ! Triangle indices (jd, kpa)
   ib_s(2)            ! Triangle indices (jd, kpa)
Real(Double) :: &
!  spv(3,NP)          ! Scalar product for vertices
   spv_s(3)           ! Scalar product for vertices
Integer      :: &
   Spk(2:3,2,6)       ! Spokes for 2nd and 3rd vertices
!----------------------------------------------------------


!----------------------------------------------------------
! 1. ICOSAHEDRAL TRIANGLE SEARCH
!----------------------------------------------------------


!--- 1.1. Interpolation points in Cartesian coordinates

Do IP=1,NP
   zr(:,IP) = (/ COS(slon(IP)*dtr)*COS(slat(IP)*dtr), &
                 SIN(slon(IP)*dtr)*COS(slat(IP)*dtr), &
                 SIN(slat(IP)*dtr) /)
End Do


!--- 1.2. Icosahedral triangle vertices

JICO(:,:,:) = Reshape &
      (Source =(/ (/    0,      1/),(/gf%ni,      1/), (/    0,gf%ni+1/),     &
                  (/gf%ni,gf%ni+1/),(/    0,gf%ni+1/), (/gf%ni,    1  /) /),  &
       Shape = (/ 2, 3, 2 /), Order = (/ 3, 2, 1 /) )

Do ipa=1,2
   Do iv=1,3
      rtv_(:,ipa,:,iv) = &
           Transpose(gf%xnglob(JICO(ipa,iv,1),JICO(ipa,iv,2),:,:))
   End Do
End Do


!--- 1.3. Icosahedral triangle search

Call Symplectic_Weights_ICO   &
  (zr(:,1:NP),                & ! <-- Point in 3D space
   rtv_(1:gf% nd,1:2,1:3,1),  & ! <-- 1st vertices
   rtv_(1:gf% nd,1:2,1:3,2),  & ! <-- 2nd vertices
   rtv_(1:gf% nd,1:2,1:3,3),  & ! <-- 3rd vertices
   SWB (1:gf% nd,1:2,1:3,1:NP)) ! --> Array of symplectic weights

SWMB(1:gf% nd,1:2,1:NP) = MinVal(SWB(1:gf% nd,1:2,:,1:NP),Dim=3)

Do IP=1,NP
   ib_s(:)  = MaxLoc(SWMB(1:gf% nd,1:2,IP))
   jd(IP)   = ib_s(1)
   kpa(IP)  = ib_s(2)
   jg(3,IP) = jd(IP)
End Do


!----------------------------------------------------------
! 2. TRIANGLE SEARCH
!----------------------------------------------------------

!--- 2.1. Initialization of search

Do IP=1,NP
   jb(:,:,IP) = JICO(kpa(IP),:,:)
End Do


!--- 2.2. Binary search

If (gf% nir > 3) call finish('Setup_Interpolation','nir > 3')

If (gf% nir == 3) then

   Call Search_Triangle &
     (NP,               & ! <-- Number of interpolation points
      zr,               & ! <-- Cartesian coordinates of point
      9,                & ! <-- Number of subtriangles
      jd,               & ! <-- Diamond index
      jb,               & ! <-- Vertex indices of big triangle
      js,               & ! --> Vertex indices of sub-triangle
      SWS)                ! --> Symplectic weights for sub-triangle

   jb(:,:,:) = js(:,:,:)

End If

Do ji = 1, gf% ni2

   Call Search_Triangle &
     (NP,               & ! <-- Number of interpolation points
      zr,               & ! <-- Cartesian coordinates of point
      4,                & ! <-- Number of subtriangles
      jd,               & ! <-- Diamond index
      jb,               & ! <-- Vertex indices of big triangle
      js,               & ! --> Vertex indices of sub-triangle
      SWS)                ! --> Symplectic weights for sub-triangle

   jb(:,:,:) = js(:,:,:)

End Do

SWM(:) = MinVal(SWS(:,:),Dim=1)

If (MinVal(SWM(:)) < -4000*Epsilon(SWM)) then
   Write(*,'(A)') 'Warning: inaccurate triangle search !!!'
   Do IP=1,NP
      If (SWM(IP) < -4000*Epsilon(SWM)) then
         Write(*,'(I6,2(A,F8.3))') &
            IP, 'Slon = ', Slon(IP), '   Slat = ', Slat(IP)
         Write(*,'(A,3(ES11.4,1X))') 'SWS = ', SWS(:,IP)
      End If
   End Do
   call message ('Setup_Interpolation (ICO_interpolation)',&
                'Warning: inaccurate triangle search ')
End If


!--- 2.3. Setting top to nearest grid point.

Do IP=1,NP

   Do iv=1,3
      spv_s(iv) = Sum(zr(:,IP)* gf% xnglob(jb(iv,1,IP),jb(iv,2,IP),:,jd(IP)))
   End Do

   jtv        = Sum(MaxLoc(spv_s(:)))
   js(:,:,IP) = CShift(js(:,:,IP), Shift=jtv-1, Dim=1)
   jb(:,:,IP) = CShift(jb(:,:,IP), Shift=jtv-1, Dim=1)
   SWS(:,IP)  = CShift(SWS(:,IP),  Shift=jtv-1)

   jg(1:2,IP) = jb(1,1:2,IP)

End Do


!--- 2.4. Determination of triangle number

Spk(2,:,:) = nspoke(:,:)
Spk(3,:,:) = CShift(nspoke(:,:), Shift=1, Dim=2)

Do IP=1,NP
   jg(4,IP) = Sum(MinLoc(nspoke(1,:), &
                Mask = ((Spk(2,1,:) == jb(2,1,IP) - jb(1,1,IP)) .and. &
                        (Spk(2,2,:) == jb(2,2,IP) - jb(1,2,IP)) .and. &
                        (Spk(3,1,:) == jb(3,1,IP) - jb(1,1,IP)) .and. &
                        (Spk(3,2,:) == jb(3,2,IP) - jb(1,2,IP)))))

End Do

End Subroutine Setup_Interpolation



!==========================================================
Subroutine Search_Triangle &
  (NP,               & ! <-- Number of interpolation points
   zr,               & ! <-- Cartesian coordinates of point
   NS,               & ! <-- Number of subtriangles
   jd,               & ! <-- Diamond index
   jb,               & ! <-- Vertex indices of big triangle
   js,               & ! --> Vertex indices of sub-triangle
   SWS)                ! --> Symplectic weights for sub-triangle
!
! Determine sub-triangle containing point. The sub-triangle
! is one of four or nine which partition big triangle.
!----------------------------------------------------------
! Method:
!   Search for big triangle with nearest center;
!   search for triangle with maximum symplectic weight.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 05 Dec 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   NP               ! Number of interpolation points
!
Real(Double), Intent(In) :: &
   zr(3,NP)         ! Cartesian coordinates of point
!
Integer, Intent(In) :: &
   NS               ! Number of subtriangles
!
Integer, Intent(In) :: &
   jd(NP)           ! Diamond index
!
Integer, Intent(In) :: &
   jb(3,2,NP)       ! Vertex indices of big triangle
                    ! (Top/Left/Right,j1/j2,Point)
!
! Output arguments:
!
Integer, Intent(Out) :: &
   js(3,2,NP)       ! Vertex indices of sub-triangle
                    ! (Top/Left/Right,j1/j2,Point)
!
Real(Double), Intent(Out) :: &
   SWS(3,NP)        ! Symplectic weights for sub-triangle
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: IP    ! Point number
Integer      :: iv    ! Vertex index
Integer      :: id    ! Diamond index
!
! Local Arrays:
!
Integer      ::  &
   ide(NP)          ! Subtriangle index
Real(Double) :: &
   SW(NS,3,NP),   & ! Symplectic weights
!  SWM(NS,NP)       ! Minimum symplectic weight
   SWM(NS)          ! Minimum symplectic weight
Integer :: &
   j(10,2,NP)       ! Subtriangle vertex indices
Real(Double) ::  &
!  rtv(NS,3,3,NP)   ! Subtriangle vertices (original dimension order)
                    ! (Sub-triangle,vertex,component,point)
   rtv_x(NS,3,NP,3) ! Subtriangle vertices (modified dimension order [ha])
                    ! (Sub-triangle,component,point,vertex)
Integer, Pointer :: &
   SI(:,:)          ! Pointer for sub-triangle index array
!----------------------------------------------------------
Integer, Target :: &
   SI9(9,3) = &     ! Sub-triangle indices
      Reshape( (/ 1,  2,  3,       &  ! D1
                  2,  4,  5,       &  ! D2
                  4,  7,  8,       &  ! D3
                  5,  8,  9,       &  ! D4
                  6,  9, 10,       &  ! D5
                  3,  5,  6,       &  ! D6
                  5,  3,  2,       &  ! D7
                  8,  5,  4,       &  ! D8
                  9,  6,  5  /),   &  ! D9
               Shape = (/ 9, 3 /), &
               Order = (/ 2, 1 /) )
!
! Poleward triangle                 1(1)
!                                   / \
!                                  / D1\
!                                 /     \
!                                2-------3
!                               / \     / \
!                              / D2\ D7/ D6\
!                             /     \ /     \
!                            4-------5-------6
!                           / \     / \     / \
!                          / D3\ D8/ D4\ D9/ D5\
!                         /     \ /     \ /     \
!                       7(2)-----8-------9-----10(3)
!
!-----------------------------------------------------------------------
Integer, Target :: &
   SI4(4,3) = &     ! Sub-triangle indices
      Reshape( (/ 1,  2,  3,       &  ! D1
                  2,  4,  5,       &  ! D2
                  3,  5,  6,       &  ! D3
                  5,  3,  2 /),    &  ! D4
               Shape = (/ 4, 3 /), &
               Order = (/ 2, 1 /) )
!
! Poleward triangle                 1(1)
!                                   / \
!                                  / D1\
!                                 /     \
!                                2------ 3
!                               / \     / \
!                              / D2\ D4/ D3\
!                             /     \ /     \
!                           4(2)-----5------6(3)
!
!-----------------------------------------------------------------------
! Antipoleward triangle: top is oriented away from pole,
!                        lower left - eastern, lower right - western.
!=======================================================================

!CDIRR ON_ADB(j)
!CDIRR ON_ADB(rtv_x)

!----------------------------------------------------------
! 1. SUB-TRIANGLE DEFINITION
!----------------------------------------------------------

!--- 1.1. Vertices

Select Case (NS)

   Case (4)

      j( 1,:,:)  = jb(1,:,:)
      j( 4,:,:)  = jb(2,:,:)
      j( 6,:,:)  = jb(3,:,:)
      j( 2,:,:)  = (j(1,:,:) + j(4,:,:))/2
      j( 3,:,:)  = (j(1,:,:) + j(6,:,:))/2
      j( 5,:,:)  = (j(4,:,:) + j(6,:,:))/2

      SI => SI4

   Case (9)

      j( 1,:,:)  = jb(1,:,:)
      j( 7,:,:)  = jb(2,:,:)
      j(10,:,:)  = jb(3,:,:)
      j( 2,:,:)  =  j(1,:,:) +   (j( 7,:,:) - j(1,:,:))/3
      j( 4,:,:)  =  j(1,:,:) + 2*(j( 7,:,:) - j(1,:,:))/3
      j( 3,:,:)  =  j(1,:,:) +   (j(10,:,:) - j(1,:,:))/3
      j( 6,:,:)  =  j(1,:,:) + 2*(j(10,:,:) - j(1,:,:))/3
      j( 8,:,:)  =  j(7,:,:) +   (j(10,:,:) - j(7,:,:))/3
      j( 9,:,:)  =  j(7,:,:) + 2*(j(10,:,:) - j(7,:,:))/3
      j( 5,:,:)  = (j(8,:,:) + j(3,:,:))/2

      SI => SI9

   Case Default

      Write(*,'(A,I0/A)') 'Search_Triangle: ns = ', ns, &
                          '   must be 4 or 9'
      Stop

End Select

!CDIRR ON_ADB(SI)

!--- 1.2. Centres

Do IP=1,NP
   Do iv=1,3
      Do id=1,NS
         rtv_x(id,:,IP,iv) = &
            gf% xnglob(j(SI(id,iv),1,IP),j(SI(id,iv),2,IP),:,jd(IP))
      End Do
   End Do
End Do


!----------------------------------------------------------
! 2. SUB-TRIANGLE SEARCH
!----------------------------------------------------------

!--- 2.1. Search for triangle with maximum symplectic weight

Call Symplectic_Weights_Small &
  (zr(:,:),         & ! <-- Point in 3D space
   rtv_x(:,:,:,1),  & ! <-- 1st vertices
   rtv_x(:,:,:,2),  & ! <-- 2nd vertices
   rtv_x(:,:,:,3),  & ! <-- 3rd vertices
   SW(:,:,:))         ! --> Array of symplectic weights

Do IP=1,NP
   SWM(:)    = MinVal(SW(:,:,IP),Dim=2)
   ide(IP)   = Sum(MaxLoc(SWM(:)))
   SWS(:,IP) = SW(ide(IP),:,IP)
End Do


!--- 2.2. Setting vertices

Do IP=1,NP
   js(:,:,IP) = j(SI(ide(IP),:),:,IP)
End Do

End Subroutine Search_Triangle



!==========================================================
Subroutine Symplectic_Weights_Small &
  (z,       & ! <-- Point in 3D space
   x1,      & ! <-- 1st vertices
   x2,      & ! <-- 2nd vertices
   x3,      & ! <-- 3rd vertices
   SW)        ! --> Array of symplectic weights
!
! Calculation of symplectic weights for
! multiple (point/small triangle set) pairs.
!----------------------------------------------------------
! Method:
!            (z,x2,x3)     (x1,z,x3)     (x1,x2,z)
! SW(1:3) = ------------, ------------, ------------
!                N             N             N
! N is norming constant: Sum(SW(1:3)) = 1
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 05 Dec 2001 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   z(:,:)     ! Points
              ! [Component,point number]
!
Real(Double), Intent(In) :: &
   x1(:,:,:)  ! 1st symplex direction
              ! [Triangle number,component,point]
!
Real(Double), Intent(In) :: &
   x2(:,:,:)  ! 2nd symplex direction
!
Real(Double), Intent(In) :: &
   x3(:,:,:)  ! 3rd symplex direction
!
! Output arguments:
!
Real(Double) :: &
   SW(:,:,:)  ! Array of symplectic weights
              ! [Triangle number,vertex,point]
!----------------------------------------------------------
! Local Scalars:
!
Integer :: IP  ! Point index
!
! Local Arrays:
!
Real(Double) :: &
!  D(Size(SW,1),Size(SW,3))   ! Symplectic denominator
   D(Size(SW,1))              ! Symplectic denominator
!----------------------------------------------------------


Do IP=1,Size(SW,3)

   SW(:,1,IP) = (  z(1,IP)*x2(:,2,IP)*x3(:,3,IP)    &
                 - z(1,IP)*x2(:,3,IP)*x3(:,2,IP) +  &
                   z(2,IP)*x2(:,3,IP)*x3(:,1,IP)    &
                 - z(2,IP)*x2(:,1,IP)*x3(:,3,IP) +  &
                   z(3,IP)*x2(:,1,IP)*x3(:,2,IP)    &
                 - z(3,IP)*x2(:,2,IP)*x3(:,1,IP))

   SW(:,2,IP) = (  x1(:,1,IP)*z(2,IP)*x3(:,3,IP)    &
                 - x1(:,1,IP)*z(3,IP)*x3(:,2,IP) +  &
                   x1(:,2,IP)*z(3,IP)*x3(:,1,IP)    &
                 - x1(:,2,IP)*z(1,IP)*x3(:,3,IP) +  &
                   x1(:,3,IP)*z(1,IP)*x3(:,2,IP)    &
                 - x1(:,3,IP)*z(2,IP)*x3(:,1,IP))

   SW(:,3,IP) = (  x1(:,1,IP)*x2(:,2,IP)*z(3,IP)    &
                 - x1(:,1,IP)*x2(:,3,IP)*z(2,IP) +  &
                   x1(:,2,IP)*x2(:,3,IP)*z(1,IP)    &
                 - x1(:,2,IP)*x2(:,1,IP)*z(3,IP) +  &
                   x1(:,3,IP)*x2(:,1,IP)*z(2,IP)    &
                 - x1(:,3,IP)*x2(:,2,IP)*z(1,IP))

#if defined (__SX__)
   ! Explicit expansion of sum() for better vectorization
   D(:)       = SW(:,1,IP) + SW(:,2,IP) + SW(:,3,IP)
#else
   D(:)       = Sum (SW(:,1:3,IP), Dim=2)
#endif

   SW(:,1,IP) = SW(:,1,IP)/D(:)
   SW(:,2,IP) = SW(:,2,IP)/D(:)
   SW(:,3,IP) = SW(:,3,IP)/D(:)

End Do


End Subroutine Symplectic_Weights_Small



!==========================================================
Subroutine Symplectic_Weights_ICO &
  (z,       & ! <-- Point in 3D space
   x1,      & ! <-- 1st vertices
   x2,      & ! <-- 2nd vertices
   x3,      & ! <-- 3rd vertices
   SW)        ! --> Array of symplectic weights
!
! Calculation of symplectic weights for icosahedral
! triangles and multiple points
!----------------------------------------------------------
! Method:
!            (z,x2,x3)     (x1,z,x3)     (x1,x2,z)
! SW(1:3) = ------------, ------------, ------------
!            (x1,x2,x3)    (x1,x2,x3)    (x1,x2,x3)
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 06 Dec 2001 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   z(:,:)       ! Points
                ! [Component,point number]
!
Real(Double), Intent(In) :: &
   x1(:,:,:)    ! 1st symplex direction
                ! [Diamond,pole-/antipoleward,component]
!
Real(Double), Intent(In) :: &
   x2(:,:,:)    ! 2nd symplex direction
!
Real(Double), Intent(In) :: &
   x3(:,:,:)    ! 3rd symplex direction
!
! Output arguments:
!
Real(Double) :: &
   SW(:,:,:,:)  ! Array of symplectic weights
                ! [Diamond,pole-/antipoleward,vertex,point]
!----------------------------------------------------------
! Local Scalars:
!
Integer :: IP  ! Point index
!
! Local Arrays:
!
Real(Double) :: &
!  D(Size(SW,1),Size(SW,2),Size(SW,4))   ! Symplectic denominator
   D(Size(SW,1),Size(SW,2))              ! Symplectic denominator
!----------------------------------------------------------

D(:,:)    =  ( x1(:,:,1)*x2(:,:,2)*x3(:,:,3)    &
             - x1(:,:,1)*x2(:,:,3)*x3(:,:,2) +  &
               x1(:,:,2)*x2(:,:,3)*x3(:,:,1)    &
             - x1(:,:,2)*x2(:,:,1)*x3(:,:,3) +  &
               x1(:,:,3)*x2(:,:,1)*x3(:,:,2)    &
             - x1(:,:,3)*x2(:,:,2)*x3(:,:,1))

Do IP=1,Size(SW,4)

   SW(:,:,1,IP) = (  z(1,IP)*x2(:,:,2)*x3(:,:,3)    &
                   - z(1,IP)*x2(:,:,3)*x3(:,:,2) +  &
                     z(2,IP)*x2(:,:,3)*x3(:,:,1)    &
                   - z(2,IP)*x2(:,:,1)*x3(:,:,3) +  &
                     z(3,IP)*x2(:,:,1)*x3(:,:,2)    &
                   - z(3,IP)*x2(:,:,2)*x3(:,:,1)) / D(:,:)

   SW(:,:,2,IP) = (  x1(:,:,1)*z(2,IP)*x3(:,:,3)    &
                   - x1(:,:,1)*z(3,IP)*x3(:,:,2) +  &
                     x1(:,:,2)*z(3,IP)*x3(:,:,1)    &
                   - x1(:,:,2)*z(1,IP)*x3(:,:,3) +  &
                     x1(:,:,3)*z(1,IP)*x3(:,:,2)    &
                   - x1(:,:,3)*z(2,IP)*x3(:,:,1)) / D(:,:)

   SW(:,:,3,IP) = (  x1(:,:,1)*x2(:,:,2)*z(3,IP)    &
                   - x1(:,:,1)*x2(:,:,3)*z(2,IP) +  &
                     x1(:,:,2)*x2(:,:,3)*z(1,IP)    &
                   - x1(:,:,2)*x2(:,:,1)*z(3,IP) +  &
                     x1(:,:,3)*x2(:,:,1)*z(2,IP)    &
                   - x1(:,:,3)*x2(:,:,2)*z(1,IP)) / D(:,:)

End Do


End Subroutine Symplectic_Weights_ICO



!==========================================================
Subroutine Interpolation &
  (NP,        & ! <-- Number of interpolation points
   px,        & ! <-- Gridded field
   SWS,       & ! <-- Symplectic weights for triangle vertices of points
   jg,        & ! <-- Index of the triangle containing point
   pxi)         ! --> Interpolated field
!
! Interpolation on icosahedral-hexagonal grid.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   0.0   | 18 Oct 2001 | Original code by Baumgartner.
!   1.0   | 05 Dec 2001 | Adaptation for use
!         |             | with ECHAM library.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   NP                    ! Number of interpolation points
!
Real(Double), Intent(In) ::   &
   px(gf% lbg(1):gf% ubg(1),  &  ! Gridded field
      gf% lbg(2):gf% ubg(2),  &
      gf% nd)
!
Real(Double), Intent(In) :: &
   SWS(3,NP)             ! Symplectic weights for triangle vertices
!
Integer, Intent(In) :: &
   jg(4,NP)              ! Index of the triangle containing point:
                         !  j1,j2,jd, mt
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   pxi(NP)               ! Interpolated field
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: zx1, zx2, zx3   ! values at three gridpoints
Integer      :: IP              ! Interpolation point index
Integer      :: j1, j2, jd      ! Gridpoint indices
Integer      :: m1, m2          ! Triangle indices
!----------------------------------------------------------


Do IP = 1, NP

   j1    = jg(1,IP)
   j2    = jg(2,IP)
   jd    = jg(3,IP)

   m1    = jg(4,IP)
   m2    = MOD (m1,6) + 1

   zx1 = px(j1,j2,jd)
   zx2 = px(j1+nspoke(1,m1), j2+nspoke(2,m1), jd)
   zx3 = px(j1+nspoke(1,m2), j2+nspoke(2,m2), jd)

   pxi(IP) = SWS(1,IP)*zx1 + SWS(2,IP)*zx2 + SWS(3,IP)*zx3

End Do

End Subroutine Interpolation



End Module ICO_interpolation
