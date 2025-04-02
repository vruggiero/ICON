!
!+ GNSS Radio occultation observation operator: Moist-air propagation model
!
MODULE Occ_MPM93
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Moist-air propagation model.
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
! Module Occ_MPM93
!
! Moist-air propagation model.
!----------------------------------------------------------
! (C) Copyright 2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 10 May 2001 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double
!----------------------------------------------------------
Implicit None
Private
Public :: N93AIR
!----------------------------------------------------------
!
Contains

!==========================================================
Subroutine N93AIR &
  (FQ,      & ! <-- Frequency [Hz]
   PP,      & ! <-- Pressure [mbar]
   T,       & ! <-- Temperature [K]
   QH,      & ! <-- Specific humidity [kg/kg]
   ZN,      & ! --> Dispersive part of refractivity [N-units]
   EN0)       ! --> Non-dispesive part of refractivity [N-units]
!
! Computation of complex refracitivity for given
! pressure, temperature, and humidity.
!----------------------------------------------------------
! Method:
!   Described in:
!   Liebe et al., "Propagation Modeling of moist air and suspended
!     water/ice particles at frequencies below 1000 GHz,
!     AGARD CP-May93, Paper 3/1-10.
!----------------------------------------------------------
! (C) Copyright 2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   |        1993 | Original Fortran-77 code by
!         |             | Dr. G. Hufford and Dr. H. Liebe,
!         |             | NTIA/ITS.S3
!   2.0   | 11 May 2001 | Fortran-90 version by
!         |             | M.E.Gorbunov.
!----------------------------------------------------------
! Modules used:
!
Use Occ_Meteoprofiles, only: &
! Imported Parameters:
    aq, bq
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
!
! Input arguments:
!
Real(Double), Intent(In) :: &
   FQ         ! Frequency [GHz]
!
Real(Double), Intent(In) :: &
   PP         ! Pressure [mbar]
!
Real(Double), Intent(In) :: &
   T          ! Temperature [K]
!
Real(Double), Intent(In) :: &
   QH         ! Specific humidity [kg/kg]
!
! Output arguments:
!
Complex(Double), Intent(Out) :: &
   ZN         ! Dispersive part of refractivity:
              !    Propagation delay is 3.336e6*Real(ZN) ps/km
              !    Specific attenuation is 0.1820e-3*FQ*Imag(ZN) dB/km
!
Real(Double), Intent(Out) :: &
   EN0        ! Non-dispesive part of refractivity
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: BFFS(79) =      &
     (/ 50.474238D0, 50.987749D0, 51.503350D0, 52.021410D0, 52.542394D0,   &
        53.066907D0, 53.595749D0, 54.130000D0, 54.671159D0, 55.221367D0,   &
        55.783802D0, 56.264775D0, 56.363389D0, 56.968206D0, 57.612484D0,   &
        58.323877D0, 58.446590D0, 59.164207D0, 59.590983D0, 60.306061D0,   &
        60.434776D0, 61.150560D0, 61.800154D0, 62.411215D0, 62.486260D0,   &
        62.997977D0, 63.568518D0, 64.127767D0, 64.678903D0, 65.224071D0,   &
        65.764772D0, 66.302091D0, 66.836830D0, 67.369598D0, 67.900867D0,   &
        68.431005D0, 68.960311D0,118.750343D0,368.498350D0,424.763124D0,   &
       487.249370D0,715.393150D0,773.839675D0,834.145330D0,   &
        22.235080D0, 67.803960D0,119.995940D0,183.310091D0,321.225644D0,   &
       325.152919D0,336.222601D0,380.197372D0,390.134508D0,437.346667D0,   &
       439.150812D0,443.018295D0,448.001075D0,470.888947D0,474.689127D0,   &
       488.491133D0,503.568532D0,504.482692D0,547.676440D0,552.020960D0,   &
       556.936002D0,620.700807D0,645.866155D0,658.005280D0,752.033227D0,   &
       841.053973D0,859.962313D0,899.306675D0,902.616173D0,906.207325D0,   &
       916.171582D0,923.118427D0,970.315022D0,987.926764D0,1780.00000D0   /)
!
Real(Double), Parameter :: A1(79) =   &
     (/    0.094D-6,   0.246D-6,   0.608D-6,   1.414D-6,   3.102D-6,   &
           6.410D-6,  12.470D-6,  22.800D-6,  39.180D-6,  63.160D-6,   &
          95.350D-6,  54.890D-6, 134.400D-6, 176.300D-6, 214.100D-6,   &
         238.600D-6, 145.700D-6, 240.400D-6, 211.200D-6, 212.400D-6,   &
         246.100D-6, 250.400D-6, 229.800D-6, 193.300D-6, 151.700D-6,   &
         150.300D-6, 108.700D-6,  73.350D-6,  46.350D-6,  27.480D-6,   &
          15.300D-6,   8.009D-6,   3.946D-6,   1.832D-6,   0.801D-6,   &
           0.330D-6,   0.128D-6,  94.500D-6,   6.790D-6,  63.800D-6,   &
          23.500D-6,   9.960D-6,  67.100D-6,  18.000D-6,               &
           0.1130D-1,  0.0012D-1,  0.0008D-1,  2.4200D-1,  0.0483D-1,  &
           1.4990D-1,  0.0011D-1, 11.5200D-1,  0.0046D-1,  0.0650D-1,  &
           0.9218D-1,  0.1976D-1, 10.3200D-1,  0.3297D-1,  1.2620D-1,  &
           0.2520D-1,  0.0390D-1,  0.0130D-1,  9.7010D-1, 14.7700D-1,  &
         487.4000D-1,  5.0120D-1,  0.0713D-1,  0.3022D-1,239.6000D-1,  &
           0.0140D-1,  0.1472D-1,  0.0605D-1,  0.0426D-1,  0.1876D-1,  &
           8.3410D-1,  0.0869D-1,  8.9720D-1,132.1000D-1,22300.0000D-1  /)
!
Real(Double), Parameter :: A2(79) =   &
     (/ 9.694_Double, 8.694_Double, 7.744_Double, 6.844_Double, 6.004_Double, &
        5.224_Double, 4.484_Double, 3.814_Double, 3.194_Double, 2.624_Double, &
        2.119_Double, 0.015_Double, 1.660_Double, 1.260_Double, 0.915_Double, &
        0.626_Double, 0.084_Double, 0.391_Double, 0.212_Double, 0.212_Double, &
        0.391_Double, 0.626_Double, 0.915_Double, 1.260_Double, 0.083_Double, &
        1.665_Double, 2.115_Double, 2.620_Double, 3.195_Double, 3.815_Double, &
        4.485_Double, 5.225_Double, 6.005_Double, 6.845_Double, 7.745_Double, &
        8.695_Double, 9.695_Double, 0.009_Double, 0.049_Double, 0.044_Double, &
        0.049_Double, 0.145_Double, 0.130_Double, 0.147_Double,   &
        2.143_Double, 8.735_Double, 8.356_Double, 0.668_Double, 6.181_Double, &
        1.540_Double, 9.829_Double, 1.048_Double, 7.350_Double, 5.050_Double, &
        3.596_Double, 5.050_Double, 1.405_Double, 3.599_Double, 2.381_Double, &
        2.853_Double, 6.733_Double, 6.733_Double, 0.114_Double, 0.114_Double, &
        0.159_Double, 2.200_Double, 8.580_Double, 7.820_Double, 0.396_Double, &
        8.180_Double, 7.989_Double, 7.917_Double, 8.432_Double, 5.111_Double, &
        1.442_Double,10.220_Double, 1.920_Double, 0.258_Double, 0.952_Double /)
!
Real(Double), Parameter :: A3(79) =   &
     (/ 0.890D-3,  0.910D-3,  0.940D-3,  0.970D-3,  0.990D-3,   &
        1.020D-3,  1.050D-3,  1.070D-3,  1.100D-3,  1.130D-3,   &
        1.170D-3,  1.730D-3,  1.200D-3,  1.240D-3,  1.280D-3,   &
        1.330D-3,  1.520D-3,  1.390D-3,  1.430D-3,  1.450D-3,   &
        1.360D-3,  1.310D-3,  1.270D-3,  1.230D-3,  1.540D-3,   &
        1.200D-3,  1.170D-3,  1.130D-3,  1.100D-3,  1.070D-3,   &
        1.050D-3,  1.020D-3,  0.990D-3,  0.970D-3,  0.940D-3,   &
        0.920D-3,  0.900D-3,  1.630D-3,  1.920D-3,  1.930D-3,   &
        1.920D-3,  1.810D-3,  1.810D-3,  1.820D-3,              &
        2.811D-3,  2.858D-3,  2.948D-3,  3.050D-3,  2.303D-3,   &
        2.783D-3,  2.693D-3,  2.873D-3,  2.152D-3,  1.845D-3,   &
        2.100D-3,  1.860D-3,  2.632D-3,  2.152D-3,  2.355D-3,   &
        2.602D-3,  1.612D-3,  1.612D-3,  2.600D-3,  2.600D-3,   &
        3.210D-3,  2.438D-3,  1.800D-3,  3.210D-3,  3.060D-3,   &
        1.590D-3,  3.060D-3,  2.985D-3,  2.865D-3,  2.408D-3,   &
        2.670D-3,  2.900D-3,  2.550D-3,  2.985D-3, 17.620D-3   /)
!
Real(Double), Parameter :: A4(45:79) =   &
     (/ 4.80_Double,  4.93_Double,  4.78_Double,  5.30_Double,  4.69_Double, &
        4.85_Double,  4.74_Double,  5.38_Double,  4.81_Double,  4.23_Double, &
        4.29_Double,  4.23_Double,  4.84_Double,  4.57_Double,  4.65_Double, &
        5.04_Double,  3.98_Double,  4.01_Double,  4.50_Double,  4.50_Double, &
        4.11_Double,  4.68_Double,  4.00_Double,  4.14_Double,  4.09_Double, &
        5.76_Double,  4.09_Double,  4.53_Double,  5.10_Double,  4.70_Double, &
        5.00_Double,  4.78_Double,  4.94_Double,  4.55_Double, 30.50_Double /)
!
Real(Double), Parameter :: A5(79) =   &
     (/ 0.240D-3,  0.220D-3,  0.197D-3,  0.166D-3,  0.136D-3,   &
        0.131D-3,  0.230D-3,  0.335D-3,  0.374D-3,  0.258D-3,   &
       -0.166D-3,  0.390D-3, -0.297D-3, -0.416D-3, -0.613D-3,   &
       -0.205D-3,  0.748D-3, -0.722D-3,  0.765D-3, -0.705D-3,   &
        0.697D-3,  0.104D-3,  0.570D-3,  0.360D-3, -0.498D-3,   &
        0.239D-3,  0.108D-3, -0.311D-3, -0.421D-3, -0.375D-3,   &
       -0.267D-3, -0.168D-3, -0.169D-3, -0.200D-3, -0.228D-3,   &
       -0.240D-3, -0.250D-3, -0.036D-3,  &
        0._Double, 0._Double, 0._Double, 0._Double, 0._Double,  0._Double,   &
        0.69_Double,  0.69_Double,  0.70_Double,  0.64_Double,  0.67_Double, &
        0.68_Double,  0.69_Double,  0.54_Double,  0.63_Double,  0.60_Double, &
        0.63_Double,  0.60_Double,  0.66_Double,  0.66_Double,  0.65_Double, &
        0.69_Double,  0.61_Double,  0.61_Double,  0.70_Double,  0.70_Double, &
        0.69_Double,  0.71_Double,  0.60_Double,  0.69_Double,  0.68_Double, &
        0.33_Double,  0.68_Double,  0.68_Double,  0.70_Double,  0.70_Double, &
        0.70_Double,  0.70_Double,  0.64_Double,  0.68_Double,  2.00_Double /)
!
Real(Double), Parameter :: A6(79) =   &
     (/ 0.790D-3,  0.780D-3,  0.774D-3,  0.764D-3,  0.751D-3,   &
        0.714D-3,  0.584D-3,  0.431D-3,  0.305D-3,  0.339D-3,   &
        0.705D-3, -0.113D-3,  0.753D-3,  0.742D-3,  0.697D-3,   &
        0.051D-3, -0.146D-3,  0.266D-3, -0.090D-3,  0.081D-3,   &
       -0.324D-3, -0.067D-3, -0.761D-3, -0.777D-3,  0.097D-3,   &
       -0.768D-3, -0.706D-3, -0.332D-3, -0.298D-3, -0.423D-3,   &
       -0.575D-3, -0.700D-3, -0.735D-3, -0.744D-3, -0.753D-3,   &
       -0.760D-3, -0.765D-3,  0.009D-3,  &
        0._Double, 0._Double, 0._Double, 0._Double, 0._Double,  0._Double,   &
        1.00_Double,  0.82_Double,  0.79_Double,  0.85_Double,  0.54_Double, &
        0.74_Double,  0.61_Double,  0.89_Double,  0.55_Double,  0.48_Double, &
        0.52_Double,  0.50_Double,  0.67_Double,  0.65_Double,  0.64_Double, &
        0.72_Double,  0.43_Double,  0.45_Double,  1.00_Double,  1.00_Double, &
        1.00_Double,  0.68_Double,  0.50_Double,  1.00_Double,  0.84_Double, &
        0.45_Double,  0.84_Double,  0.90_Double,  0.95_Double,  0.53_Double, &
        0.78_Double,  0.80_Double,  0.67_Double,  0.90_Double,  5.00_Double /)
!
! Local Scalars:
!
Complex(Double) :: ZNN    ! Dispersive refractivity normalized to frequency
Real(Double)    :: TH     ! 300/T
Real(Double)    :: THV    ! Log(TH)
Real(Double)    :: E      ! Water vapour partial pressure [mb]
Real(Double)    :: P      ! Dry air partial pressure [mb]
Real(Double)    :: GSP    !
Real(Double)    :: AP1    !
Real(Double)    :: AP2    !
Real(Double)    :: Q1     !
Real(Double)    :: Q2     !
Real(Double)    :: Q3     !
Real(Double)    :: BFQ    ! Frequency [GHz]
Integer         :: J      ! Work index
Real(Double)    :: Q      !
Real(Double)    :: AGD    !
!
! Local Arrays:
!
Complex(Double) :: ZH(38)
Real(Double)    :: AH(39:79)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. NON-DISPERSIVE REFRACTIVITY
!----------------------------------------------------------

ZNN = (0._Double, 0._Double)
TH  = 300._Double/T
THV = Log(TH)

E   = PP*QH/(aq + bq*QH)

P   = PP - E
EN0 = 1e-6_Double*(0.2588_Double*P + (0.239_Double + 4.163_Double*TH)*E)*TH


!----------------------------------------------------------
! 2. OXYGEN
!----------------------------------------------------------

GSP = (P + E)*Exp(0.8_Double*THV)*0.56D-3
AP1 = 6.14D-5*P*TH**2
AP2 = P**2*Exp(3.5_Double*THV - 27.2945_Double)
Q1  = P*Exp(3._Double*THV)
Q3  = PP*Exp(0.8_Double*THV)
Q2  = P*Exp(0.8_Double*THV) + 1.1_Double*E*TH

BFQ  = 1e-9_Double*FQ

Do J=1,38
   Q     = Q1*A1(J)*Exp(A2(J)*(1._Double - TH))
   ZH(J) = (Q/BFFS(J))*Cmplx(1._Double, -(A5(J) + A6(J)*TH)*Q3, Double)
   Q     = Sqrt((A3(J)*Q2)**2 + 2.25D-6)
   ZNN   = ZNN + ZH(J)/Cmplx(BFFS(J) - BFQ, -Q, Double)   &
           -Conjg(ZH(J)/Cmplx(BFFS(J) + BFQ, -Q, Double))
End Do

Q2 = P*Exp(0.2_Double*THV) + 1.1_Double*E*TH

DO J=39,44
   Q     = Q1*A1(J)*Exp(A2(J)*(1._Double - TH))
   AH(J) = Q/BFFS(J)
   Q     = Sqrt((A3(J)*Q2)**2 + 2.25D-6)
   ZNN   = ZNN+AH(J)/Cmplx(BFFS(J) - BFQ, -Q, Double)   &
           -AH(J)/Cmplx(BFFS(J) + BFQ, Q, Double)
End Do

If (aimag(ZNN) < 0._Double) then
   ZNN = 0._Double
End If

ZNN = ZNN - AP1/Cmplx(BFQ, GSP, Double) + &
      Cmplx(0._Double, AP2/(1._Double + 1.93D-5*BFQ**1.5_Double), Double)


!----------------------------------------------------------
! 3. WATER VAPOR
!----------------------------------------------------------

AGD = 2.1316D-12/TH

Q1  = E*Exp(3.5_Double*THV)

Do J=45,79
   Q     = Q1*A1(J)*Exp(A2(J)*(1._Double - TH))
   AH(J) = Q/BFFS(J)
   Q     = A3(J)*(P*Exp(A5(J)*THV) + E*A4(J)*Exp(A6(J)*THV))
   Q     = 0.535_Double*Q + Sqrt(0.217_Double*Q**2 + AGD*BFFS(J)**2)
   ZNN   = ZNN + AH(J)/Cmplx(BFFS(J) - BFQ, -Q, Double)   &
         -AH(J)/Cmplx(BFFS(J) + BFQ, Q, Double)
End Do


ZN = 1e-6_Double*BFQ*ZNN


Return


End Subroutine N93AIR



End Module Occ_MPM93
