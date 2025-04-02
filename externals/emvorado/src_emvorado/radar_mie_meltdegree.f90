! Source module for the radar forward operator EMVORADO
!
! ---------------------------------------------------------------
! Copyright (C) 2005-2024, DWD, KIT
! Contact information: ulrich.blahak (at) dwd.de 
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE radar_mie_meltdegree

!------------------------------------------------------------------------------
!
! Description:  Functions for parameterization of the degree
!               of melting of hydrometeors for
!               the model specific interface(s) (at the moment only COSMO)
!               to the EMVORADO libraries for radar reflectivity
!               calculations based on Mie-Theory
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:
!

  USE radar_kind    , ONLY: dp

  USE radar_data_mie, ONLY:  &
       t_mgd_params, particle, quasi_zero

  USE radar_mie_utils, ONLY :  &
       D_of_X_average_1M_exp

!===============================================================================

  IMPLICIT NONE

  PUBLIC

!===============================================================================


!===============================================================================
!===============================================================================


CONTAINS

!===============================================================================
!===============================================================================


  !----------------------------------------------------------------------------------------!
  ! Simple parameterisation of degree of melting of ice and snow particles
  ! as a function of temperature and diameter of particle Dvec(i),
  ! which consists of a simple parameter model, with all parameter set heuristicly
  ! and needed to be validated by more sophisticated melting models.
  !
  ! Returned value of degree of melting at point im,jm,km is estimated from temperature at im,jm,km
  ! as well as a maximum temperature Tmax_XX, where the hydrometeor XX (ice, snow, graupel, hail)
  ! is present within a neighbourhood im +/- neigh, jm +/- neigh.
  !
  ! The temperature Tmax_XX is 2D (i,j) on the model grid and has to be provided
  !  by the interface module radar_interface.f90. It does not appear explicitly in the present
  !  module belwo, because it is input via parameter list into the below procedures.
  !
  ! Riming induced melting at temperatures below Tmin can be simulated by accounting for
  ! temperature Tmeltbegin (K).
  ! As far Tmeltbegin < Tmin, the degree of melting is a linear function of temperature in the interval
  ! Tmeltbegin < T < T_freeze. Its assumed to be independent of diameter and riming rate.
  ! At T = Tmin the degree of melting gets the value meltdegTmin.
  ! At Tmeltbegin > Tmin nothing is done
  ! Please note, that in reality the effect will strongly depend on riming rate, which is not considered here.
  !
  !----------------------------------------------------------------------------------------!

  ! First a master function, later wrapper-functionen for single hydro meteors
  ! (Version with 2D-output-array fmelt(anzt, anzD))
  ! temperatur(anzt), Dref(anzt) and Tmax(anzt) will vary along the x dimension of the model.

  SUBROUTINE degree_of_melting_fun_vec(temperatur,anzt,Dvec,anz,Dref,Dexpo,Tmeltbegin,meltdegTmin,Tmin,Tmax,fmelt)

    IMPLICIT NONE

    INTEGER,       INTENT(in)  :: anzt, anz
    REAL(KIND=dp), INTENT(out) :: fmelt(anzt,anz)
    REAL(KIND=dp), INTENT(in)  :: Tmeltbegin, meltdegTmin, Dvec(anz), temperatur(anzt)

    REAL(KIND=dp), INTENT(in)  :: Dref(anzt), Dexpo
    REAL(KIND=dp), INTENT(in)  :: &
         Tmin, &        ! temperature, when smallest particle starts to melt
         Tmax(anzt)     ! temperature, when all particles are already melted
    REAL(KIND=dp) :: &
         Dfullmelted(anzt,anz), tmpmeltdegTmin, Dmat(anzt,anz), tempmat(anzt,anz), &
         Drefmat(anzt,anz), Tmaxmat(anzt,anz)
    INTEGER :: i

    fmelt = 0d0

    ! Assign dummy variable Dmat and tempmat in a way, which can be verctorized
    IF (anzt == 1) THEN
      DO i=1, anz
        Dmat(1,i) = Dvec(i)
        tempmat(1,i) = temperatur(1)
        Drefmat(1,i) = Dref(1)
        tmaxmat(1,i) = Tmax(1)
      END DO
    ELSE
      IF (anz == 1) THEN
        DO i=1, anzt
          Dmat(i,1) = Dvec(1)
          tempmat(i,1) = temperatur(i)
          Drefmat(i,1) = Dref(i)
          Tmaxmat(i,1) = Tmax(i)
        END DO
      ELSE
        DO i=1, anz
          Dmat(:,i) = Dvec(i)
        END DO
        DO i=1, anzt
          tempmat(i,:) = temperatur(i)
          Drefmat(i,:) = Dref(i)
          Tmaxmat(i,:) = Tmax(i)
        END DO
      END IF
    END IF

    ! Start of degree-of-melting-parameterization
    IF (Tmeltbegin < Tmin) THEN
      tmpmeltdegTmin = meltdegTmin
      WHERE (tempmat > Tmeltbegin)
        fmelt = &
             MIN(MAX(meltdegTmin/(Tmin-Tmeltbegin)*(tempmat-Tmeltbegin), 0d0), 1d0)
      ELSEWHERE
        fmelt = 0d0
      END WHERE
    ELSE
      tmpmeltdegTmin = 0d0
    END IF

    Dfullmelted = 0d0
    WHERE (tempmat >= Tmaxmat)
      fmelt = 1d0
    END WHERE
    WHERE (tempmat < Tmaxmat .AND. tempmat >= Tmin .AND. Dmat >= quasi_zero)
      Dfullmelted = Drefmat**Dexpo * LOG((Tmaxmat-Tmin)/(Tmaxmat-tempmat))
      fmelt = MIN((1d0-tmpmeltdegTmin) * Dfullmelted / Dmat**Dexpo + tmpmeltdegTmin, 1d0)
    END WHERE
    WHERE (tempmat < Tmaxmat .AND. tempmat >= Tmin .AND. Dmat < quasi_zero)
      fmelt = 1d0
    END WHERE

    RETURN
  END SUBROUTINE degree_of_melting_fun_vec


  !=================================================================================
  ! Routine, which can only be applied at single grid points
  !=================================================================================

  SUBROUTINE degree_of_melting_fun(temperatur,Dvec,anz,Dref,Dexpo,Tmeltbegin,meltdegTmin,Tmin,Tmax,fmelt)

    IMPLICIT NONE

    INTEGER,       INTENT(in)  :: anz
    REAL(KIND=dp), intent(out) :: fmelt(anz)
    REAL(KIND=dp), INTENT(in)  :: Tmeltbegin, meltdegTmin, Dvec(anz), temperatur

    REAL(KIND=dp), INTENT(in)  :: Dref, Dexpo
    REAL(KIND=dp), INTENT(in)  :: &
         Tmin, &  ! temperature, when smallest particle starts to melt
         Tmax     ! temperature, when all particles are already melted
    REAL(KIND=dp) :: Dfullmelted, tmpmeltdegTmin

    fmelt = 0d0

    IF (Tmeltbegin < Tmin) THEN
      tmpmeltdegTmin = meltdegTmin
      IF (temperatur > Tmeltbegin) THEN
!!$        fmelt = &
!!$             MIN(MAX(meltdegTmin/(Tmin-Tmeltbegin)*(temperatur-Tmeltbegin), 0d0), 1d0)
        ! the max(...,1d-20) is necessary because of strange behaviour when this code is inlined ...
        fmelt = &
             MIN(meltdegTmin/MAX(Tmin-Tmeltbegin,quasi_zero)*(temperatur-Tmeltbegin), 1d0)
      ELSE
        fmelt = 0.0
      END IF
    ELSE
      tmpmeltdegTmin = 0.0
    END IF

    IF (temperatur >= Tmax) THEN
      fmelt = 1d0
    ELSE IF (temperatur < Tmax .AND. temperatur >= Tmin) THEN
!!$      Dfullmelted = Dref**Dexpo * LOG((Tmax-Tmin)/(Tmax-temperatur))
      ! the max(...,1d-20) is necessary because of strange behaviour when this code is inlined ...
      Dfullmelted = Dref**Dexpo * LOG((Tmax-Tmin)/MAX(Tmax-temperatur,quasi_zero))
      WHERE (Dvec**Dexpo >= Dfullmelted .AND. Dvec >= quasi_zero)
        fmelt = (1d0-tmpmeltdegTmin) * Dfullmelted / Dvec**Dexpo + tmpmeltdegTmin
      ELSEWHERE
        fmelt = 1d0
      END WHERE
    END IF

    RETURN
  END SUBROUTINE degree_of_melting_fun

  !=================================================================================
  ! Routine, which can only be applied at one single grid points and only for one D
  !=================================================================================

  SUBROUTINE degree_of_melting_fun_single(temperatur,D,Dref,Dexpo,Tmeltbegin,meltdegTmin,Tmin,Tmax,fmelt)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(out) :: fmelt
    REAL(KIND=dp), INTENT(in)  :: Tmeltbegin, meltdegTmin, D, temperatur

    REAL(KIND=dp), INTENT(in)  :: Dref, Dexpo
    REAL(KIND=dp), INTENT(in)  :: &
         Tmin, &  ! temperature, when smallest particle starts to melt
         Tmax     ! temperature, when all particles are already melted
    REAL(KIND=dp) :: Dfullmelted, tmpmeltdegTmin, D_Dexpo

    fmelt = 0d0

    IF (Tmeltbegin < Tmin) THEN
      tmpmeltdegTmin = meltdegTmin
      IF (temperatur > Tmeltbegin) THEN
!!$        fmelt = &
!!$             MIN(MAX(meltdegTmin/(Tmin-Tmeltbegin)*(temperatur-Tmeltbegin), 0d0), 1d0)
        ! the max(...,1d-20) is necessary because of strange behaviour on the NEC when this code is inlined ...
        fmelt = &
             MIN(meltdegTmin/MAX(Tmin-Tmeltbegin,quasi_zero)*(temperatur-Tmeltbegin), 1d0)
      ELSE
        fmelt = 0.0
      END IF
    ELSE
      tmpmeltdegTmin = 0.0
    END IF

    IF (temperatur >= Tmax) THEN
      fmelt = 1d0
    ELSE IF (temperatur < Tmax .AND. temperatur >= Tmin) THEN
!!$      Dfullmelted = Dref**Dexpo * LOG((Tmax-Tmin)/(Tmax-temperatur))
      ! the max(...,1d-20) is necessary because of strange behaviour when this code is inlined ...
      Dfullmelted = EXP(LOG(Dref)*Dexpo) * LOG((Tmax-Tmin)/MAX(Tmax-temperatur,quasi_zero))
      IF (D >= quasi_zero) THEN
        D_Dexpo = EXP(LOG(D)*Dexpo)
        fmelt = MERGE((1d0-tmpmeltdegTmin) * Dfullmelted / D_Dexpo + tmpmeltdegTmin, 1d0, &
                      D_Dexpo >= Dfullmelted)
      ELSE
        fmelt = 1d0
      END IF
    END IF

    RETURN
  END SUBROUTINE degree_of_melting_fun_single


  !=============================================================================================
  ! Routines, which can be applied only at a single grid point
  !=============================================================================================

  ! Generalized (1- & 2-mom) version of a common routine:
  SUBROUTINE degree_of_melting_xxx(mgd,T_a,Dvec,anz,Dexpo,Dref_max, &
                                   parti,Tmeltbegin,meltdegTmin,Tmin,Tmax,fmelt)

    IMPLICIT NONE

    INTEGER,            INTENT(in)  :: anz
    REAL(KIND=dp),      INTENT(in)  :: Tmax, Tmin, Tmeltbegin, meltdegTmin, Dvec(anz), &
                                       T_a, Dexpo, Dref_max
    TYPE(t_mgd_params), INTENT(in)  :: mgd
    CLASS(particle),    INTENT(in)  :: parti

    REAL(KIND=dp),      INTENT(out) :: fmelt(anz)

    REAL(KIND=dp) :: xref, Dref

    ! Dref is the diameter of a particle that is just melted completely at T of
    ! 1/e of (T-Tmin)/(Tmax-Tmin).
    ! It is set here to the medial diameter of the size spectrum.

    xref = MIN(MAX(mgd%q/(mgd%qn+quasi_zero),parti%x_min),parti%x_max)
    ! In case we don't want to depend on x_min settings (for 1mom, they are not
    ! from native model, but set from 2mom values which are not consistent with
    ! Dref(D_of_X_average_1M_exp) low-q bounds
    !xref = MIN(MAX(mgd%q/(mgd%qn+quasi_zero),quasi_zero),parti%x_max)
    ! most accurate, but probably not numerically stable
    !xref = mgd%q/(mgd%qn+quasi_zero)

    Dref = MIN( (xref / parti%a_geo)**(1d0/parti%b_geo), Dref_max)

    CALL degree_of_melting_fun(T_a,Dvec,anz,Dref,Dexpo,Tmeltbegin,meltdegTmin, &
                               Tmin,Tmax,fmelt)

    RETURN
  END SUBROUTINE degree_of_melting_xxx

  ! Generalized (1- & 2-mom) version of a common routine with fixed Dref:
  SUBROUTINE degree_of_melting_xxx_Dref(T_a,Dvec,anz,Dexpo,Dref    , &
                                        Tmeltbegin,meltdegTmin,Tmin,Tmax,fmelt)

    IMPLICIT NONE

    INTEGER,       INTENT(in)  :: anz
    REAL(KIND=dp), INTENT(in)  :: Tmax, Tmin, Tmeltbegin, meltdegTmin, Dvec(anz), &
                                  T_a, Dexpo, Dref

    REAL(KIND=dp), INTENT(out) :: fmelt(anz)

    ! Dref is the diameter of a particle that is just melted completely at T of
    ! 1/e of (T-Tmin)/(Tmax-Tmin).
    ! It is set as a constant in the interface, therefore fmelt does not depend
    !  on the PSD. For example, this makes it possible to speed up the generation
    !  of the lookup tables, because for a given set of T_a and Tmax, the
    !  scattering amplitudes do not depend on the PSD parameters and have
    !  to be computed only once for all D_s and can be recycled for any PSD parameters.

    CALL degree_of_melting_fun(T_a,Dvec,anz,Dref,Dexpo,Tmeltbegin,meltdegTmin, &
                               Tmin,Tmax,fmelt)

    RETURN
  END SUBROUTINE degree_of_melting_xxx_Dref

  ! Version only for one diameter:
  SUBROUTINE degree_of_melting_xxx_single(xref,T_a,D,Dexpo,Dref_max, &
                                          parti,Tmeltbegin,meltdegTmin,Tmin,Tmax,fmelt)

    IMPLICIT NONE

    REAL(KIND=dp),   INTENT(in)  :: xref, Tmax, Tmin, Tmeltbegin, meltdegTmin, D, &
                                    T_a, Dexpo, Dref_max
    CLASS(particle), INTENT(in)  :: parti

    REAL(KIND=dp),   INTENT(out) :: fmelt

    REAL(KIND=dp) :: Dref

    ! Dref is the diameter of a particle that is just melted completely at T of
    ! 1/e of (T-Tmin)/(Tmax-Tmin).
    ! It is set here to the medial diameter of the size spectrum.

    !xref = MIN(MAX(mgd%q/(mgd%qn+quasi_zero),parti%x_min),parti%x_max)
    ! In case we don't want to depend on x_min settings (for 1mom, they are not
    ! from native model, but set from 2mom values which are not consistent with
    ! Dref(D_of_X_average_1M_exp) low-q bounds
    !xref = MIN(MAX(mgd%q/(mgd%qn+quasi_zero),quasi_zero),parti%x_max)
    ! most accurate, but probably not numerically stable
    !xref = mgd%q/(mgd%qn+quasi_zero)

    Dref = MIN( (xref / parti%a_geo)**(1d0/parti%b_geo), Dref_max)

    CALL degree_of_melting_fun_single(T_a,D,Dref,Dexpo,Tmeltbegin,meltdegTmin, &
                                      Tmin,Tmax,fmelt)

    RETURN
  END SUBROUTINE degree_of_melting_xxx_single

  ! Version only for one diameter and fixed Dref:
  SUBROUTINE degree_of_melting_xxx_single_Dref(T_a,D,Dexpo,Dref, &
                                               Tmeltbegin,meltdegTmin,Tmin,Tmax,fmelt)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)  :: Tmax, Tmin, Tmeltbegin, meltdegTmin, D, &
                                  T_a, Dexpo, Dref

    REAL(KIND=dp), INTENT(out) :: fmelt

    ! Dref is the diameter of a particle that is just melted completely at T of
    ! 1/e of (T-Tmin)/(Tmax-Tmin).
    ! It is set as a constant in the interface, therefore fmelt does not depend
    !  on the PSD.

    CALL degree_of_melting_fun_single(T_a,D,Dref,Dexpo,Tmeltbegin,meltdegTmin, &
                                      Tmin,Tmax,fmelt)

    RETURN
  END SUBROUTINE degree_of_melting_xxx_single_Dref

END MODULE radar_mie_meltdegree
