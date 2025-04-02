! Source module for the radar forward operator EMVORADO
!
! ---------------------------------------------------------------
! Copyright (C) 2014-2024, DWD
! Contact information: ulrich.blahak (at) dwd.de 
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE radar_bubblegen

!------------------------------------------------------------------------------
!
! Description:
!   This module provides routines for the warm bubble generator which is
!   tighlty connected to the radar forward operator EMVORADO.
!
! Method:
!   See subroutines below
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

  USE radar_kind, ONLY : dp

  USE radar_data, ONLY :  &
       cmaxlen,           &
       my_radar_id_dom,   &
       radario_master_dom,&
       degrad,            & ! factor for transforming degree to rad
       composite_meta_type,    &
       bubble_list, nautobubbles_max, composite_meta_type, &
       pi

  USE radar_parallel_utilities, ONLY :  &
       global_values_radar

  USE radar_composites, ONLY : &
       comp_cellindex2geo

!==============================================================================

#ifndef NOMPI
  USE mpi
#endif

!==============================================================================

  IMPLICIT NONE

!==============================================================================

! default private
  PRIVATE

!==============================================================================

#ifdef NOMPI
  INCLUDE "nompi_mpif.h"
#endif

!==============================================================================

  ! define parameters for array dimensions for object detection
  INTEGER, PARAMETER :: NMAX_OBJECTS = 1000

!==============================================================================

PUBLIC detect_missing_cells

!==============================================================================
!==============================================================================

CONTAINS

!==============================================================================
!==============================================================================

  !============================================================================
  !
  ! Subroutines related to the warm bubble detector (ldo_bubbles=.true.)
  !
  !============================================================================


  !============================================================================
  ! Subroutine that checks if some point is in certain object. It includes two thresholds:
  !  the lower one determines IF the point is in the object. The second one
  ! counts how many points are above a second thresshold (with the ide of
  !  identifying high precipitation cores)
  !============================================================================

  SUBROUTINE check_pixel_inobject (inobject, p_not_checked, queue, ni,nj, &
       npoints_queue, object, npoints_obj, iloc,jloc,intensity, thresholds)

     IMPLICIT NONE

    !INPUT
     INTEGER,             INTENT(in)            :: ni,nj                               ! Size of the array
     INTEGER,             INTENT(in)            :: iloc,jloc           ! Coordinates of the point to check
     REAL(kind=dp),       INTENT(in)            :: intensity           ! Rw data (in this case, reflectivities of the checked point.
     REAL(kind=dp),       INTENT(in)            :: thresholds(2)       ! Thresholds that define the object

    !OUTPUT
     LOGICAL,             INTENT(out)           :: inobject            ! The point initiates a new object
     LOGICAL,             INTENT(inout)         :: p_not_checked         ! True if the point has not been checked and False when it has been checked

     INTEGER,             INTENT(inout)         :: queue(ni*nj,2)          ! Queue of points in this object that need to be checked. The second index is x-y
     INTEGER,             INTENT(inout)         :: npoints_queue             ! Number of points in the queue

     INTEGER,             INTENT(inout)         :: object(ni*nj,2)         ! Points that belong tp the current object. The second index is x-y
     INTEGER,             INTENT(inout)         :: npoints_obj(2)         ! Number of points above certain threshold. The first number is also the number of points in the object.


     p_not_checked = .FALSE.                                ! The point is checked

     IF ( intensity > thresholds(1) ) THEN
       inobject = .TRUE.


       ! Check the second threshol and add the counter
       IF ( intensity .GT. thresholds(2) ) THEN
         npoints_obj(2) = npoints_obj(2) + 1
       END IF

       ! Point added to the queue
       npoints_queue = npoints_queue +1
       queue(npoints_queue,1) = iloc
       queue(npoints_queue,2) = jloc

       ! Point added to the object
       npoints_obj(1) = npoints_obj(1) +1
       object(npoints_obj(1),1) = iloc
       object(npoints_obj(1),2) = jloc

     END IF

   END SUBROUTINE check_pixel_inobject

!==============================================================================
!==============================================================================

   ! Find connected objects. Find connected objects in comp_array defined by certain threshold.
   ! The object have to have nmin_thres(1) number of point above the minimum threshold and nmin_thres(2) above a second threshold
   ! The outpout is one array with the coordinates of i-j position of the estimated center abd a second array with the maximum i-j length of the object.
   ! This routine can be imporved by using elipsoides or a splitting algorithm for large objects (this is already developed in NCL)
   SUBROUTINE find_connected_objects( ellipse_arr,nobjects,comp_array, not_checked, &
        queue, object,ni,nj,thresholds,nmin,mult_dist,add_dist)

     IMPLICIT NONE

     !INPUT
     INTEGER,             INTENT(in)      :: ni, nj                           ! Size of the array
     REAL(kind=dp),       INTENT(in)      :: comp_array(ni,nj)                ! Array with the raw data (in this case, composite reflectivities)
     INTEGER ,            INTENT(inout)   :: queue(ni*nj,2), object(ni*nj,2)  ! Object and queue memory space (allocated previously)
     LOGICAL,             INTENT(inout)   :: not_checked(ni+2,nj+2)
     REAL(kind=dp),       INTENT(in)      :: thresholds(2)                    ! Thresholds for the connecting algorithm
     INTEGER,             INTENT(in)      :: nmin(2)                          ! Mninimum size criteria for the connecting algorithm
     REAL (kind=dp),      INTENT(in)      :: mult_dist                        !  Multiplicative axis factor for the object. If it is high, bubbles are more difficult
     REAL (kind=dp),      INTENT(in)      :: add_dist                         !  Additive axis increase for the object. If it is high, bubbles are more difficult

     !OUTPUT
     REAL (kind=dp),      INTENT(out)              :: ellipse_arr(NMAX_OBJECTS,10) ! Array with properties of the objects
     INTEGER                                       :: nobjects                     ! Number of objects

     ! LOCAL VARIABLES
     INTEGER                                       :: np_obj(2), np_queue          ! Number of points in current objects and queue
     INTEGER                                       :: ii,jj                        ! Counters
     INTEGER                                       :: iloc,jloc                    ! Local coordianates from queue
     LOGICAL                                       :: in_object                    ! Tells if point is in object

     !Initiate counters and logical arrays
     np_obj(1:2) = 0
     np_queue = 0
     nobjects = 0
     in_object = .FALSE.
     !Initiate not checked array. Set borders as already checked.
     not_checked(: ,: ) = .TRUE.
     not_checked(1 ,: ) = .FALSE.
     not_checked(ni+2,: ) = .FALSE.
     not_checked(: ,1 ) = .FALSE.
     not_checked(: ,nj+2) = .FALSE.

     DO ii = 1,ni
       DO jj = 1,nj
         !Check first point
         IF  ( not_checked(ii+1,jj+1) )  CALL check_pixel_inobject (in_object, not_checked(ii+1,jj+1), &
              queue, ni,nj, np_queue, object, np_obj, ii,jj,comp_array(ii,jj), thresholds)

         IF (in_object) THEN ! A new object is started
           DO WHILE ( np_queue .GT. 0)
             iloc = queue(np_queue,1)
             jloc = queue(np_queue,2)
             np_queue = np_queue -1

             ! Check all neighbouhs (different neighbours for different connectivities).
             ! These functions send a non-allocated address (but not used) for the borders.
             ! It might produce segmentation fault.
             !Check North
             IF  ( not_checked(iloc+2,jloc+1) )  CALL  check_pixel_inobject (in_object, not_checked(iloc+2,jloc+1), &
                  queue, ni,nj, np_queue, object, np_obj, iloc+1,jloc  ,comp_array(iloc+1,jloc  ), thresholds)
             !Check South
              IF  ( not_checked(iloc  ,jloc+1) ) CALL  check_pixel_inobject (in_object, not_checked(iloc  ,jloc+1), &
                   queue, ni,nj, np_queue, object, np_obj, iloc-1,jloc  ,comp_array(iloc-1,jloc  ), thresholds)
             !Check East
              IF  ( not_checked(iloc+1,jloc+2) ) CALL  check_pixel_inobject (in_object, not_checked(iloc+1,jloc+2), &
                   queue, ni,nj, np_queue, object, np_obj, iloc  ,jloc+1,comp_array(iloc  ,jloc+1), thresholds)
             !Check West
              IF  ( not_checked(iloc+1,jloc  ) ) CALL  check_pixel_inobject (in_object, not_checked(iloc+1,jloc  ), &
                   queue, ni,nj, np_queue, object, np_obj, iloc  ,jloc-1,comp_array(iloc  ,jloc-1), thresholds)
             !Check North-East
              IF  ( not_checked(iloc+2,jloc+2) ) CALL  check_pixel_inobject (in_object, not_checked(iloc+2,jloc+2), &
                   queue, ni,nj, np_queue, object, np_obj, iloc+1,jloc+1,comp_array(iloc+1,jloc+1), thresholds)
             !Check North-West
              IF  ( not_checked(iloc+2,jloc  ) ) CALL  check_pixel_inobject (in_object, not_checked(iloc+2,jloc  ), &
                   queue, ni,nj, np_queue, object, np_obj, iloc+1,jloc-1,comp_array(iloc+1,jloc-1), thresholds)
             !Check South-East
              IF  ( not_checked(iloc  ,jloc+2) ) CALL  check_pixel_inobject (in_object, not_checked(iloc  ,jloc+2), &
                   queue, ni,nj, np_queue, object, np_obj, iloc-1,jloc+1,comp_array(iloc-1,jloc+1), thresholds)
             !Check South-West
              IF  ( not_checked(iloc  ,jloc  ) ) CALL  check_pixel_inobject (in_object, not_checked(iloc  ,jloc  ), &
                   queue, ni,nj, np_queue, object, np_obj, iloc-1,jloc-1,comp_array(iloc-1,jloc-1), thresholds)
           END DO


           IF ( (np_obj(1) .GT. nmin(1) ) .AND. ( np_obj(2) .GT. nmin(2) ) ) THEN ! Only consider large objects with nmin(2) large reflectivity values
             nobjects = nobjects + 1

             ! Find the properties of elliptic approximation of the object:
             CALL  object_ellipse( ellipse_arr(nobjects,:), object(1:np_obj(1),:), np_obj(1), mult_dist, add_dist )

           END IF
           !Reinitialitate counters for new object
           np_obj(1:2) = 0
           np_queue = 0
           in_object = .FALSE.

         END IF ! End object
       END DO
     END DO

   END SUBROUTINE find_connected_objects

!==============================================================================
!==============================================================================

   SUBROUTINE detect_missing_cells(idom, time_mod_sec, composite_obs, composite_mod, comp_meta, &
        dt_bubble_search, prob_bubble, maxdim_obs, lbub_isolated, &
        threshold_obs, threshold_mod, areamin_mod, areamin_obs, &
        mult_dist_obs, add_dist_obs, mult_dist_mod, add_dist_mod, ldebug)

     IMPLICIT NONE

     !INPUT
     !-------------------

     INTEGER, INTENT(in)          :: idom         ! domain identifier in the hosting model [1-ndoms]
     REAL (kind=dp), INTENT(in)   :: time_mod_sec ! model time in seconds since model start

     REAL (kind=dp), INTENT(in)   :: composite_obs(:,:), composite_mod(:,:) ! Composite fields from observations and radar operator
     TYPE(composite_meta_type), INTENT(in) :: comp_meta

     ! Parameters for detection and downstream advection of automatic bubbles
     REAL (kind=dp), INTENT(in)   :: dt_bubble_search  ! Time interval from one automatic bubble search to the next [seconds] (SYNCHRONIZE WITH COMPOSITE TIMES!!!)
     REAL (kind=dp), INTENT(in)   :: prob_bubble       ! Probability of triggering a bubble when it is found  [0-1]
     REAL (kind=dp), INTENT(in)   :: maxdim_obs        ! Maximum dimension of not found cell that can trigger a bubble (very large objects are not targeted) [meters]
     LOGICAL, INTENT(in)          :: lbub_isolated        ! Check that the bubbles are isolated from other objects in observations
     REAL (kind=dp), INTENT(in)   :: threshold_obs(2), threshold_mod(2)  ! Thresholds that define the minimum dBZ in an object- and the high intensity region [dBZ]
     REAL (kind=dp), INTENT(in)   :: areamin_mod(2), areamin_obs(2)     ! Minimum area of the object- and high intensity region [m^2]
     REAL (kind=dp), INTENT(in)   :: mult_dist_obs     ! Multiplicative axis factor for the obs object. If it is high, bubbles are more difficult [-]
     REAL (kind=dp), INTENT(in)   :: mult_dist_mod     ! Multiplicative axis factor for the sim object. If it is high, bubbles are more difficult [-]
     REAL (kind=dp), INTENT(in)   :: add_dist_obs      ! Additive axis increase for the obs object. If it is high, bubbles are more difficult [meters]
     REAL (kind=dp), INTENT(in)   :: add_dist_mod      ! Additive axis increase for the sim object. If it is high, bubbles are more difficult [meters]
     LOGICAL, INTENT(in)          :: ldebug            ! Print debug messages

     ! LOCAL Variables
     !-------------------

     ! Output of finding algorithm
     REAL (kind=dp)               :: ellipse_obs(NMAX_OBJECTS,10),  & ! Ellipse found by the object algorithm
                                     ellipse_mod(NMAX_OBJECTS,10)
     INTEGER                      :: nobj_obs, nobj_mod
     ! Memory arrays for finding algorithm
     LOGICAL, DIMENSION(:,:), ALLOCATABLE :: not_checked
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: queue(:,:), object(:,:)
     INTEGER                              :: ni, nj

     ! Warm bubbles data
     REAL(KIND=dp)             :: i_newbubble_comp, j_newbubble_comp ! i and j position in composite coordinates of the new warm bubble (PROBABLY TOO LARGE ARRAYS)
     REAL(KIND=dp)             :: i_newbubble_mod, j_newbubble_mod   ! i and j position in model coordinates of the new warm bubble (PROBABLY TOO LARGE ARRAYS)
     REAL(KIND=dp)             :: x_newbubble_mod, y_newbubble_mod   ! x and y position in model coordinates of the new warm bubble (PROBABLY TOO LARGE ARRAYS)
     INTEGER                   :: nbubbles, nbubbles_reject
     LOGICAL                   :: is_inside_domain

     ! Internal variables
     INTEGER                   :: nmin_mod(2), nmin_obs(2) ! Minimum area of the object- and high intensity region [pixels]
     INTEGER                   :: maxsize_obs       ! Maximum dimension of not founded cell that can trigger a bubble (very large objects are not targeted) [pixels]
     REAL (kind=dp)            :: add_pix_obs       ! Additive axis increase for the obs object. If it is high, bubbles are more difficult [pixels]
     REAL (kind=dp)            :: add_pix_mod       ! Additive axis increase for the sim object. If it is high, bubbles are more difficult [pixels]
     INTEGER                   :: ccmod,ccobs,ccobs2,i,cc  ! Counters
     REAL (kind=dp)            :: dcenter           ! Distance beteween objects
     INTEGER                   :: ierr              ! Error from MPI
     REAL (kind=dp)            :: zdx,zdy                ! Grid spacing (in m)
     LOGICAL                   :: a_isolated(NMAX_OBJECTS) ! Check if objects are isolated
     REAL (KIND=dp)            :: lon_geo, lat_geo       ! geogr. positions of bubble centers

     CHARACTER(LEN=50)         :: yzroutine        ! Error Handling
     CHARACTER(LEN=cmaxlen)    :: yerrmsg
     INTEGER                   :: ierror

     yzroutine(:) = ' '
     yzroutine    = 'detect_missing_cells'
     yerrmsg(:)   = ' '

     ! grid lengths (without metrical terms):
     zdx   = comp_meta%r_earth * comp_meta%dlon * degrad
     zdy   = comp_meta%r_earth * comp_meta%dlat * degrad

     ! Initialize bubble positions and times with neutral values on all PEs:
     bubble_list%nbubbles         = 0
     bubble_list%nbubbles_reject  = 0
     bubble_list%bub_centlon(:)   = -HUGE(1.0_dp)
     bubble_list%bub_centlat(:)   = -HUGE(1.0_dp)
     bubble_list%bub_timestamp(:) = -HUGE(1.0_dp)

     maxsize_obs = CEILING(maxdim_obs/SQRT(zdx*zdy)) + 1
     add_pix_obs = add_dist_obs / SQRT(zdx*zdy)
     add_pix_mod = add_dist_mod / SQRT(zdx*zdy)
     nmin_mod(1:2) = CEILING (areamin_mod / (zdx*zdy))
     nmin_obs(1:2) = CEILING (areamin_obs / (zdx*zdy))

     ! .. bubble search on proc radario_master_dom in scope icomm_radar_dom (=0 if synchr. radar IO, >= num_compute_fwo if asynchr. IO):
     IF (my_radar_id_dom(idom) == radario_master_dom(idom)) THEN

       ! Allocate memory for arrays
       ni = SIZE(composite_obs, dim=1)
       nj = SIZE(composite_obs, dim=2)
       ALLOCATE(not_checked(ni+2,nj+2))
       ALLOCATE(queue(ni*nj,2))
       ALLOCATE(object(ni*nj,2))

       ! Find objects in observation composite
       CALL  find_connected_objects( ellipse_obs,nobj_obs,composite_obs, &
            not_checked, queue, object,ni,nj,threshold_obs,nmin_obs,mult_dist_obs,add_pix_obs)
       ! Find objects in model composite
       CALL  find_connected_objects( ellipse_mod,nobj_mod,composite_mod, &
            not_checked, queue, object,ni,nj,threshold_mod,nmin_mod,mult_dist_mod,add_pix_mod)

       DEALLOCATE(not_checked)
       DEALLOCATE(queue)
       DEALLOCATE(object)

       ! Print information on identified objects:
       IF (ldebug) THEN
         WRITE(*,'(a,i4,a,i4,a)') 'INFO '//TRIM(yzroutine)//': ', &
              nobj_obs,' objects in observations and ', nobj_mod, ' objects in model:'
         DO ccobs = 1, nobj_obs
           i_newbubble_comp = ellipse_obs(ccobs,1)
           j_newbubble_comp = ellipse_obs(ccobs,2)
           CALL comp_cellindex2geo(i_newbubble_comp, j_newbubble_comp, comp_meta, lon_geo, lat_geo)
           WRITE (*,'(a,i3,2(a,f0.5),2(a,f0.1),a,f0.1,a)') '  \_ in obs Nr. ', ccobs, &
                ' lon = ', lon_geo, ' lat = ', lat_geo, &
                ' a = ', ellipse_obs(ccobs,3), ' b = ', ellipse_obs(ccobs,4), &
                ' time = ', time_mod_sec, ' s'
         END DO
         DO ccmod = 1, nobj_mod
           i_newbubble_comp = ellipse_mod(ccmod,1)
           j_newbubble_comp = ellipse_mod(ccmod,2)
           CALL comp_cellindex2geo(i_newbubble_comp, j_newbubble_comp, comp_meta, lon_geo, lat_geo)
           WRITE (*,'(a,i3,2(a,f0.5),2(a,f0.1),a,f0.1,a)') '  \_ in model Nr. ', ccmod, &
                ' lon = ', lon_geo, ' lat = ', lat_geo, &
                ' a = ', ellipse_mod(ccmod,3), ' b = ', ellipse_mod(ccmod,4), &
                ' time = ', time_mod_sec, ' s'
         END DO
       END IF

       bubble_list%nbubbles = 0
       bubble_list%nbubbles_reject = 0
       a_isolated(:) = .TRUE.
       !Bubble candidates are those that objects in obs space, which are far from objects in model space
       DO ccobs = 1, nobj_obs ! Check every observation object

         ! Check if the bubble candidate is isolated from other candidates in observation space (this avoids disturbing large structures)
         IF(lbub_isolated) THEN
           DO ccobs2 = ccobs+1,nobj_obs
             IF (  overlapping_ellipses( ellipse_obs(ccobs,:),ellipse_obs(ccobs2,:) ) ) THEN !Second check If some object in model is too close, then it is not a candidate
               a_isolated(ccobs)  = .FALSE.
               a_isolated(ccobs2) = .FALSE.
               IF (ldebug) THEN
                 WRITE(*,'(a,i4,a,i4,a,f0.1,a)') 'INFO '//TRIM(yzroutine)//': Observed bubble candidates no. ', &
                      ccobs,' and ', ccobs2, ' removed due to not beeing isolated '// &
                      ' time = ', time_mod_sec, ' s'
               END IF
             END IF
           END DO
         END IF

         IF (a_isolated(ccobs)) THEN
           IF ( ellipse_obs(ccobs,3) .LT. maxsize_obs) THEN ! First check, only relative small objects in observation space are checked
             IF ( ellipse_obs(ccobs,1) .GT. 20 .AND.  ellipse_obs(ccobs,1) .LT. (ni-20) .AND. &
                  ellipse_obs(ccobs,2) .GT. 20 .AND.  ellipse_obs(ccobs,2) .LT. (nj-20)  ) THEN ! Check bubbles not in boundary
               ccmod = 0

               DO WHILE (ccmod .LT. nobj_mod)
                 ccmod = ccmod +1
                 IF (  overlapping_ellipses( ellipse_obs(ccobs,:),ellipse_mod(ccmod,:) ) ) THEN
                   ccmod = NMAX_OBJECTS + 10 ! Exit loop, this one cannot be a bubble
                 END IF
               END DO

               IF (ccmod .LT. ( NMAX_OBJECTS + 2) ) THEN ! Second check If it did go through the previous loop without triggering, then it is an object candidate
                 ! IF (Random) THEN ! Final Check Random bubble cna be done here for maybe little performance or at the end for better readability
                 !Finally, we have a bubble !!!

                 i_newbubble_comp = ellipse_obs(ccobs,1)
                 j_newbubble_comp = ellipse_obs(ccobs,2)

                 ! Convert bubble's composite coordinates into geographic coordinates:
                 CALL comp_cellindex2geo(i_newbubble_comp, j_newbubble_comp, comp_meta, lon_geo, lat_geo)

                 IF ( bubble_list%nbubbles < nautobubbles_max ) THEN

                   bubble_list%nbubbles = bubble_list%nbubbles + 1
                   bubble_list%bub_centlon(bubble_list%nbubbles)   =  lon_geo
                   bubble_list%bub_centlat(bubble_list%nbubbles)   =  lat_geo
                   bubble_list%bub_timestamp(bubble_list%nbubbles) =  time_mod_sec
                   IF (ldebug) THEN
                     WRITE (*,'(a,i3,2(a,f0.5),2(a,f0.1),a,f0.1,a)') &
                          'INFO '//TRIM(yzroutine)//': Accepted bubble candidate ', ccobs, &
                          ' at lon = ', lon_geo, ' lat = ', lat_geo, &
                          ' a = ', ellipse_obs(ccobs,3), ' b = ', ellipse_obs(ccobs,4), &
                          ' time = ', time_mod_sec, ' s'
                   END IF

                 ELSE

                   bubble_list%nbubbles_reject = bubble_list%nbubbles_reject + 1
                   IF (ldebug) THEN
                     WRITE(*,'(a,i3,a,f0.4,a,f0.4,a,f0.1,a)') &
                          'WARNING '//TRIM(yzroutine)//': bubble candidate ', ccobs, &
                          ' rejected at lon = ' , lon_geo , ' lat = ', lat_geo , &
                          ' time = ',time_mod_sec,' s because nautobubbles_max is too small'
                   END IF

                 END IF

               END IF

             END IF   ! Boundary check
           END IF  ! IF first check small objects
         END IF  ! Isolated check

       END DO

       WRITE(*,'(a,i3,a,f0.1,a)') 'INFO SUMMARY '//TRIM(yzroutine)//': ', &
            bubble_list%nbubbles, ' missing cells detected in model '// &
            ' at time = ', time_mod_sec, ' s'

       !Check that we do not have too many bubbles
       IF  (bubble_list%nbubbles_reject .GT. 0) THEN
         WRITE(*,'(a,f0.1,a,i4,a,i4,a)') &
              'WARNNING '//TRIM(yzroutine)//': Too many missing cells found at time = ', time_mod_sec,': ',&
              bubble_list%nbubbles+bubble_list%nbubbles_reject , &
              'were found, but only ', bubble_list%nbubbles, ' will be used. Increase nautobubbles_max!'
       END IF

     END IF ! If proccesor radario_master

   END SUBROUTINE detect_missing_cells

!==============================================================================
!==============================================================================

!  This function provides the best ellipse that fits a continous region provided by obj
!  The output is in a ellise in array: (consider an object for the future)
!                          index 1: x position of the center
!                          index 2: y position of the center
!                          index 3: length of larger axis
!                          index 4: length of shorter axis
!                          index 5: cos of angle of larger axis with x (same as x-coord of larger eigenvector)
!                          index 6: sin of angle of larger axis with x (same as x-coord of larger eigenvector)
!                          index 7: x position of focal point 1
!                          index 8: y position of focal point 1
!                          index 9: x position of focal point 2
!                          index 10: y position of focal point 3

! Provides an array of x,y points on the ellipse from the ellipse object described above

   SUBROUTINE xy_ellipse(xelip, yelip, Npoints, ellipse)

     IMPLICIT NONE

      ! INPUT
      REAL(kind=dp), INTENT(in)     ::   ellipse(:)
      INTEGER,       INTENT(in)     ::   Npoints

      ! OUTPUT
      REAL(kind=dp), INTENT (out)   ::   xelip(Npoints),yelip(Npoints)

      !LOCAL VARIABLES
      REAL (kind=dp)                ::   vec_large(2), vec_short(2) ! large and short vector of the ellipse
      INTEGER                       ::   ii                         ! counter
      REAL (kind=dp)                ::   theta                      ! angle of parametric representation of the ellipse

      vec_large(1) =  ellipse(3)*ellipse(5)
      vec_large(2) =  ellipse(3)*ellipse(6)
      vec_short(1) = -ellipse(4)*ellipse(6)
      vec_short(2) =  ellipse(4)*ellipse(5)

      DO ii = 1,Npoints
        theta = 2.0_dp*pi*real(ii-1)/real(Npoints)
        xelip(ii) = ellipse(1) + vec_large(1)*cos(theta) + vec_short(1)*sin(theta)
        yelip(ii) = ellipse(2) + vec_large(2)*cos(theta) + vec_short(2)*sin(theta)
      END DO
   END SUBROUTINE xy_ellipse

   ! Detects if two ellipses do overlap (only aproximated)
   FUNCTION overlapping_ellipses( ellipse1, ellipse2 ) RESULT (overlap)

     IMPLICIT NONE

     ! INPUT
     REAL(kind=dp), INTENT(in)    ::   ellipse1(:), ellipse2(:)

     ! Output
     LOGICAL                      ::   overlap ! True if they overlap, False if they do not

     ! Local variables
     REAL(kind=dp)                ::   distance_center    ! Distance from center to center of both ellipses
     INTEGER                      ::   Npoints_check      ! Points in the smaller ellipse to check
     REAL (kind=dp), ALLOCATABLE  ::   xelip(:), yelip(:) ! Points in the ellipse
     INTEGER                      ::   ii                ! Counter

     distance_center = SQRT ( (ellipse1(1) - ellipse2(1) ) * (ellipse1(1) - ellipse2(1) )  +  &
                              (ellipse1(2) - ellipse2(2) ) * (ellipse1(2) - ellipse2(2) ) )

     IF ( distance_center .GT. ( ellipse1(3) + ellipse2(3)) ) THEN        !.. Easy case ellipses too far away too overlap
       overlap = .FALSE.
     ELSE IF ( distance_center .LT. ( ellipse1(4) + ellipse2(4)) ) THEN   !.. Easy case ellipses too far away too overlap
       overlap = .TRUE.
     ELSE                                                                 !.. Difficult case, intermediate state
       IF ( ellipse1(3) .GT. ellipse2(3) ) THEN     ! 2 is the smaller ellipse. Check if some points are inside the first ellipse
         Npoints_check = FLOOR(ellipse2(3)/8.0)*4    ! Check a number of points given by the large axis. Must be a multiple of 4.
         IF (Npoints_check .LT. 8) THEN              ! Smaller number of points
           Npoints_check = 8
         END IF
         ALLOCATE ( xelip( Npoints_check) )
         ALLOCATE ( yelip( Npoints_check) )
         CALL  xy_ellipse(xelip, yelip, Npoints_check, ellipse2) ! Get Npoints_check of ellipse 2

         ! Check if the points are inside the larger ellipse using the focal points
         overlap = .FALSE.
         DO ii=1,Npoints_check
           IF ( ( SQRT( (xelip(ii) - ellipse1(7 ) ) * (xelip(ii) - ellipse1(7 ) ) + &
                        (yelip(ii) - ellipse1(8 ) ) * (yelip(ii) - ellipse1(8 ) ) )  &
                + SQRT( (xelip(ii) - ellipse1(9 ) ) * (xelip(ii) - ellipse1(9 ) ) + &
                        (yelip(ii) - ellipse1(10) ) * (yelip(ii) - ellipse1(10) ) ) ) < 2.0*ellipse1(3) ) THEN
             overlap = .TRUE.
             EXIT ! Exits the do loop
           END IF
         END DO
       ELSE                                          ! 1 is the smaller ellipse. Check if some points are inside the first ellipse
         Npoints_check = FLOOR(ellipse1(3)/8.0_dp)*4    ! Check a number of points given by the large axis. Must be a multiple of 4.
         IF (Npoints_check .LT. 8) THEN              ! Smaller number of points
           Npoints_check = 8
         END IF
         ALLOCATE ( xelip( Npoints_check) )
         ALLOCATE ( yelip( Npoints_check) )
         CALL  xy_ellipse(xelip, yelip, Npoints_check, ellipse1) ! Get Npoints_check of ellipse 2

         ! Check if the points are inside the larger ellipse using the focal points
         overlap = .FALSE.
         DO ii=1,Npoints_check
           IF ( ( SQRT( (xelip(ii) - ellipse2(7 ) ) * (xelip(ii) - ellipse2(7 ) ) + &
                        (yelip(ii) - ellipse2(8 ) ) * (yelip(ii) - ellipse2(8 ) ) )  &
                + SQRT( (xelip(ii) - ellipse2(9 ) ) * (xelip(ii) - ellipse2(9 ) ) + &
                        (yelip(ii) - ellipse2(10) ) * (yelip(ii) - ellipse2(10) ) ) ) < 2.0*ellipse2(3) ) THEN
             overlap = .TRUE.
             EXIT ! Exits the do loop
           END IF
         END DO
       END IF
     END IF

   END FUNCTION overlapping_ellipses

!==============================================================================
!==============================================================================

   ! This function provides the best ellipse that fits a continous region provided by obj
   SUBROUTINE  object_ellipse(ellipse, obj, np_obj, mult_dist, add_dist)

     IMPLICIT NONE

     ! INPUT
     INTEGER,        INTENT(in)    :: np_obj          ! Number of points of object
     INTEGER,        INTENT(inout) :: obj(np_obj,2)   ! Pixel coordinates in object. The memory from this array is
                                                      !   used in the function, as the object is not longer needed.
     REAL (kind=dp), INTENT(in)    :: mult_dist       ! Multiplicative axis factor for the object. If it is high, bubbles are more difficult
     REAL (kind=dp), INTENT(in)    :: add_dist        ! Additive axis increase for the object. If it is high, bubbles are more difficult
     ! OUTPUT
     REAL (kind=dp), INTENT(out)   :: ellipse(10)     ! Best ellipse that covers the object

     ! LOCAL VARIABLES
     REAL (kind=dp)  ::  cov_matr(2,2), eigen_val(2), eigen_vec(2,2)
     REAL (kind=dp)  ::   np_obj_real, trace, deter, dummy, c_focal

     CHARACTER(len=*), PARAMETER :: routine='radar_bubblegen.f90::object_ellipse'

     np_obj_real = REAL( np_obj,dp )
     ! The center is the mean value of the object
     ellipse(1) = REAL( SUM( obj(:,1) ),dp) /  np_obj_real
     ellipse(2) = REAL( SUM( obj(:,2) ),dp) /  np_obj_real

     ! Center the object
     obj(:,1) =  obj(:,1) - NINT(ellipse(1))
     obj(:,2) =  obj(:,2) - NINT(ellipse(2))

     ! Calculate covariance matrix
     cov_matr(1,1) = DOT_PRODUCT( obj(:,1),obj(:,1))
     cov_matr(2,2) = DOT_PRODUCT( obj(:,2),obj(:,2))
     cov_matr(1,2) = DOT_PRODUCT( obj(:,1),obj(:,2))

     cov_matr = cov_matr/ np_obj_real

     ! Matrix eigenvalues
     trace = cov_matr(1,1)+ cov_matr(2,2)
     deter =  cov_matr(1,1)*cov_matr(2,2) -  cov_matr(1,2)* cov_matr(1,2)

     eigen_val(1) = trace/2.0_dp + SQRT(trace*trace/4.0_dp-deter)
     eigen_val(2) = trace/2.0_dp - SQRT(trace*trace/4.0_dp-deter)

     ! Matrix eigenvectors
     IF   ( ABS(cov_matr(1,2)) > 1d-10 ) THEN
       eigen_vec(1,1) =  eigen_val(1) - cov_matr(2,2)
       eigen_vec(2,1) =  eigen_val(2) - cov_matr(2,2)
       eigen_vec(:,2) =  cov_matr(1,2)

       eigen_vec(1,:) =  eigen_vec(1,:)/SQRT(  eigen_vec(1,1)* eigen_vec(1,1) +  eigen_vec(1,2)* eigen_vec(1,2) )
       eigen_vec(2,:) =  eigen_vec(2,:)/SQRT(  eigen_vec(2,1)* eigen_vec(2,1) +  eigen_vec(2,2)* eigen_vec(2,2) )
     ELSE
       eigen_vec = 0.0
       eigen_vec(1,1) = 1.0
       eigen_vec(2,2) = 1.0
     END IF

     ! Find axis length
     ellipse(3) = mult_dist*MAXVAL( ABS( eigen_vec(1,1)* obj(:,1) + eigen_vec(1,2)* obj(:,2) ) ) + add_dist
     ellipse(4) = mult_dist*MAXVAL( ABS( eigen_vec(2,1)* obj(:,1) + eigen_vec(2,2)* obj(:,2) ) ) + add_dist


     IF ( ellipse(4) > ellipse(3) ) THEN ! large and small axis are interchanged (rare case as the larger eigenvalue is almost always the larger axis)
       dummy = ellipse(3)
       ellipse(3) = ellipse(4)
       ellipse(4) = dummy
       dummy = eigen_vec(1,1)
       eigen_vec(1,1) = eigen_vec(2,1)
       eigen_vec(2,1) = dummy
       dummy = eigen_vec(1,2)
       eigen_vec(1,2) = eigen_vec(2,2)
       eigen_vec(2,2) = dummy
     END IF

     ! Set cosinus and sinus of larger axis
     ellipse(5) = eigen_vec(1,1)
     ellipse(6) = eigen_vec(1,2)

     ! Determine the focal points (needed by the overlapping algorithm)
     c_focal = ellipse(3)*ellipse(3) - ellipse(4)*ellipse(4)
     ! the MAX( ,0.0) is needed on the cray to prevent a floating point exception in rare cases:
     c_focal = SQRT( MAX(c_focal, 0.0_dp) )

     ellipse(7)  = ellipse(1) + c_focal*eigen_vec(1,1)
     ellipse(8)  = ellipse(2) + c_focal*eigen_vec(1,2)
     ellipse(9)  = ellipse(1) - c_focal*eigen_vec(1,1)
     ellipse(10) = ellipse(2) - c_focal*eigen_vec(1,2)

   END SUBROUTINE object_ellipse

!==============================================================================
!==============================================================================

END MODULE radar_bubblegen
