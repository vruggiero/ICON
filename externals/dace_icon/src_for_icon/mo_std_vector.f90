!
!+ utility routines for Slant Total Delay operator
!
! $Id$
!
Module mo_std_vector
!
! Description:
!   utility routines for Slant Total Delay operator
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_22        2013-02-13 Michael Bender
!  utility routines for Slant Total Delay operator
! V1_27        2013-11-08 Michael Bender
!  tl/adjoint routines
! V1_43        2015-08-19 Michael Bender
!  Raytracer added to STD operator; GNSS bias correction.
! V1_47        2016-06-06 Harald Anlauf
!  LinePos_ad: fix argument intent
! V1_50        2017-01-09 Michael Bender
!  extended namelist STD_OBS; modules merged with the COSMO code.
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================
!
!> @file mo_std_vector.f90
!> Routines used by the  GNSS STD observation operator
!> (@see mo_std_operator.f90  and @see mo_std_coord.f90)
!>
!> This module provides some routines to define a straight line as a directed
!> line segment (L = u + t*v), to create points on that line, to compute
!> intersection points between different lines or between a line and a sphere
!> or ellipsoid, ...
!>
!> List of module procedures:
!>
!> LineDef           Defines a directed line segment between two points
!> LineFind          Finds point with a given coordinate on a line
!> LinePos           Computes coordinates of point on line with lambda
!> LinePrintParam    Prints parameter defining a directed line segment
!> LinePrintPoint    Prints parameter of point on a directed line segment
!> LinePointDist     Computes the distance of a point from a straight line.
!> LineDist          Computes the distance between skew lines.
!> LineCrossHeight   Computes point with minimum distance to a second line
!> eAngle            Computes the angle between skew lines
!> CrossProduct      Computes the cross product of two vectors
!> QuadEq            Solves the quadratic equation
!> PointDist         Computes the distance between two points

!=====================================================================

#ifdef __COSMO__
use kind_parameters, only: wp
#else
use mo_kind, only: wp
#endif


implicit none

private
public :: line        ! straight line derived type in Cartesian coordinates
public :: LineDef     ! defines a directed line segment between two points
public :: QuadEq      ! solves the quadratic equation
public :: QuadEq_tl   ! tangent-linear of QuadEq
public :: QuadEq_ad   ! adjoint of QuadEq
public :: LinePos     ! Computes coordinates of point on line with lambda
public :: LinePos_tl  ! tangent-linear of LinePos
public :: LinePos_ad  ! adjoint of LinePos
public :: PointDist   ! distance between two points in cartesian coordinates


!---------------------------------------------------------------------
! directed line segment
!---------------------------------------------------------------------
!>
!> @brief Definition of a directed line segment
!>
!> A directed line segment between a start point and an end point can be
!> defined as L = u + t*v,
!> where u is a vector pointing from the origin of the reference frame
!> so some point on the line, v is an unit vector parallel to the line and
!> t is a free parameter which defines a point on that line.
!> It is not necessary to define a start or end point but it is more
!> convenient to compute a straight line through two given points than to
!> compute the vectors u and v.
!> The routine LineDef takes two points and computes all variables within
!> this deived type:
!> L = u + t*v = start + t*v,  v = unitvec = (end - start ) / len
!>
!> All vectors need to be defined in a Cartesian frame of reference!
!> <table>
!> <tr><th>variable  <th>description
!> <tr><td><b> dim </b>  <td>
!>    Dimension of the vectors: 2 and 3 dimensions are supported
!> <tr><td><b> start </b>  <td>
!>    Start point of the directed line segment in Cartesian coordinates
!> <tr><td><b> end </b>  <td>
!>    End point of the directed line segment in Cartesian coordinates
!> <tr><td><b> unitvec </b>  <td>
!>    Unit vector pointing from "start" to "end", i.e. (end - start ) / len
!> <tr><td><b> len </b>  <td>
!>    Length of the directed line segment, the unit of "len" depends on the
!>    units of "start" to "end".
!> </table>
!>
!---------------------------------------------------------------------
type Line
  integer                  :: dim     ! dimension (2 or 3)
  real(wp), dimension(1:3) :: start   ! start point of line
  real(wp), dimension(1:3) :: end     ! end point of line
  real(wp), dimension(1:3) :: unitvec ! unit vector from start to end
  real(wp)                 :: len     ! length of the line segment
end type Line


!---------------------------------------------------------------------
! point on a directed line segment
!---------------------------------------------------------------------
!>
!> @brief Definition of a point on a directed line segment
!>
!> <table>
!> <tr><th>variable  <th>description
!> <tr><td><b> point </b>  <td>
!>    Coordinates of a point on a line.
!> <tr><td><b> cindex </b>  <td>
!>    Cell indices, integer indices of a grid cell containing this point
!> <tr><td><b> length </b>  <td>
!>    Distance between the point and the start of the directed line segment.
!> </table>
!>
!---------------------------------------------------------------------
type LinePoint
   real(wp), dimension(1:3) :: point  ! Punkt auf der Geraden
   integer , dimension(1:3) :: cindex ! Gitterzelle, in der dieser Punkt liegt
   real(wp)                 :: length ! Abstand zum Startpunkt oder
                                      ! Wegstrecke in der Zelle
end type LinePoint

!> Some small real number used to decide if two reals are equal or different
real (wp), parameter :: epsilon = 1.0E-10_wp

!=====================================================================
contains
!=====================================================================

!---------------------------------------------------------------------
! subroutine LineDef
!---------------------------------------------------------------------
!
!> @brief Defines a directed line segment between two points.
!>
!> <b> call  LineDef(thisline, start, end, dim) </b>
!>
!> This routine takes two points and computes all variables within
!> the deived type "line":
!> L = u + t*v = start + t*v,  v = unitvec = (end - start ) / len
!>
!> @param[in,out] thisline  derived type holding all variables defining
!>                          the line
!> @param[in]     start     start point of the line segment
!> @param[in]     end       end point of the line segment
!> @param[in]     dim       dimension of points/vectors, dim = 2 or 3
!>
!---------------------------------------------------------------------
subroutine LineDef(thisline, start, end, dim)

implicit none

! List of calling arguments:
type (Line),              intent(inout) :: thisline
real(wp), dimension(1:3), intent(in)    :: start ! start of line
real(wp), dimension(1:3), intent(in)    :: end   ! end of line
integer,                  intent(in)    :: dim   ! dimension

! List of local variables:
real(wp), dimension(1:3) :: UnitVec      ! Vektor von Start zu Ende
!---------------------------------------------------------------------

! save parameter
thisline%start = start
thisline%end = end
thisline%dim = dim

! compute vector pointing from start to end
UnitVec = end - start

! distance between start to end
thisline%len = sqrt(UnitVec(1)**2 + UnitVec(2)**2 + UnitVec(3)**2)

! compute unit vector
thisline%unitvec = UnitVec / thisline%len

end subroutine LineDef


!---------------------------------------------------------------------
! subroutine LineFind
!---------------------------------------------------------------------
!
!> @brief Finds point with a given coordinate on a line
!>
!> <b> call LineFind(thisline, val, coord, point, lambda, online) </b>
!>
!> This routine tries to find a point on the line which has a given
!> coordinate x, y or z, e.g. if x0 is given a point (x0,y,z) on the line
!> is computed if such a point exists. Alternatively y0 or z0 can be provided.
!> The parameter lambda which is equivalent to the length of the line
!> segment between start and that point is also given.
!> ( p0 = (x0,y,z) = start + lambda_0 * UnitVec )
!>
!> @param[in]  thisline  derived type holding all variables defining the line
!> @param[in]  val       coordinate, x0 or y0 or z0
!> @param[in]  coord     dimension of "val": \n
!>                       (x0 - coord=1, y0 - coord=2, z0 - coord=3)
!> @param[out] point     coordinates of point on the line
!> @param[out] lambda    distance between "start" and "point"
!> @param[out] online    Could a point be found? \n
!>                       online = .true.  - a point on the line with "val"
!>                                          could be found \n
!>                       online = .false. - there is no such point. \n
!>                              This happens if the line is parallel to the
!>                              corresponding axis of the reference frame.
!>
!---------------------------------------------------------------------
subroutine LineFind(thisline, val, coord, point, lambda, online)

implicit none

! List of calling arguments:
type (Line),              intent(in)  :: thisline
real(wp),                 intent(in)  :: val
integer,                  intent(in)  :: coord
real(wp), dimension(1:3), intent(out) :: point
real(wp),                 intent(out) :: lambda
logical,                  intent(out) :: online

! List of local variables:
!---------------------------------------------------------------------

online = .true.
select case (coord)
case(1)
   if ( abs(thisline%unitvec(1)) > epsilon ) then
      lambda = (val - thisline%start(1))/thisline%unitvec(1)
      point(1) = val
      point(2) = thisline%start(2) + lambda*thisline%unitvec(2)
      point(3) = thisline%start(3) + lambda*thisline%unitvec(3)
   else
      online = .false.
   end if
case(2)
   if ( abs(thisline%unitvec(2)) > epsilon ) then
      lambda = (val - thisline%start(2))/thisline%unitvec(2)
      point(1) = thisline%start(1) + lambda*thisline%unitvec(1)
      point(2) = val
      point(3) = thisline%start(3) + lambda*thisline%unitvec(3)
   else
      online = .false.
   end if
case(3)
   if ( abs(thisline%unitvec(3)) > epsilon ) then
      lambda = (val - thisline%start(3))/thisline%unitvec(3)
      point(1) = thisline%start(1) + lambda*thisline%unitvec(1)
      point(2) = thisline%start(2) + lambda*thisline%unitvec(2)
      point(3) = val
   else
      online = .false.
   end if
case default
   write(*,*) 'Coordinate seems not to exist: ', coord
end select

end subroutine LineFind


!---------------------------------------------------------------------
! subroutine LinePos
!---------------------------------------------------------------------
!
!> @brief Computes coordinates of point on line with given parameter lambda
!>
!> <b> call LinePos(thisline, lambda, point) </b>
!>
!> For a given parameter lambda the coordinates of the corresponding point
!> on the line are computed: \n
!> lambda_0 is given: (x,y,z) = start + lambd_0 * UnitVec \n
!> Lambda is equivalent to the length of the line segment between
!> "start" and the point (x,y,z).
!>
!> @param[in]  thisline  derived type holding all variables defining the line
!> @param[in]  lambda    parameter lambda used to compute the coordinates
!> @param[out] point     point on the line with  parameter lambda
!>
!---------------------------------------------------------------------
subroutine LinePos(thisline, lambda, point)

implicit none

! List of calling arguments:
type (Line),              intent(in)  :: thisline
real(wp),                 intent(in)  :: lambda
real(wp), dimension(1:3), intent(out) :: point

! List of local variables:
integer :: i
!---------------------------------------------------------------------

do i=1, thisline%dim
   point(i) = thisline%start(i) + lambda*thisline%unitvec(i)
end do

end subroutine LinePos

!---------------------------------------------------------------------
! subroutine LinePrintParam
!---------------------------------------------------------------------
!
!> @brief Prints parameter defining a directed line segment to standard out.
!>
!> <b> call LinePrintParam( thisline ) </b>
!>
!> This routine can be used to print the current state of a directed line
!> segment to the screen. All variables of the derived type "line" are
!> printed.
!>
!> @param[in]  thisline  derived type holding all variables defining the line
!>
!---------------------------------------------------------------------
subroutine LinePrintParam( thisline )

implicit none

! List of calling arguments:
type (Line), intent(inout) :: thisline

! List of local variables:
!---------------------------------------------------------------------

write(*,'(a)')                        &
        '-------- directed line segment -----------------------------'
write(*,'(a,tr6,i1)')                 &
        'Dimension of the line  : ', thisline%dim
write(*,'(a,3(f18.6,tr4))')           &
        'Start  (x,y,z)         : ', thisline%start
write(*,'(a,3(f18.6,tr4))')           &
        'End    (x,y,z)         : ', thisline%end
write(*,'(a,3(f18.6,tr4))')           &
        'Unit vector (x,y,z)    : ', thisline%unitvec
write(*,'(a,(f18.6,tr4))')            &
        'Length of the segment  : ', thisline%len
write(*,'(a)')                        &
        '------------------------------------------------------------'

end subroutine LinePrintParam

!---------------------------------------------------------------------
! subroutine LinePrintPoint
!---------------------------------------------------------------------
!
!> @brief Prints parameter defining a point on a directed line segment
!>        to standard out.
!>
!> <b> call LinePrintPoint( lpoint ) </b>
!>
!> This routine can be used to print the current state of a point on a
!> directed line  segment to the screen.
!> All variables of the derived type "LinePoint" are
!> printed.
!>
!> @param[in]  lpoint  derived type holding all variables defining the point
!>
!---------------------------------------------------------------------
subroutine LinePrintPoint( lpoint )

implicit none

! List of calling arguments:
type (LinePoint), intent(in) :: lpoint

! List of local variables:
!---------------------------------------------------------------------

write(*,'(a,3(f7.3,tr2),a,3(i3,tr2),a,f5.3)')   &
        'Pos. (x,y,z)= ', lpoint%point,         &
        'Index (i,j,k)= ', lpoint%cindex,         &
        'Length= ', lpoint%length

end subroutine LinePrintPoint


!---------------------------------------------------------------------
! function LinePointDist
!---------------------------------------------------------------------
!
!> @brief Computes the distance of a point from a straight line.
!>
!> <b> dist = LinePointDist(g1, p) </b>
!>
!> For any given point the shortest disance to the line is computed.
!>
!> @param[in]  g1             derived type defining the line
!> @param[in]  p              point
!> @param[out] LinePointDist  shortest distance between p and g1
!>                            (return value)
!>
!---------------------------------------------------------------------
real(wp) function LinePointDist(g1, p)

implicit none

! List of calling arguments:
type(Line),               intent(in) :: g1
real(wp), dimension(1:3), intent(in) :: p

! List of local variables:
real(wp), dimension(1:3) :: c, d
!---------------------------------------------------------------------

c = p - g1%start
call CrossProduct(c, g1%unitvec, d)

LinePointDist = sqrt( DOT_PRODUCT(d, d) )

end function LinePointDist


!---------------------------------------------------------------------
! function LineDist
!---------------------------------------------------------------------
!
!> @brief Computes the distance between skew lines.
!>
!> <b> dist =  LineDist(g1, g2) </b>
!>
!> The minimum distance between skew lines is computed.
!>
!> @param[in]  g1             derived type defining the first line
!> @param[in]  g2             derived type defining the second line
!> @param[out] LineDist       shortest distance between g1 and g2
!>                            (return value)
!>
!---------------------------------------------------------------------
real(wp) function LineDist(g1, g2)

implicit none

! List of calling arguments:
type (Line), intent(in) :: g1, g2

! List of local variables:
real(wp), dimension(1:3) :: d, c
!---------------------------------------------------------------------

d = g1%start - g2%start
call CrossProduct(g1%unitvec, g2%unitvec, c)

!write(*,*) 'a = ', d
!write(*,*) 'b = ', c
!write(*,*) 'dotprod = ', DOT_PRODUCT(d, c)

LineDist = abs( DOT_PRODUCT(d, c) )

end function LineDist


!---------------------------------------------------------------------
! subroutine  LineCrossHeight
!---------------------------------------------------------------------
!
!> @brief Computes the coordinates of a point on a line which has a
!>        minimum distance to a second line.
!>
!> <b> call LineCrossHeight(g1, g2, dist, point, range, lambda) </b>
!>
!> The function "LineDist" provides the minimum distance between skew lines
!> but doesn't compute the corresponding points on the lines. This routine
!> computes the coordinates of one point on the first line which has the
!> minimum distance to a second line. To obtain the corresponding point on
!> the second line this routine has to be called again with the lines in
!> reverse order.\n
!> Example: \n
!> call LineCrossHeight(g1, g2, dist, point2) \n
!>      => Provides the coordinates of point2 on g2 which has the minimum
!>         distance to g1. \n
!> call LineCrossHeight(g2, g1, dist, point1) \n
!>      => Provides the coordinates of point1 on g1 which has the minimum
!>         distance to g2.
!>
!> @param[in]  g1   first line
!> @param[in]  g2   second line, point is located on this line
!> @param[in]  dist distance between skew lines g1 and g2
!> @param[out] point
!> @param[out] range   Is "point" between "start" and "end" of g2? \n
!>                     range=.true.  - point is between start and end of g2 \n
!>                     range=.false. - point is somewhere else \n
!>                     (optional)
!> @param[out] lambda  parameter defining the position of "point" on the
!>                     second line. (optional)
!>
!---------------------------------------------------------------------
subroutine LineCrossHeight(g1, g2, dist, point, range, lambda)

implicit none

! List of calling arguments:
type (Line),              intent(in)  :: g1, g2
real(wp),                 intent(in)  :: dist
real(wp), dimension(1:3), intent(out) :: point
logical,  optional,       intent(out) :: range
real(wp), optional ,      intent(out) :: lambda

! List of local variables:
real(wp), dimension(1:3) :: a, x, y
real(wp)                 :: p, q, t, t1, t2, ysquared
logical                  :: IsReal
!---------------------------------------------------------------------

a = g2%start - g1%start
call CrossProduct(a, g1%unitvec, x)
call CrossProduct(g2%unitvec, g1%unitvec, y)

!!$WRITE(*,*) 'a2 = ', g2%start
!!$WRITE(*,*) 'a1 = ', g1%start
!!$WRITE(*,*) 'a2-a1 = ', a
!!$WRITE(*,*) 'a2-a1 = ', g2%start(1)-g1%start(1),g2%start(2)-g1%start(2),  &
!!$                       g2%start(3)-g1%start(3)
!!$WRITE(*,*) 'm = ', x
!!$WRITE(*,*) 'n = ', y, DOT_PRODUCT(y, y)
!!$WRITE(*,*) 'x*y = ', DOT_PRODUCT(x, y)

ysquared = DOT_PRODUCT(y, y)
p = 2.0D0*DOT_PRODUCT(x, y)/ysquared
q = (DOT_PRODUCT(x, x)-dist**2)/ysquared

call  QuadEq(p, q, t1, t2, IsReal)

! Durch Rundungsfehler kann die Diskriminante der quadratischen Gleichung
! negativ werden, so dass formal keine reelle Loesung existiert. Sie existiert
! aber doch: p etwas vergroessern, damit die Diskriminate positiv wird und
! nochmal versuchen:
do while (.not. IsReal)
   p = p + 0.001_wp * p
   call  QuadEq(p, q, t1, t2, IsReal)
end do

! Die quadratische Gleichung liefert meist zwei Loesungen (im Idealfall einer
! exakten Rechnung fallen beide Loesungen zusammen), der "wahre" Wert liegt
! genau in der Mitte: Mittelwert bilden.
t = 0.5_wp*(t1 + t2)

!!$WRITE(*,*) 'p, q = ', p, q
!!$WRITE(*,*) 'Loesung der quadratischen Gleichung fuer t: ',  t1, t2, IsReal

call LinePos(g2, t, point)
!write(*,'(a,3(f8.3,tr3))') 'Punkt auf der Geraden fuer t1 : ' , point

if (present(range)) then
   ! Der Parameter "range" wurde uebergeben
   if (t .ge. 0.0_wp .and. t .le. g2%len) then
      ! Der Punkt liegt zwischen g2%start und g2%end
      range = .true.
   else
      ! Der Punkt liegt nicht zwischen g2%start und g2%end
      range = .false.
   end if
end if

if (present(lambda)) then
   lambda = t
end if

end subroutine LineCrossHeight


!---------------------------------------------------------------------
! function LineAngle
!---------------------------------------------------------------------
!
!> @brief Computes the angle between skew lines.
!>
!> <b> angle = LineAngle(g1, g2) </b>
!>
!> This routine computes the angle between skew lines, i.e. the angle
!> between the unit vectors g1%unitvec and g2%unitvec.
!> In general both lines have no common point!
!>
!> @param[in]  g1         derived type defining the first line
!> @param[in]  g2         derived type defining the second line
!> @param[out] LineAngle  angle between the lines g1 and g2 in radian
!>                        (return value)
!>
!---------------------------------------------------------------------
real(wp) function LineAngle(g1, g2)

implicit none

! List of calling arguments:
type (Line), intent(in) :: g1, g2

! List of local variables:
real(wp)    :: CosPhi
!---------------------------------------------------------------------

! b1*b2 = abs(b1)*abs(b2)*cos(phi)
! b1 und b2 sind Einheitsvektoren, daher gilt b1*b2 = cos(phi)
CosPhi = DOT_PRODUCT(g1%unitvec, g2%unitvec)
LineAngle = acos(CosPhi)

end function LineAngle


!---------------------------------------------------------------------
! subroutine CrossProduct
!---------------------------------------------------------------------
!
!> @brief Computes the cross product of two vectors.
!>
!> <b> call CrossProduct(a, b, c) </b>
!>
!> This routine computes the cross product of two vectors in three dimensions.
!>
!> @param[in]  a   first vector
!> @param[in]  b   second vector
!> @param[out] c   c = a x b, cross product of a and b
!>
!---------------------------------------------------------------------
subroutine CrossProduct(a, b, c)

implicit none

! List of calling arguments:
real(wp), dimension(1:3), intent(in)  :: a, b
real(wp), dimension(1:3), intent(out) :: c

! List of local variables:
!---------------------------------------------------------------------

c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = a(3)*b(1) - a(1)*b(3)
c(3) = a(1)*b(2) - a(2)*b(1)

end subroutine CrossProduct


!---------------------------------------------------------------------
! subroutine QuadEq
!---------------------------------------------------------------------
!
!> @brief Solves the quadratic equation.
!>
!> <b> call QuadEq(p, q, x1, x2, IsReal) </b>
!>
!> This routine solves the reduced quadratic equation with real coefficients. \n
!> x**2 + p*x + q = 0 ,  p, q are real numbers
!>
!> @param[in]  p       real coefficient of the quadratic equation
!> @param[in]  q       real coefficient of the quadratic equation
!> @param[out] x1      first real solution of the quadratic equation
!> @param[out] x2      second real solution of the quadratic equation
!> @param[out] IsReal  Are there real solutions of the quadratic equation? \n
!>                     IsReal = .true.  - a real solution exists \n
!>                     IsReal = .false.  - there is no real solution
!>
!---------------------------------------------------------------------
subroutine QuadEq(p, q, x1, x2, IsReal)

implicit none

! List of calling arguments:
real(wp), intent(in)  :: p, q
real(wp), intent(out) :: x1, x2
logical               :: IsReal

! List of local variables:
real(wp)  :: D  ! Diskriminante
real(wp)  :: W  ! Wurzel
!---------------------------------------------------------------------

! Berechnung der Diskriminante
D = (0.5_wp*p)**2 - q
!write(*,*) 'Diskri ', D

! Berechnung der Loesung
if ( abs(D) .lt. 1.0E-8_wp ) then
   ! D = 0, Wurzel = 0, beide reelle Loesungen sind gleich
   x1 = -0.5_wp*p
   x2 = x1
   IsReal = .true.
else if ( D .gt. 0.0_wp ) then
   ! D > 0, Wurzel ist reell, zwei unterschiedliche reelle Loesungen
   W = sqrt(D)
   x1 = -0.5_wp*p + W
   x2 = -0.5_wp*p - W
   IsReal = .true.
else
   ! D < 0, keine reelle Loesung
   IsReal = .false.
end if

end subroutine QuadEq


!---------------------------------------------------------------------
! function PointDist
!---------------------------------------------------------------------
!
!> @brief Computes the distance between two points in 3 dimensions
!>
!> <b> dist = PointDist(a, b) </b>
!>
!> Computes the distance between two points a = (x1, y1, z1) and
!> b = (x2, y2, z2) for Cartesian coordinates in 3 dimensions.
!>
!> @param[in]  a              first point a = (x1, y1, z1)
!> @param[in]  b              second point b = (x2, y2, z2)
!> @param[out] PointDist      euclidean distance between a and b
!>                            (return value)
!>
!---------------------------------------------------------------------
function PointDist(a, b)

implicit none

! List of calling arguments:
real(wp)                              :: PointDist
real(wp), dimension(1:3), intent(in)  :: a, b

! List of local variables:
!---------------------------------------------------------------------

PointDist = sqrt( (b(1)-a(1))**2 + (b(2)-a(2))**2 + (b(3)-a(3))**2 )

end function PointDist


!=====================================================================
! TL / ADJOINT routines follow
!=====================================================================

subroutine LineDef_tl(thisline, start, end, dim,  &
                      thisline_tl, start_tl, end_tl)

implicit none

! List of calling arguments:
type (Line),              intent(inout) :: thisline
real(wp), dimension(1:3), intent(in)    :: start    ! Startpunkt
real(wp), dimension(1:3), intent(in)    :: end      ! Endpunkt
integer,                  intent(in)    :: dim      ! Dimension
real(wp), dimension(1:3), intent(in)    :: start_tl ! Startpunkt
real(wp), dimension(1:3), intent(in)    :: end_tl   ! Endpunkt
type (Line),              intent(out)   :: thisline_tl

! List of local variables:
real(wp), dimension(1:3) :: UnitVec      ! Vektor von Start zu Ende
real(wp), dimension(1:3) :: dUnitVec
!---------------------------------------------------------------------

! Start- und Endpunkt der Geraden in Struktur speichern
thisline%start = start
thisline%end = end
thisline%dim = dim

thisline_tl%start = start_tl
thisline_tl%end = end_tl
thisline_tl%dim = dim

! Richtungsvektor
UnitVec = end - start
dUnitVec = end_tl - start_tl

! Laenge des Strahls
thisline%len = sqrt(UnitVec(1)**2 + UnitVec(2)**2 + UnitVec(3)**2)
thisline_tl%len = ( UnitVec(1)*dUnitVec(1) + UnitVec(2)*dUnitVec(2) +  &
                    UnitVec(3)*dUnitVec(3) ) / thisline%len

! Einheitsvektor in diese Richtung
thisline%unitvec = UnitVec / thisline%len
thisline_tl%unitvec = dUnitVec/thisline%len -                      &
                      UnitVec*thisline_tl%len/thisline%len**2

end subroutine LineDef_tl


subroutine LineDef_ad(thisline, start, end, dim,   &
                      thisline_ad, start_ad, end_ad)

implicit none

! List of calling arguments:
type (Line),              intent(inout) :: thisline
real(wp), dimension(1:3), intent(in)    :: start    ! Startpunkt
real(wp), dimension(1:3), intent(in)    :: end      ! Endpunkt
integer,                  intent(in)    :: dim      ! Dimension
real(wp), dimension(1:3), intent(out)   :: start_ad ! Startpunkt
real(wp), dimension(1:3), intent(out)   :: end_ad   ! Endpunkt
type (Line),              intent(inout) :: thisline_ad

! List of local variables:
real(wp), dimension(1:3) :: UnitVec      ! Vektor von Start zu Ende
real(wp), dimension(1:3) :: DUnitVec
!---------------------------------------------------------------------

! Dieser Teil koennte entfallen, wenn man voraussetzt, dass "thisline"
! die aktuellen Daten enthaelt. Geht das ???

! Start- und Endpunkt der Geraden in Struktur speichern
thisline%start = start
thisline%end = end
thisline%dim = dim

! Richtungsvektor
UnitVec = end - start

! Laenge des Strahls
thisline%len = sqrt(UnitVec(1)**2 + UnitVec(2)**2 + UnitVec(3)**2)

! Einheitsvektor in diese Richtung
thisline%unitvec = UnitVec / thisline%len

! 0000000000000000000000000000000000000000000000000

! Konstante speichern:
thisline_ad%dim = dim

! thisline%unitvec = UnitVec / thisline%len
DUnitVec = thisline_ad%unitvec / thisline%len
thisline_ad%len = thisline_ad%len -                                      &
                  thisline_ad%unitvec(3) * UnitVec(3) / thisline%len**2
thisline_ad%len = thisline_ad%len -                                      &
                  thisline_ad%unitvec(2) * UnitVec(2) / thisline%len**2
thisline_ad%len = thisline_ad%len -                                      &
                  thisline_ad%unitvec(1) * UnitVec(1) / thisline%len**2

! thisline%len = sqrt(UnitVec(1)**2 + UnitVec(2)**2 + UnitVec(3)**2)
DUnitVec = DUnitVec +  thisline_ad%len * UnitVec / thisline%len

! UnitVec = end - start
end_ad = DUnitVec
start_ad = -DUnitVec

thisline_ad%start = start_ad
thisline_ad%end   = end_ad

thisline_ad%len = 0.0_wp
thisline_ad%unitvec = 0.0_wp

end subroutine LineDef_ad


subroutine LinePos_tl(thisline, lambda, point,   &
                      thisline_tl, lambda_tl, point_tl )

implicit none

! List of calling arguments:
type(Line),               intent(in)  :: thisline
real(wp),                 intent(in)  :: lambda
real(wp), dimension(1:3), intent(out) :: point
type(Line),               intent(in)  :: thisline_tl
real(wp),                 intent(in)  :: lambda_tl
real(wp), dimension(1:3), intent(out) :: point_tl

! List of local variables:
integer :: i
!---------------------------------------------------------------------

do i=1, thisline%dim
   point(i) = thisline%start(i) + lambda*thisline%unitvec(i)
   point_tl(i) = thisline_tl%start(i)             +    &
                 thisline%unitvec(i) * lambda_tl  +    &
                 lambda * thisline_tl%unitvec(i)
end do

end subroutine LinePos_tl


subroutine LinePos_ad(thisline, lambda, point,              &
                      thisline_ad, lambda_ad, point_ad)

implicit none

! List of calling arguments:
type(Line),               intent(in)    :: thisline
real(wp),                 intent(in)    :: lambda
real(wp), dimension(1:3), intent(out)   :: point
type(Line),               intent(out)   :: thisline_ad
real(wp),                 intent(inout) :: lambda_ad
real(wp), dimension(1:3), intent(inout) :: point_ad

! List of local variables:
integer :: i
!---------------------------------------------------------------------

do i=thisline%dim, 1, -1
   point(i) = thisline%start(i) + lambda*thisline%unitvec(i)

   thisline_ad%start(i)   = point_ad(i)
   lambda_ad              = lambda_ad + point_ad(i) * thisline%unitvec(i)
   thisline_ad%unitvec(i) =  point_ad(i) * lambda
end do

point_ad = 0.0_wp

end subroutine LinePos_ad


subroutine QuadEq_tl(p, q, x1, x2, IsReal,    &
                     p_tl, q_tl, x1_tl, x2_tl)

implicit none

! List of calling arguments:
real(wp), intent(in)  :: p, q
real(wp), intent(out) :: x1, x2
logical               :: IsReal
real(wp), intent(in)  :: p_tl, q_tl
real(wp), intent(out) :: x1_tl, x2_tl

! List of local variables:
real(wp)  :: D  ! Diskriminante
real(wp)  :: W  ! Wurzel
real(wp)  :: dD, dW
!---------------------------------------------------------------------

! Berechnung der Diskriminante
D = (0.5_wp*p)**2 - q
!write(*,*) 'Diskri ', D
dD = 0.5_wp*p * p_tl - q_tl

! Berechnung der Loesung
if ( abs(D) .lt. 1.0E-8_wp ) then
   ! D = 0, Wurzel = 0, beide reelle Loesungen sind gleich
   x1 = -0.5_wp*p
   x1_tl = -0.5_wp * p_tl
   x2 = x1
   x2_tl = x1_tl
   IsReal = .true.
else if ( D .gt. 0.0_wp ) then
   ! D > 0, Wurzel ist reell, zwei unterschiedliche reelle Loesungen
   W = sqrt(D)
   dW = dD*0.5_wp/sqrt(D)
   x1 = -0.5_wp*p + W
   x1_tl =  -0.5_wp * p_tl + dW
   x2 = -0.5_wp*p - W
   x2_tl = -0.5_wp * p_tl - dW
   IsReal = .true.
else
   ! D < 0, keine reelle Loesung
   IsReal = .false.
   x1_tl = 0.0_wp
   x2_tl = 0.0_wp
end if

end subroutine QuadEq_tl


subroutine QuadEq_ad(p, q, IsReal,              &
                     p_ad, q_ad, x1_ad, x2_ad)

implicit none

! List of calling arguments:
real(wp),  intent(in)   :: p, q
!real(wp), intent(out)  :: x1, x2
logical                 :: IsReal
real(wp), intent(out)   :: p_ad, q_ad
real(wp), intent(inout) :: x1_ad, x2_ad

! List of local variables:
real(wp)  :: D  ! Diskriminante
!real(wp)  :: W  ! Wurzel
real(wp)  :: dW, dD
!---------------------------------------------------------------------

! Berechnung der Diskriminante
D = (0.5_wp*p)**2 - q

! Berechnung der Loesung
if ( abs(D) .lt. 1.0E-8_wp ) then
   ! D = 0, Wurzel = 0, beide reelle Loesungen sind gleich
   IsReal = .true.
   ! x2 = x1
   x1_ad = x1_ad + x2_ad
   ! x1 = -0.5D0*p
   p_ad = -0.5_wp*x1_ad
   dD = 0.0_wp
   x1_ad = 0.0_wp
   x2_ad = 0.0_wp
else if ( D .gt. 0.0_wp ) then
   ! D > 0, Wurzel ist reell, zwei unterschiedliche reelle Loesungen
   IsReal = .true.
   ! x2 = -0.5D0*p - W
   p_ad =  -0.5_wp*x2_ad
   dW   = - x2_ad
   ! x1 = -0.5D0*p + W
   p_ad = p_ad - 0.5_wp*x1_ad
   dW   = dW + x1_ad
   ! W = sqrt(D)
   dD   = 0.50_wp * dW / sqrt(D)
   x1_ad = 0.0_wp
   x2_ad = 0.0_wp
else
   ! D < 0, keine reelle Loesung
   IsReal = .false.
   p_ad = 0.0_wp
   dD = 0.0_wp
end if

! D = (0.5D0*p)**2 - q
p_ad = p_ad + 0.5_wp*p*dD
q_ad = -dD

end subroutine QuadEq_ad

!=====================================================================
end module mo_std_vector
!=====================================================================
