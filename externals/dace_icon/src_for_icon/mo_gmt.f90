!
!+ Write files to be read by GMT plot routines
!
MODULE mo_gmt
!
! Description:
!   Write files to be read by GMT plot routines
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_7         2009/08/24 Andreas Rhodin
!  adapt to GMT version 4.3.1
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  cleanup (remove unused variables)
! V1_20        2012-06-18 Andreas Rhodin
!  change format for temporary plot-file
! V1_26        2013/06/27 Jaison Ambadan
!  option to plot filled polygons
! V1_27        2013-11-08 Jaison Ambadan
!  use colormap (instead of black) for un-filled contour plot
! V1_37        2014-12-23 Harald Anlauf
!  replace iarea (integer area boundary) by real darea
! V1_42        2015-06-08 Andreas Rhodin
!  adapt to COSMO grid; option to prescribe scaling (min/max value)
! V1_48        2016-10-06 Harald Anlauf
!  ICON grid coordinates: use radians instead of degrees
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2002  original code
! Oliver Schmid   DWD  2006  changes for GMT4
! Harald Anlauf   DWD  2007  fixes for GMT 4.2
!==============================================================================
  use mo_kind,          only: wp, sp                 ! kind parameters
  use mo_system       !,only: system                 ! platform specific
  use mo_physics,       only: pi,                   &! 3.14159...
                              r2d                    ! factor radians->degree
  use mo_atm_state,     only: t_atm                  ! atmospheric state type
  use mo_atm_grid,      only: t_grid,               &! grid meta data type
                              phirot2phi,           &! rotated longitude -> lon
                              rlarot2rla             ! rotated latitude  -> lat
  use mo_time,          only: t_time,               &! date&time data type
                              iyyyy, imm, idd, ihh, &! time conversion routines
                              cdate, ctime           !
  use mo_mpi_dace,      only: dace                   ! MPI group info
  use mo_atm_transp,    only: gather_level           !
  use mo_wmo_tables,    only: DWD6_ICOSAHEDRON,     &
                              DWD6_ICON,            &
                              WMO6_GAUSSIAN,        &
                              WMO6_LATLON,          &
                              WMO6_ROTLL
  use mo_fortran_units, only: get_unit_number, return_unit_number
  use mo_ico_grid,      only: nspoke2

  implicit none

  private
  public :: write_gmt_level
  public :: write_script
  public :: atm_to_gmt
  public :: run_script

  ! constants

  REAL(wp), PARAMETER :: pi2      = 2.0_wp*pi
  REAL(wp), PARAMETER :: pid2     = 0.5_wp*pi
  REAL(wp), PARAMETER :: pid180i  = 180.0_wp/pi

  ! module variables

  real(wp)           :: minv  = 0._wp
  real(wp)           :: maxv  = 0._wp
  real(wp)           :: minc  = 0._wp
  real(wp)           :: maxc  = 0._wp
  real(wp)           :: mean  = 0._wp
  real(wp)           :: delt  = 0._wp
  real(wp)           :: stdev = 0._wp
  real(wp)           :: rms   = 0._wp
  logical            :: shade = .false.
  integer            :: level = 0
  character(len=16)  :: nam   = ''
  character(len=10)  :: date  = '0000-00-00'
  character(len=8)   :: time  = '00:00:00'
  character(len=1)   :: ptype = 'b'
  logical            :: lgme
! logical            :: lico
! logical            :: lgau
! logical            :: lreg
  logical            :: ltr
  logical            :: leps = .true.

contains
!==============================================================================
  subroutine atm_to_gmt (atm, basename)
  type (t_atm)               ,intent(in) :: atm
  character(len=*) ,optional ,intent(in) :: basename

    integer           :: i, k
    character(len=12) :: name
    character(len=64) :: fullname
    character(len=3)  :: clev

    if (atm% grid% gridtype /= DWD6_ICOSAHEDRON) return
    do i=1,size(atm% m)
      if (.not. atm% m(i)% i% alloc) cycle
      name = atm% m(i)% i% name
      level = 0
      if (atm% m(i)%i% lb(3) == atm% m(i)%i% ub(3)) then
        fullname = name
        if(present(basename)) fullname = trim(basename)//'_'//trim(name)
        call write_gmt_level (atm% m(i)% ptr(:,:,1,:), name, atm% grid, &
                              atm% time,                                &
                              basename=fullname, code=0, lgather=.true.,&
                              ptype=ptype, lev=level                    )
        call write_script (fullname, ptype)
      else
        name = atm% m(i)% i% name
        do k=atm% m(i)%i% lb(3), atm% m(i)%i% ub(3)
          level = k
          write(clev,'(i3.3)') level
          fullname = name
          if(present(basename)) fullname = trim(basename)//'_'//&
                                           trim(name)//'_'//clev
          call write_gmt_level (atm% m(i)% ptr(:,:,k,:), name, atm% grid, &
                                atm% time,                                &
                                basename=fullname, code=0, lgather=.true.,&
                                ptype=ptype, lev=level                    )
          call write_script (fullname, ptype)
        end do
      endif
    end do
  end subroutine atm_to_gmt
!------------------------------------------------------------------------------
  subroutine write_gmt_level (xx, name, grid, ver_time, basename, code,   &
                              lgather, ptype, lev, minmax, delta, shaded, &
                              ltri, area, mask)

  character(len=*) ,intent(in)           :: name
  type(t_grid)     ,intent(in)           :: grid
  real(wp)         ,intent(in) ,target   :: xx (grid% lbg(1):,:,:)
  type(t_time)     ,intent(in)           :: ver_time
  character(len=*) ,intent(in)           :: basename
  integer          ,intent(in)           :: code
  logical          ,intent(in)           :: lgather
  character(len=*) ,intent(in)           :: ptype
  integer          ,intent(in) ,optional :: lev
  real(wp)         ,intent(in) ,optional :: minmax (2)
  real(wp)         ,intent(in) ,optional :: delta
  logical          ,intent(in) ,optional :: shaded
  logical          ,intent(in) ,optional :: ltri
  real(wp)         ,intent(in) ,optional :: area (4)
  logical ,pointer ,intent(in) ,optional :: mask (:,:,:)

    real(wp) ,pointer :: xg    (:,:,:)
    integer           :: ni
    integer           :: nd
    integer           :: day,month,year,hour,minute
    integer           :: l, i, j, n, nin, ic
    integer           :: i2, j2, i3, j3, l3, il, ib, iv
    real(wp)          :: x, y, z, xd
    real(wp)          :: xpole
!   real(sp)          :: xsp,ysp,zsp
    integer           :: lbg(4), ubg(4)
    integer           :: iunit, iunitt !, iunitx
!   integer           :: irecl
    real(wp)          :: s, s2
    integer           :: ij   (grid% lbg(1):grid% ubg(1),&
                               grid% lbg(2):grid% ubg(2),&
                               grid% lbg(4):grid% ubg(4))
    logical           :: lmask(grid% lbg(1):grid% ubg(1),&
                               grid% lbg(2):grid% ubg(2),&
                               grid% lbg(4):grid% ubg(4))
!   real(sp)          :: eps = 0.01_sp
    real(wp)          :: darea (4)      ! Plotting area
    character(len=18) :: ch
    real(wp)          :: cx(4), cy(4), cz(4), vx(4), vy(4)
    real(wp)          :: dx, dy
    !--------------------
    ! derive array bounds
    !--------------------
    ni     = grid% ni
    nd     = grid% nd
    lbg    = grid% lbg
    ubg    = grid% ubg
    !---------------------
    ! set module variables
    !---------------------
    ltr    = .false.; if (present(ltri)) ltr = ltri
    lgme   = grid% gridtype == DWD6_ICOSAHEDRON
!   lico   = grid% gridtype == DWD6_ICON
!   lgau   = grid% gridtype == WMO6_GAUSSIAN
!   lreg   = grid% gridtype == WMO6_LATLON
    date   = cdate (ver_time)
    time   = ctime (ver_time)

    year   = iyyyy (ver_time)
    month  = imm   (ver_time)
    day    = idd   (ver_time)
    hour   = ihh   (ver_time)
    minute = 0

    level = 0;                     if (present (lev )) level = lev
    darea = (/-180,+180,-90,+90/); if (present (area)) darea = area
    lmask = .true.
    if   (present    (mask)) then
      if (associated (mask)) lmask = mask
    endif

    if (lgather) then
      nullify (xg)
      if(dace% lpio) allocate (xg (lbg(1):ubg(1),lbg(2):ubg(2),lbg(4):ubg(4)))
      call gather_level (xg, xx, grid% dc, dace% pio)
    else
      xg => xx
    endif
    nam  = name
    shade = .false.; if (present(shaded)) shade = shaded

    if(dace% lpio) then

      iunit  = get_unit_number()
!     iunitx = get_unit_number()
      iunitt = get_unit_number()
!     inquire (iolength = irecl) xsp, ysp, zsp
      OPEN(iunit,  FILE=trim(basename)//'.xyz')
!     OPEN(iunitx, FILE=trim(basename)//'.bin',access='direct',recl=irecl)
      if (ltr) OPEN(iunitt, FILE=trim(basename)//'.ijk')

      WRITE (iunit,'(80("#"))')
      WRITE (iunit,'(a,2i5)')     '# Triangular grid: ', ni, nd
      WRITE (iunit,'(a,i0,1x,a)') '# Parameter: ', code, name
      WRITE (iunit,'(a,i2.2,a,i2.2,a,i4.4,a,i2.2,a,i2.2)') &
        '# Date: ', day,'.',month,'.',year,' ', hour, ':', minute
      WRITE (iunit,'(80("#"))')

      n   = 0
      nin = 0
      s   = 0._wp
      s2  = 0._wp
      ij  = 0
      DO l = lbg(4), ubg(4)
        if (lgme) xpole = pid180i * grid% rlon(ni,ni+1,1,l)
        DO i = lbg(1), ubg(1)
          DO j = lbg(2), ubg(2)
            !--------------------------
            ! index (line) in plot file
            !--------------------------
            n = n + 1
            ij (i,j,l) = n-1
            !---------------------------------
            ! coordinates and value to contour
            !---------------------------------
            x = pid180i * grid% rlon(i,j,1,l)
            y = pid180i * grid% rlat(i,j,1,l)
            z = xg(i,j,l)                      ! store the value

            if (x > darea(2)) x = x - 360._wp
            if (x < darea(1)) x = x + 360._wp
            if (x < darea(1)) lmask(i,j,l) = .false.
            if (x > darea(2)) lmask(i,j,l) = .false.
            if (y < darea(3)) lmask(i,j,l) = .false.
            if (y > darea(4)) lmask(i,j,l) = .false.

            if (lmask(i,j,l)) then
              nin = nin + 1
              s  = s  + z
              s2 = s2 + z * z
            endif
            !---------------------------------------------------
            ! special handling of the poles for triangular grids
            !---------------------------------------------------
            if (lgme) then
               if (grid%marr(4,i,j,l)/=l) lmask(i,j,l) = .false.
               if (i==0.and.j==1)         lmask(i,j,l) = .true.
               if (i==0.and.j==1)         x = xpole
               if ((x-xpole)>180._wp)     x = x - 360._wp
               if ((xpole-x)>180._wp)     x = x + 360._wp
            endif

            !----------------
            ! write plot file
            !----------------

!           xsp=x; ysp=y; zsp=z

              if (ptype == 'b' .or. ptype == 'p') then
                if (lmask(i,j,l)) then
                  WRITE (ch ,'(F18.8)') z
                  WRITE (iunit ,'(A4,F18.8)') "> -Z", z
                  cx(1) = grid% xnglob(i,j,1,l)
                  cy(1) = grid% xnglob(i,j,2,l)
                  cz(1) = grid% xnglob(i,j,3,l)
                  select case (grid% gridtype)
                  case (DWD6_ICOSAHEDRON)
                    do ic=1,7
                      i2 = i + nspoke2 (1,mod(ic-1,6)+1)
                      j2 = j + nspoke2 (2,mod(ic-1,6)+1)

                      i3 = grid% marr (2,i2,j2,l)
                      j3 = grid% marr (3,i2,j2,l)
                      l3 = grid% marr (4,i2,j2,l)

                      cx(2) = grid% xnglob(i3,j3,1,l3)
                      cy(2) = grid% xnglob(i3,j3,2,l3)
                      cz(2) = grid% xnglob(i3,j3,3,l3)

                      i2 = i + nspoke2 (1,mod(ic,6)+1)
                      j2 = j + nspoke2 (2,mod(ic,6)+1)

                      i3 = grid% marr (2,i2,j2,l)
                      j3 = grid% marr (3,i2,j2,l)
                      l3 = grid% marr (4,i2,j2,l)

                      cx(3) = grid% xnglob(i3,j3,1,l3)
                      cy(3) = grid% xnglob(i3,j3,2,l3)
                      cz(3) = grid% xnglob(i3,j3,3,l3)

                      if (cx(2)==cx(3).and.cy(2)==cy(3).and.cz(2)==cz(3)) cycle

                      cx(4) = sum (cx(1:3))/3
                      cy(4) = sum (cy(1:3))/3
                      cz(4) = sum (cz(1:3))/3

                      x = pid180i * atan2 (cy(4), cx(4))
                      y = pid180i * asin  (cz(4))

                      if (x > darea(2)) x = x - 360._wp
                      if (x < darea(1)) x = x + 360._wp

                      WRITE (iunit ,'(2(F18.8))') x,  y
!                     WRITE (iunitx, rec=n) xsp, ysp, zsp               ! binary
                    enddo
                  case (DWD6_ICON)
                    xd = pid180i * grid% rlon(i,j,1,l)
                    if (xd > darea(2)) xd = xd - 360._wp
                    if (xd < darea(1)) xd = xd + 360._wp
                    do iv = 1,3
                       il = grid% icongrid% patch% cells% vertex_idx(i,j,iv)
                       ib = grid% icongrid% patch% cells% vertex_blk(i,j,iv)
                       x  = grid% icongrid% patch% verts% vertex    (il,ib)% lon * r2d
                       y  = grid% icongrid% patch% verts% vertex    (il,ib)% lat * r2d
                       if (x > darea(2)) x = x - 360._wp
                       if (x < darea(1)) x = x + 360._wp
                       ! Avoid GMT plotting strange triangles
                       if      ((x - xd) >  180._wp) then
                          x = x - 360._wp
                       else if ((x - xd) < -180._wp) then
                          x = x + 360._wp
                       end if
                       vx(iv) = x
                       vy(iv) = y
                    end do
                    vx(4) = vx(1); vy(4) = vy(1)
!                   x = minval (vx); y = maxval (vx)
!                   if (abs (x - y) > 180) then
!                      print *, xd,':', real (vx(1:3))
!                   end if
                    do iv = 1,4
                       WRITE (iunit ,'(2(F16.8))') vx(iv), vy(iv)
                    end do
                  case (WMO6_LATLON, WMO6_ROTLL)
                    x  = grid% dlon (i)
                    y  = grid% dlat (j)
                    dx = grid% di / 2._wp
                    dy = grid% dj / 2._wp
                    vx(1) = x - dx; vy(1) = y - dy
                    vx(2) = x + dx; vy(2) = y - dy
                    vx(3) = x + dx; vy(3) = y + dy
                    vx(4) = x - dx; vy(4) = y + dy
                    if (grid% gridtype == WMO6_ROTLL) then
                      do iv = 1,4
                        x = rlarot2rla(vy(iv),vx(iv),grid% dlatr,grid% dlonr,0._wp)
                        y = phirot2phi(vy(iv),vx(iv),grid% dlatr,grid% dlonr,0._wp)
                        vx(iv) = x
                        vy(iv) = y
                      end do
                    endif
                    WRITE (iunit ,'(2(F16.8))') vx(1), vy(1)
                    WRITE (iunit ,'(2(F16.8))') vx(2), vy(2)
                    WRITE (iunit ,'(2(F16.8))') vx(3), vy(3)
                    WRITE (iunit ,'(2(F16.8))') vx(4), vy(4)
                    WRITE (iunit ,'(2(F16.8))') vx(1), vy(1)
                  case default
                    print *, 'mo_gmt: unsupprted grid type: ', grid% gridtype
                  end select
                end if
              else
                if (ptype == 'c') then
                    if (lmask(i,j,l)) then
                      WRITE (iunit ,'(3(F18.8:1X),L1)') x,  y,  z, lmask(i,j,l) ! ASCII
!                     WRITE (iunitx, rec=n)        xsp, ysp, zsp               ! binary
                      if (i==lbg(1).or.j==lbg(2)) cycle
                    endif
                endif
              endif
            !
            if (i==lbg(1).or.j==lbg(2)) cycle
            !------------------------------------------------
            ! add a point at possible plot boundaries for GMT
            !------------------------------------------------
!!           if (lgme.and.j-i==1) then
!            if (lgme) then
!              n = n + 1
!              x = pid180i *  grid% rlon(i  ,j,1,l)
!              y = pid180i *  grid% rlat(i-1,j,1,l)
!              z =           (       xg (i-1,j,  l) + xg (i,j-1,  l))/2._wp
!              xsp=x; ysp=y; zsp=z
!              WRITE (iunit , '(3(F18.8:1X),a)') x+eps,   y,   z, '-'  ! ASCII
!              WRITE (iunit , '(3(F18.8:1X),a)') x-eps,   y,   z, '-'  ! ASCII
!!             WRITE (iunit , '(3(F18.8:1X),a)') x,       y,   z, '-'  ! ASCII
!              WRITE (iunitx, rec=n)        xsp+eps, ysp, zsp          ! binary
!              WRITE (iunitx, rec=n)        xsp-eps, ysp, zsp          ! binary
!!             WRITE (iunitx, rec=n)        xsp,     ysp, zsp          ! binary
!            endif
            !--------------------
            ! write triangulation
            !--------------------
            if (ltr) then
!            if (lgme.and.j-i==1) then
!              !------------------------------------------------
!              ! add a point at possible plot boundaries for GMT
!              !------------------------------------------------
!              WRITE (iunitt, '(3I6,a)') ij(i-1,j-1,l), n-2        , ij(i-1,j,l),'-'
!              WRITE (iunitt, '(3I6,a)') ij(i-1,j-1,l), ij(i,j-1,l), n-1        ,'-'
!              WRITE (iunitt, '(3I6,a)') ij(i  ,j  ,l), n-2        , ij(i-1,j,l),'-'
!              WRITE (iunitt, '(3I6,a)') ij(i  ,j  ,l), ij(i,j-1,l), n-1        ,'-'
!            else
              if (any (lmask(i-1:i,j-1:j,l))) then
                WRITE (iunitt, '(3I6)') ij(i-1,j-1,l), ij(i,j-1,l), ij(i-1,j,l)
                WRITE (iunitt, '(3I6)') ij(i  ,j  ,l), ij(i,j-1,l), ij(i-1,j,l)
              endif
!            endif
            endif
            !-----------
            ! statistics
            !-----------
          END DO
        END DO
      END DO
      !-----------
      ! statistics
      !-----------

      ! TODO: adjust/round interval boundaries minv, maxv so that only
      !       a limited but still useful number of digits is retained
      minv = minval (xg, mask=lmask); minc = minv
      maxv = maxval (xg, mask=lmask); maxc = maxv
      if (present (minmax)) then
        if (minmax(2) > minmax(1)) then
          minc = minmax(1)
          maxc = minmax(2)
        endif
      endif
      if        (minc == maxc) then
        if      (minc > 0._wp) then
          minc =  0._wp
        else if (minc < 0._wp) then
          maxc =  0._wp
        else
          maxc =  1._wp
          minc = -1._wp
        endif
      endif
      delt = (maxc-minc)/12
      if (present(delta)) then
        if(delta>0._wp) delt = delta
      endif

      mean  = s / n
      rms   = sqrt (s2 / n)
      stdev = sqrt (max(0._wp, s2/n -mean**2))
      !------------
      ! close files
      !------------
      CLOSE (iunit)
!     CLOSE (iunitx)
      if (ltr) CLOSE (iunitt)
      call return_unit_number (iunit)
!     call return_unit_number (iunitx)
      call return_unit_number (iunitt)
      !--------
      ! cleanup
      !--------
      if (lgather) deallocate (xg)

    endif

  end subroutine write_gmt_level
!------------------------------------------------------------------------------
  subroutine write_script (basename, ptype, area)
  character(len=*) ,intent(in)           :: basename
  character(len=*) ,intent(in)           :: ptype
  real(wp)         ,intent(in) ,optional :: area(4)

    integer           :: iunit
    integer           :: i
    character(len=12) :: clev
    real(wp)          :: anno
!   character(len=13) :: blue  = '   0   0 255 '
!   character(len=13) :: red   = ' 255   0   0 '
!   character(len=13) :: green = '   0 255   0 '
!   character(len=13) :: color
    character(len=32) :: bopt
    integer           :: iarea(4)
!   character(len=15) :: ctri  = '' ! or ' -T${BASE}.ijk '

    if(.not.dace% lpio) return

    iunit  = get_unit_number()
    open  (iunit,file=trim(basename)//'.sh')

    i = index (basename, '/', back=.true.)              ! Find last '/'
    write (clev,"(""'level "",i4,""'"")") level

    if (delt == 0._wp) delt = 1._wp
    anno = 5 * delt

    write (iunit,'()')

!   bopt =' -R-180/+180/-90/+90'
    bopt =' -Rg'
    if (present(area)) then
      iarea = nint (area)
      if (all (abs (area - iarea) < 1.e-2_wp)) then
         ! Old format (rounded to nearest integers)
         write(bopt,"(' -R',SP,I4.3,'/',I4.3,'/',I3.2,'/',I3.2)") iarea
      else
         write(bopt,"(' -R',SP,F0.2,'/',F0.2,'/',F0.2,'/',F0.2)") area
      end if
!!print *, "bopt='",trim(bopt),"' len=", len_trim (bopt)
    endif

!   ctri  = ''; if (ltr) ctri  = ' -T${BASE}.ijk '

    call set_var_r ("MIN"   ,minv)
    call set_var_r ("MAX"   ,maxv)
    call set_var_r ("MINC"  ,minc)
    call set_var_r ("MAXC"  ,maxc)
    call set_var_r ("MEAN"  ,mean)
    call set_var_r ("STDEV" ,stdev)
    call set_var_r ("RMS"   ,rms)
    call set_var_r ("DELT"  ,delt)
    call set_var_r ("ANNO"  ,anno)
    call set_var_s ("BASE"  ,trim(basename(i+1:)))
    call set_var_s ("NAME"  ,trim(nam))
    call set_var_s ("DATE"  ,date)
    call set_var_s ("TIME"  ,time)
    call set_var_s ("LEVEL" ,clev)

    write (iunit,'(a)') "gmt gmtset PS_MEDIA A4                           "
    write (iunit,'(a)') "gmt gmtset PS_PAGE_ORIENTATION portrait          "
    write (iunit,'(a)') "gmt gmtset MAP_GRID_PEN_PRIMARY   0.05p,gray,.-  "
    write (iunit,'(a)') "gmt gmtset MAP_GRID_PEN_SECONDARY 0.03p,gray,.-  "
    write (iunit,'(a)') "gmt gmtset FONT_ANNOT_PRIMARY 5p                 "
    write (iunit,'(a)') "gmt gmtset FONT_ANNOT_SECONDARY 4p               "
    write (iunit,'(a)') "gmt gmtset MAP_ANNOT_OFFSET_PRIMARY   0.1c       "
    write (iunit,'(a)') "gmt gmtset MAP_ANNOT_OFFSET_SECONDARY -2p        "
    write (iunit,'(a)') "gmt gmtset FONT_TITLE 6p MAP_ANNOT_OBLIQUE 0     "
    write (iunit,'(a)') "gmt gmtset FONT_LABEL 5p                         "
    write (iunit,'(a)') "gmt gmtset MAP_FRAME_PEN 0.2p                    "
    write (iunit,'(a)') "gmt gmtset MAP_TITLE_OFFSET -0.08i               "
    write (iunit,'(a)') "gmt gmtset FORMAT_FLOAT_OUT %4.1f                "
    write (iunit,'(a)') "gmt gmtset MAP_FRAME_TYPE plain                  "
    write (iunit,'(a)') "gmt gmtset MAP_TICK_LENGTH 0.1c                  "
!   write (iunit,'(a)') "gmt gmtset DOTS_PR_INCH 150                      "
   !
    write (iunit,'()')
    write (iunit,'(a)') "echo 2.0 7.25 6 0 1 BL ${NAME}        > ${BASE}.txt"
    write (iunit,'(a)') "echo 4.0 7.25 6 0 1 BL ${LEVEL}      >> ${BASE}.txt"
    write (iunit,'(a)') "echo 6.0 7.25 6 0 1 BL ${DATE}       >> ${BASE}.txt"
    write (iunit,'(a)') "echo 8.0 7.25 6 0 1 BL ${TIME}       >> ${BASE}.txt"
    write (iunit,'(a)') "echo 2.0 7.00 6 0 8 BL min: ${MIN}   >> ${BASE}.txt"
    write (iunit,'(a)') "echo 2.0 6.75 6 0 8 BL max: ${MAX}   >> ${BASE}.txt"
    write (iunit,'(a)') "echo 4.0 7.00 6 0 8 BL mean :${MEAN} >> ${BASE}.txt"
    write (iunit,'(a)') "echo 4.0 6.75 6 0 8 BL delt :${DELT} >> ${BASE}.txt"
    write (iunit,'(a)') "echo 6.0 7.00 6 0 8 BL stdev:${STDEV}>> ${BASE}.txt"
    write (iunit,'(a)') "echo 6.0 6.75 6 0 8 BL rms  :${RMS}  >> ${BASE}.txt"
    write (iunit,'()')

    !======
    ! COLOR
    !======
    if (ptype == 'b' .or. ptype == 'p') then
       write (iunit,'(a)') "gmt gmtset FORMAT_FLOAT_OUT %13.6f"

       write (iunit,'(a)') &
            "gmt makecpt -D -Crainbow -T${MINC}/${MAXC}/${DELT} > ${BASE}.cpt"

       write (iunit,'(a)') "gmt gmtset FORMAT_FLOAT_OUT %4.1f"

       if (shade) then
         write (iunit,'(a,a,a,a)') &
           "gmt psxy ", trim (bopt), " ${BASE}.xyz -C${BASE}.cpt", &
            " -JQ0/9 -Ba60/a30/WeSn -L -P -Xc -Yc -K > TMP.ps"
       else
         write (iunit,'(a,a,a,a)') &
           "gmt psxy ", trim (bopt), " ${BASE}.xyz -C${BASE}.cpt", &
            " -JQ0/9 -Ba60/a30/WeSn    -P -Xc -Yc -K > TMP.ps"
       endif

       write (iunit,'(a,a,a)') &
            "gmt pscoast ", trim (bopt), " -JQ0/9 -P -Dc -W0.01p,gray -O -K >> TMP.ps"

       write (iunit,'(a)') "gmt gmtset FORMAT_FLOAT_OUT %4.2f"

       write (iunit,'(a)') &
            "gmt psscale -D1.8i/-0.2i/3i/0.13ih -S -C${BASE}.cpt -O -K >> TMP.ps"

       write (iunit,'(a)') "gmt gmtset FORMAT_FLOAT_OUT %4.2f"

       write (iunit,'(a)') &
            "gmt pstext ${BASE}.txt -R0/12/0/8 -Jx1 -N -Y-0.85i -O >> TMP.ps"

    elseif (ptype == 'c') then

       write (iunit,'(a)') "gmt gmtset FORMAT_FLOAT_OUT %13.6f"

       write (iunit,'(a)') &
            "gmt makecpt -D -Crainbow -T${MINC}/${MAXC}/${DELT} > ${BASE}.cpt"

       write (iunit,'(a)') "gmt gmtset FORMAT_FLOAT_OUT %4.1f"

       write (iunit,'(a,a,a)') &
             "gmt psbasemap ", trim (bopt), " -JQ0/9 -Ba60/a30/WeSn -Xc -Yc -K > TMP.ps"

       if (shade) then
          write (iunit,'(a,a,a,a,a)') &
             "gmt pscontour ${BASE}.xyz ", trim (bopt), " -JQ0/9 -I -C${BASE}.cpt ", &
                   " -A+a0+f3+jBB+c0.5p/1.25p+gwhite+p0.01p,white+r1.0p -W0.2p",     &
                   " -O -K >> TMP.ps"
       else

           write (iunit,'(a,a,a,a,a)') &
             "gmt pscontour ${BASE}.xyz ", trim (bopt), " -JQ0/9    -C${BASE}.cpt ", &
                   " -A+a0+f3+jBB+c0.5p/1.25p+gwhite+p0.01p,white+r1.0p -W+0.2p",    &
                   " -O -K >> TMP.ps"
       endif

       write (iunit,'(a,a,a)') &
            "gmt pscoast ", trim (bopt), " -JQ0/9 -P -Dc -W0.01p,gray -O -K >> TMP.ps"

       write (iunit,'(a)') "gmt gmtset FORMAT_FLOAT_OUT %4.2f"

!      if (shade) write (iunit,'(a)') &
       write (iunit,'(a)') &
            "gmt psscale -D1.8i/-0.2i/3i/0.13ih -S -C${BASE}.cpt -O -K >> TMP.ps"

       write (iunit,'(a)') &
            "gmt pstext ${BASE}.txt -R0/12/0/8 -Jx1 -N -Y-0.85i -O >> TMP.ps"

!!! -T${BASE}.ijk       Name of file with network information
!!! -I                  Color the triangles using the color palette table.
!!! -L                  Draw the underlying triangular mesh
!!! -W+0.5p             Select contouring and set contour pen attributes.
!!! -B30g15/B30g15WESN  Sets map boundary annotation and tickmark intervals
    end if

    write (iunit,'(a)') "\mv TMP.ps ${BASE}.ps"
    if (leps) write (iunit,'(a)') "gmt psconvert -A0.06i -Te ${BASE}.ps 2>/dev/null"
    write (iunit,'()')
    close (iunit)
    call return_unit_number (iunit)

    call system ('chmod +x '//trim(basename)//'.sh')

  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine set_var_s (name, string)
    character(len=*) ,intent(in) :: name
    character(len=*) ,intent(in) :: string

      write (iunit,'(a,"=",a)') name,string

    end subroutine set_var_s
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine set_var_r (name, value)
    character(len=*) ,intent(in) :: name
    real(wp)         ,intent(in) :: value

      character(len=10) :: c
      write (c,'(g10.3)') value
      write (iunit,'(a,"=",a)') name,trim(adjustl(c))
    end subroutine set_var_r
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine write_script
!------------------------------------------------------------------------------
  subroutine run_script (basename, call_gv)
  character(len=*)  ,intent(in) :: basename
  logical, optional ,intent(in) :: call_gv

    integer :: i
    logical :: gv
    i = index (basename, '/', back=.true.)              ! Find last '/'
    gv = .false.; if (present(call_gv)) gv = call_gv

    if(.not.dace% lpio) return

    if (i/=0) then
      call system ("(cd "//trim(basename(1:i))//";./"//&
                           trim(basename(i+1:))//".sh)")
    else
      call system ("./"//trim(basename(i+1:))//".sh")
    endif
    if (gv) call system ("gv -noantialias -watch "//trim(basename)//".eps &")

  end subroutine run_script
!==============================================================================
end module mo_gmt
