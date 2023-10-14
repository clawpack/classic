module abl_module

  implicit none
  integer, parameter :: NO_ABL = 0
  integer, parameter :: TRIG = 1
  integer, parameter :: AC = 2 ! Appelo & Colonius
  integer, parameter :: PS = 3 ! Petersson & Sjogreen

  integer, parameter :: MAX_PARAM = 3

  real(kind=8), parameter :: PIH = 2.d0*datan(1.d0)


  integer :: abl_type = 0
  integer :: p_AC, q_AC
  real(kind=8) :: eps_AC, eps_PS
  real(kind=8) :: xdepth(2), ydepth(2), xpos(2), ypos(2)
  logical :: initialized = .false.

contains

  subroutine initialize()

    implicit none

    real(kind=8) :: xlower, xupper, ylower, yupper

    ! Read in ABL data
    ! open the unit with new routine from Clawpack 4.4 to skip over
    ! comment lines starting with #:
    call opendatafile(7, 'abl.data')

    ! Read in parameters:
    read(7,*) abl_type
    if (abl_type == TRIG) then
      print *, "ABL is trigonometric"
    else if (abl_type == AC) then
      read(7,*) eps_AC
      read(7,*) p_AC
      read(7,*) q_AC
      print *, "ABL is Appelo/Colonius with (eps,p,q)=", eps_AC, p_AC, q_AC
    else if (abl_type == PS) then
      read(7,*) eps_PS
      print *, "ABL is Petersson/Sjogreen with eps=", eps_PS
    end if
    read(7,*) xdepth(1), ydepth(1)
    read(7,*) xdepth(2), ydepth(2)
    close(7)

    ! Read in grid parameters
    ! open the unit with new routine from Clawpack 4.4 to skip over
    ! comment lines starting with #:
    call opendatafile(7, 'claw.data')

    read(7,*) xlower ! num_dim
    read(7,*) xlower, ylower ! lower
    read(7,*) xupper, yupper ! upper
    close(7)

    ! Compute ABL positions
    xpos(1) = xlower + xdepth(1)
    xpos(2) = xupper - xdepth(2)
    ypos(1) = ylower + ydepth(1)
    ypos(2) = yupper - ydepth(2)

    initialized = .true.

  end subroutine initialize

  subroutine set_scaling_factor(xlower, ylower, dx, dy, mx, my, num_ghost, &
                                ablX_M, ablX_J, ablX_R, ablY_M, ablY_J, ablY_R, &
                                i_abl_lower, i_abl_upper, j_abl_lower, j_abl_upper)

    implicit none

    integer, intent(in) :: mx, my, num_ghost
    real(kind=8), intent(in) :: xlower, ylower, dx, dy
    integer, intent(out) :: i_abl_lower, i_abl_upper, j_abl_lower, j_abl_upper
    real(kind=8), intent(out) :: ablX_M(1-num_ghost:mx+num_ghost)
    real(kind=8), intent(out) :: ablX_J(1-num_ghost:mx+num_ghost)
    real(kind=8), intent(out) :: ablX_R(1-num_ghost:mx+num_ghost)
    real(kind=8), intent(out) :: ablY_M(1-num_ghost:my+num_ghost)
    real(kind=8), intent(out) :: ablY_J(1-num_ghost:my+num_ghost)
    real(kind=8), intent(out) :: ablY_R(1-num_ghost:my+num_ghost)

    integer :: i, j
    real(kind=8) :: xcenter, ycenter, xedge, yedge

    ! initialize module if needed
    if (.not. initialized) then
      call initialize()
    end if

    ! x direction

    ! determine which indices are in the layer and loop over those
    i_abl_lower = int((1.d0 + 1.d-14)*(xpos(1)-xlower)/dx)
    i_abl_upper = 1 + int((1.d0 + 1.d-14)*(xpos(2)-xlower)/dx)

    do i=1-num_ghost,i_abl_lower

      xcenter = xlower + (i-0.5d0)*dx

      ! Calculate inverse of Jacobian at cell edges
      xedge = max(xcenter-0.5d0*dx, xpos(1)-xdepth(1))
      xedge = min(xedge, xpos(2)+xdepth(2))
      ablX_J(i) = 1.d0/(1.d0 + gp(xedge,1))

      ! Calculate grid map scaling and ratio at cell centers
      xcenter = max(xcenter, xpos(1)-xdepth(1)+0.5d0*dx)
      xcenter = min(xcenter, xpos(2)+xdepth(2)-0.5d0*dx)
      if (abl_type == TRIG) then
        xedge = min(xedge, xpos(2)+xdepth(2)-dx)
        ablX_M(i) = 1.d0/(1.d0 + (g(xedge+dx,1)-g(xedge,1))/dx)
        ablX_R(i) = ablX_M(i)*(1.d0 + gp(xcenter,1))
      else
        ablX_M(i) = 1.d0/(1.d0 + gp(xcenter,1))
        ablX_R(i) = 1.d0
      end if

    end do

    do i=i_abl_upper,mx+num_ghost

      xcenter = xlower + (i-0.5d0)*dx

      ! Calculate inverse of Jacobian at cell edges
      xedge = max(xcenter-0.5d0*dx, xpos(1)-xdepth(1))
      xedge = min(xedge, xpos(2)+xdepth(2))
      ablX_J(i) = 1.d0/(1.d0 + gp(xedge,1))

      ! Calculate grid map scaling and ratio at cell centers
      xcenter = max(xcenter, xpos(1)-xdepth(1)+0.5d0*dx)
      xcenter = min(xcenter, xpos(2)+xdepth(2)-0.5d0*dx)
      if (abl_type == TRIG) then
        xedge = min(xedge, xpos(2)+xdepth(2)-dx)
        ablX_M(i) = 1.d0/(1.d0 + (g(xedge+dx,1)-g(xedge,1))/dx)
        ablX_R(i) = ablX_M(i)*(1.d0 + gp(xcenter,1))
      else
        ablX_M(i) = 1.d0/(1.d0 + gp(xcenter,1))
        ablX_R(i) = 1.d0
      end if

    end do


    ! y direction

    ! determine which indices are in the layer and loop over those
    j_abl_lower = int((1.d0 + 1.d-14)*(ypos(1)-ylower)/dy)
    j_abl_upper = 1 + int((1.d0 + 1.d-14)*(ypos(2)-ylower)/dy)

    do j=1-num_ghost,j_abl_lower

      ycenter = ylower + (j-0.5d0)*dy

      ! Calculate inverse of Jacobian at cell edges
      yedge = max(ycenter-0.5d0*dy, ypos(1)-ydepth(1))
      yedge = min(yedge, ypos(2)+ydepth(2))
      ablY_J(j) = 1.d0/(1.d0 + gp(yedge,2))

      ! Calculate grid map scaling and ratio at cell centers
      ycenter = max(ycenter, ypos(1)-ydepth(1)+0.5d0*dy)
      ycenter = min(ycenter, ypos(2)+ydepth(2)-0.5d0*dy)
      if (abl_type == TRIG) then
        yedge = min(yedge, ypos(2)+ydepth(2)-dy)
        ablY_M(j) = 1.d0/(1.d0 + (g(yedge+dy,2)-g(yedge,2))/dy)
        ablY_R(j) = ablY_M(j)*(1.d0 + gp(ycenter,2))
      else
        ablY_M(j) = 1.d0/(1.d0 + gp(ycenter,2))
        ablY_R(j) = 1.d0
      end if

    end do

    do j=j_abl_upper,my+num_ghost

      ycenter = ylower + (j-0.5d0)*dy

      ! Calculate inverse of Jacobian at cell edges
      yedge = max(ycenter-0.5d0*dy, ypos(1)-ydepth(1))
      yedge = min(yedge, ypos(2)+ydepth(2))
      ablY_J(j) = 1.d0/(1.d0 + gp(yedge,2))

      ! Calculate grid map scaling and ratio at cell centers
      ycenter = max(ycenter, ypos(1)-ydepth(1)+0.5d0*dy)
      ycenter = min(ycenter, ypos(2)+ydepth(2)-0.5d0*dy)
      if (abl_type == TRIG) then
        yedge = min(yedge, ypos(2)+ydepth(2)-dy)
        ablY_M(j) = 1.d0/(1.d0 + (g(yedge+dy,2)-g(yedge,2))/dy)
        ablY_R(j) = ablY_M(j)*(1.d0 + gp(ycenter,2))
      else
        ablY_M(j) = 1.d0/(1.d0 + gp(ycenter,2))
        ablY_R(j) = 1.d0
      end if

    end do


  end subroutine set_scaling_factor

  function gp(z,direction)

    implicit none

    real(kind=8) :: gp, z
    integer :: direction

    real(kind=8) :: pos(2), depth(2), sig, val

    gp = 0.d0
    sig = 0.d0
    if (direction == 1) then
      pos = xpos
      depth = xdepth
    else if (direction == 2) then
      pos = ypos
      depth = ydepth
    end if

    if (z < pos(1)) then
      val = (pos(1)-z)/depth(1)
    else if (z > pos(2)) then
      val = (z-pos(2))/depth(2)
    else
      return
    end if

    if (abl_type == TRIG) then
      gp = dtan(PIH*val)**2
    else if (abl_type == AC) then
      sig = (1.0 - eps_AC)*(1.0 - (1.0 - val)**p_AC)**q_AC
      gp = sig/(1.0-sig)
    else if (abl_type == PS) then
      sig = (1.0 - eps_PS)* &
      val**6*(462.0 - 1980.0*val + 3465.0*val**2 - 3080.0*val**3 + 1386.0*val**4 - 252.0*val**5)
      gp = sig/(1.0-sig)
    end if

  end function gp

  function g(z,direction)

    implicit none

    real(kind=8) :: g, z
    integer :: direction

    real(kind=8) :: pos(2), depth(2), val, depth_scalar

    g = 0.d0
    if (direction == 1) then
      pos = xpos
      depth = xdepth
    else if (direction == 2) then
      pos = ypos
      depth = ydepth
    end if

    if (z < pos(1)) then
      val = (z-pos(1))/depth(1)
      depth_scalar = depth(1)
    else if (z > pos(2)) then
      val = (z-pos(2))/depth(2)
      depth_scalar = depth(2)
    else
      return
    end if

    if (abl_type == TRIG) then
      g = depth_scalar*(dtan(PIH*val)/PIH - val)
    end if

  end function g



end module abl_module
