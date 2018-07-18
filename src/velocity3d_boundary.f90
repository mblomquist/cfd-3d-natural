! velocity3d_boundary Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-16 (YYYY-MM-DD)
!
! This subroutine resets boundary conditions and coefficients for
! u-, v-, or w-velocity.

subroutine velocity3d_boundary(direction)

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define input variables
  character :: direction

  ! Define internal variables
  integer :: i, j, k

  ! ====================== U-Velocity ====================== !
  if (direction .eq. "u") then
    
	! Set Boundary Condition :: West Plane
    Aw_u(1,2:n-2,2:l-2) = 0.
    Ae_u(1,2:n-2,2:l-2) = 0.
    As_u(1,2:n-2,2:l-2) = 0.
    An_u(1,2:n-2,2:l-2) = 0.
    Ab_u(1,2:n-2,2:l-2) = 0.
    At_u(1,2:n-2,2:l-2) = 0.
    Ap_u(1,2:n-2,2:l-2) = 1.
    b_u(1,2:n-2,2:l-2) = 0.

	! Check Boundary Condition Type
	if (u_bc_wc .eq. 1) then
	  ! Given Input Condition
	  b_u(1,2:n-2,2:l-2) = u_bc_wv
	elseif (u_bc_wc .eq. 2) then
	  ! Symmetry Condition
	  Ae_u(1,2:n-2,2:l-2) = 1.
	  b_u(1,2:n-2,2:l-2) = 0.
	else
	  ! Wall Condition
	  b_u(1,2:n-2,2:l-2) = 0.
	end if

    ! Set Boundary Condition :: East Plane
    Aw_u(m,2:n-2,2:l-2) = 0.
    Ae_u(m,2:n-2,2:l-2) = 0.
    As_u(m,2:n-2,2:l-2) = 0.
    An_u(m,2:n-2,2:l-2) = 0.
    Ab_u(m,2:n-2,2:l-2) = 0.
    At_u(m,2:n-2,2:l-2) = 0.
    Ap_u(m,2:n-2,2:l-2) = 1.
    b_u(m,2:n-2,2:l-2) = 0.

	! Check Boundary Condition Type
	if (u_bc_ec .eq. 1) then
	  ! Given Input Condition
	  b_u(m,2:n-2,2:l-2) = u_bc_ev
	elseif (u_bc_ec .eq. 2) then
	  ! Symmetry Condition
	  Aw_u(m,2:n-2,2:l-2) = 1.
	  b_u(m,2:n-2,2:l-2) = 0.
	else
	  ! Wall Condition
	  b_u(m,2:n-2,2:l-2) = 0.
	end if

    ! Set Boundary Condition :: South Plane
    Aw_u(2:m-1,1,2:l-2) = 0.
    Ae_u(2:m-1,1,2:l-2) = 0.
    As_u(2:m-1,1,2:l-2) = 0.
    An_u(2:m-1,1,2:l-2) = 0.
    Ab_u(2:m-1,1,2:l-2) = 0.
    At_u(2:m-1,1,2:l-2) = 0.
    Ap_u(2:m-1,1,2:l-2) = 1.
    b_u(2:m-1,1,2:l-2) = 0.

	! Check Boundary Condition Type
	if (u_bc_sc .eq. 1) then
	  ! Given Input Condition
	  b_u(2:m-1,1,2:l-2) = u_bc_sv
	elseif (u_bc_sc .eq. 2) then
	  ! Symmetry Condition
	  An_u(2:m-1,1,2:l-2) = 1.
	  b_u(2:m-1,1,2:l-2) = 0.
	else
	  ! Wall Condition
	  b_u(2:m-1,1,2:l-2) = 0.
	end if

    ! Set Boundary Condition :: North Plane
    Aw_u(2:m-1,n-1,2:l-2) = 0.
    Ae_u(2:m-1,n-1,2:l-2) = 0.
    As_u(2:m-1,n-1,2:l-2) = 0.
    An_u(2:m-1,n-1,2:l-2) = 0.
    Ab_u(2:m-1,n-1,2:l-2) = 0.
    At_u(2:m-1,n-1,2:l-2) = 0.
    Ap_u(2:m-1,n-1,2:l-2) = 1.
    b_u(2:m-1,n-1,2:l-2) = 0.

	! Check Boundary Condition Type
	if (u_bc_nc .eq. 1) then
	  ! Given Input Condition
	  b_u(2:m-1,n-1,2:l-2) = u_bc_nv
	elseif (u_bc_nc .eq. 2) then
	  ! Symmetry Condition
	  As_u(2:m-1,n-1,2:l-2) = 1.
	  b_u(2:m-1,n-1,2:l-2) = 0.
	else
	  ! Wall Condition
	  b_u(2:m-1,n-1,2:l-2) = 0.
	end if

    ! Set Boundary Condition :: Bottom Plane
    Aw_u(2:m-1,2:n-2,1) = 0.
    Ae_u(2:m-1,2:n-2,1) = 0.
    As_u(2:m-1,2:n-2,1) = 0.
    An_u(2:m-1,2:n-2,1) = 0.
    Ab_u(2:m-1,2:n-2,1) = 0.
    At_u(2:m-1,2:n-2,1) = 0.
    Ap_u(2:m-1,2:n-2,1) = 1.
    b_u(2:m-1,2:n-2,1) = 0.

	! Check Boundary Condition Type
	if (u_bc_bc .eq. 1) then
	  ! Given Input Condition
	  b_u(2:m-1,2:n-2,1) = u_bc_bv
	elseif (u_bc_bc .eq. 2) then
	  ! Symmetry Condition
	  At_u(2:m-1,2:n-2,1) = 1.
	  b_u(2:m-1,2:n-2,1) = 0.
	else
	  ! Wall Condition
	  b_u(2:m-1,2:n-2,1) = 0.
	end if

    ! Set Boundary Condition :: Top Plane
    Aw_u(2:m-1,2:n-2,l-1) = 0.
    Ae_u(2:m-1,2:n-2,l-1) = 0.
    As_u(2:m-1,2:n-2,l-1) = 0.
    An_u(2:m-1,2:n-2,l-1) = 0.
    Ab_u(2:m-1,2:n-2,l-1) = 0.
    At_u(2:m-1,2:n-2,l-1) = 0.
    Ap_u(2:m-1,2:n-2,l-1) = 1.
    b_u(2:m-1,2:n-2,l-1) = 0.

	! Check Boundary Condition Type
	if (u_bc_tc .eq. 1) then
	  ! Given Input Condition
	  b_u(2:m-1,2:n-2,l-1) = u_bc_tv
	elseif (u_bc_tc .eq. 2) then
	  ! Symmetry Condition
	  Ab_u(2:m-1,2:n-2,l-1) = 1.
	  b_u(2:m-1,2:n-2,l-1) = 0.
	else
	  ! Wall Condition
	  b_u(2:m-1,2:n-2,l-1) = 0.
	end if

    ! Set Boundary Condition :: West-South Line
	Aw_u(1,1,2:l-2) = 0.
    Ae_u(1,1,2:l-2) = 1.
    As_u(1,1,2:l-2) = 0.
    An_u(1,1,2:l-2) = 1.
    Ab_u(1,1,2:l-2) = 0.
    At_u(1,1,2:l-2) = 0.
    Ap_u(1,1,2:l-2) = 2.
    b_u(1,1,2:l-2) = 0.

    ! Set Boundary Condition :: West-North Line
    Aw_u(1,n-1,2:l-2) = 0.
    Ae_u(1,n-1,2:l-2) = 1.
    As_u(1,n-1,2:l-2) = 1.
    An_u(1,n-1,2:l-2) = 0.
    Ab_u(1,n-1,2:l-2) = 0.
    At_u(1,n-1,2:l-2) = 0.
    Ap_u(1,n-1,2:l-2) = 2.
    b_u(1,n-1,2:l-2) = 0.

    ! Set Boundary Condition :: West-Bottom Line
	Aw_u(1,2:n-2,1) = 0.
    Ae_u(1,2:n-2,1) = 1.
    As_u(1,2:n-2,1) = 0.
    An_u(1,2:n-2,1) = 0.
    Ab_u(1,2:n-2,1) = 0.
    At_u(1,2:n-2,1) = 1.
    Ap_u(1,2:n-2,1) = 2.
    b_u(1,2:n-2,1) = 0.

    ! Set Boundary Condition :: West-Top Line
	Aw_u(1,2:n-2,l-1) = 0.
    Ae_u(1,2:n-2,l-1) = 1.
    As_u(1,2:n-2,l-1) = 0.
    An_u(1,2:n-2,l-1) = 0.
    Ab_u(1,2:n-2,l-1) = 1.
    At_u(1,2:n-2,l-1) = 0.
    Ap_u(1,2:n-2,l-1) = 2.
    b_u(1,2:n-2,l-1) = 0.

    ! Set Boundary Condition :: East-South Line
	Aw_u(m,1,2:l-2) = 1.
    Ae_u(m,1,2:l-2) = 0.
    As_u(m,1,2:l-2) = 0.
    An_u(m,1,2:l-2) = 1.
    Ab_u(m,1,2:l-2) = 0.
    At_u(m,1,2:l-2) = 0.
    Ap_u(m,1,2:l-2) = 2.
    b_u(m,1,2:l-2) = 0.

    ! Set Boundary Condition :: East-North Line
    Aw_u(m,n-1,2:l-2) = 1.
    Ae_u(m,n-1,2:l-2) = 0.
    As_u(m,n-1,2:l-2) = 1.
    An_u(m,n-1,2:l-2) = 0.
    Ab_u(m,n-1,2:l-2) = 0.
    At_u(m,n-1,2:l-2) = 0.
    Ap_u(m,n-1,2:l-2) = 2.
    b_u(m,n-1,2:l-2) = 0.

    ! Set Boundary Condition :: East-Bottom Line
	Aw_u(m,2:n-2,1) = 1.
    Ae_u(m,2:n-2,1) = 0.
    As_u(m,2:n-2,1) = 0.
    An_u(m,2:n-2,1) = 0.
    Ab_u(m,2:n-2,1) = 0.
    At_u(m,2:n-2,1) = 1.
    Ap_u(m,2:n-2,1) = 2.
    b_u(m,2:n-2,1) = 0.

    ! Set Boundary Condition :: East-Top Line
	Aw_u(m,2:n-2,l-1) = 1.
    Ae_u(m,2:n-2,l-1) = 0.
    As_u(m,2:n-2,l-1) = 0.
    An_u(m,2:n-2,l-1) = 0.
    Ab_u(m,2:n-2,l-1) = 1.
    At_u(m,2:n-2,l-1) = 0.
    Ap_u(m,2:n-2,l-1) = 2.
    b_u(m,2:n-2,l-1) = 0.

    ! Set Boundary Condition :: South-Bottom Line
	Aw_u(2:m-1,1,1) = 0.
    Ae_u(2:m-1,1,1) = 0.
    As_u(2:m-1,1,1) = 0.
    An_u(2:m-1,1,1) = 1.
    Ab_u(2:m-1,1,1) = 0.
    At_u(2:m-1,1,1) = 1.
    Ap_u(2:m-1,1,1) = 2.
    b_u(2:m-1,1,1) = 0.

    ! Set Boundary Condition :: South-Top Line
	Aw_u(2:m-1,1,l-1) = 0.
    Ae_u(2:m-1,1,l-1) = 0.
    As_u(2:m-1,1,l-1) = 0.
    An_u(2:m-1,1,l-1) = 1.
    Ab_u(2:m-1,1,l-1) = 1.
    At_u(2:m-1,1,l-1) = 0.
    Ap_u(2:m-1,1,l-1) = 2.
    b_u(2:m-1,1,l-1) = 0.

    ! Set Boundary Condition :: North-Bottom Line
	Aw_u(2:m-1,n-1,1) = 0.
    Ae_u(2:m-1,n-1,1) = 0.
    As_u(2:m-1,n-1,1) = 1.
    An_u(2:m-1,n-1,1) = 0.
    Ab_u(2:m-1,n-1,1) = 0.
    At_u(2:m-1,n-1,1) = 1.
    Ap_u(2:m-1,n-1,1) = 2.
    b_u(2:m-1,n-1,1) = 0.

    ! Set Boundary Condition :: North-Top Line
	Aw_u(2:m-1,n-1,l-1) = 0.
    Ae_u(2:m-1,n-1,l-1) = 0.
    As_u(2:m-1,n-1,l-1) = 1.
    An_u(2:m-1,n-1,l-1) = 0.
    Ab_u(2:m-1,n-1,l-1) = 1.
    At_u(2:m-1,n-1,l-1) = 0.
    Ap_u(2:m-1,n-1,l-1) = 2.
    b_u(2:m-1,n-1,l-1) = 0.

    ! Set Boundary Condition :: West-South-Bottom Corner
	Aw_u(1,1,1) = 0.
    Ae_u(1,1,1) = 1.
    As_u(1,1,1) = 0.
    An_u(1,1,1) = 1.
    Ab_u(1,1,1) = 0.
    At_u(1,1,1) = 1.
    Ap_u(1,1,1) = 3.
    b_u(1,1,1) = 0.

    ! Set Boundary Condition :: West-North-Bottom Corner
	Aw_u(1,n-1,1) = 0.
    Ae_u(1,n-1,1) = 1.
    As_u(1,n-1,1) = 1.
    An_u(1,n-1,1) = 0.
    Ab_u(1,n-1,1) = 0.
    At_u(1,n-1,1) = 1.
    Ap_u(1,n-1,1) = 3.
    b_u(1,n-1,1) = 0.

    ! Set Boundary Condition :: East-South-Bottom Corner
	Aw_u(m,1,1) = 1.
    Ae_u(m,1,1) = 0.
    As_u(m,1,1) = 0.
    An_u(m,1,1) = 1.
    Ab_u(m,1,1) = 0.
    At_u(m,1,1) = 1.
    Ap_u(m,1,1) = 3.
    b_u(m,1,1) = 0.

    ! Set Boundary Condition :: East-North-Bottom Corner
	Aw_u(m,n-1,1) = 1.
    Ae_u(m,n-1,1) = 0.
    As_u(m,n-1,1) = 1.
    An_u(m,n-1,1) = 0.
    Ab_u(m,n-1,1) = 0.
    At_u(m,n-1,1) = 1.
    Ap_u(m,n-1,1) = 3.
    b_u(m,n-1,1) = 0.

    ! Set Boundary Condition :: West-South-Top Corner
	Aw_u(1,1,l-1) = 0.
    Ae_u(1,1,l-1) = 1.
    As_u(1,1,l-1) = 0.
    An_u(1,1,l-1) = 1.
    Ab_u(1,1,l-1) = 1.
    At_u(1,1,l-1) = 0.
    Ap_u(1,1,l-1) = 3.
    b_u(1,1,l-1) = 0.

    ! Set Boundary Condition :: West-North-Top Corner
	Aw_u(1,n-1,l-1) = 0.
    Ae_u(1,n-1,l-1) = 1.
    As_u(1,n-1,l-1) = 1.
    An_u(1,n-1,l-1) = 0.
    Ab_u(1,n-1,l-1) = 1.
    At_u(1,n-1,l-1) = 0.
    Ap_u(1,n-1,l-1) = 3.
    b_u(1,n-1,l-1) = 0.

    ! Set Boundary Condition :: East-South-Top Corner
	Aw_u(m,1,l-1) = 1.
    Ae_u(m,1,l-1) = 0.
    As_u(m,1,l-1) = 0.
    An_u(m,1,l-1) = 1.
    Ab_u(m,1,l-1) = 1.
    At_u(m,1,l-1) = 0.
    Ap_u(m,1,l-1) = 3.
    b_u(m,1,l-1) = 0.

    ! Set Boundary Condition :: East-North-Top Corner
	Aw_u(m,n-1,l-1) = 1.
    Ae_u(m,n-1,l-1) = 0.
    As_u(m,n-1,l-1) = 1.
    An_u(m,n-1,l-1) = 0.
    Ab_u(m,n-1,l-1) = 1.
    At_u(m,n-1,l-1) = 0.
    Ap_u(m,n-1,l-1) = 3.
    b_u(m,n-1,l-1) = 0.

  end if

  ! ====================== V-Velocity ====================== !
  if (direction .eq. "v") then
    ! Set Boundary Condition :: West Plane
    ! Set Boundary Condition :: East Plane
    ! Set Boundary Condition :: South Plane
    ! Set Boundary Condition :: North Plane
    ! Set Boundary Condition :: Bottom Plane
    ! Set Boundary Condition :: Top Plane

    ! Set Boundary Condition :: West-South Line
    ! Set Boundary Condition :: West-North Line
    ! Set Boundary Condition :: West-Bottom Line
    ! Set Boundary Condition :: West-Top Line
    ! Set Boundary Condition :: East-South Line
    ! Set Boundary Condition :: East-North Line
    ! Set Boundary Condition :: East-Bottom Line
    ! Set Boundary Condition :: East-Top Line
    ! Set Boundary Condition :: South-Bottom Line
    ! Set Boundary Condition :: South-Top Line
    ! Set Boundary Condition :: North-Bottom Line
    ! Set Boundary Condition :: North-Top Line

    ! Set Boundary Condition :: West-South-Bottom Corner
    ! Set Boundary Condition :: West-North-Bottom Corner
    ! Set Boundary Condition :: East-South-Bottom Corner
    ! Set Boundary Condition :: East-North-Bottom Corner
    ! Set Boundary Condition :: West-South-Top Corner
    ! Set Boundary Condition :: West-North-Top Corner
    ! Set Boundary Condition :: East-South-Top Corner
    ! Set Boundary Condition :: East-North-Top Corner

  end if

  ! ====================== W-Velocity ====================== !
  if (direction .eq. "w") then
    ! Set Boundary Condition :: West Plane
    ! Set Boundary Condition :: East Plane
    ! Set Boundary Condition :: South Plane
    ! Set Boundary Condition :: North Plane
    ! Set Boundary Condition :: Bottom Plane
    ! Set Boundary Condition :: Top Plane

    ! Set Boundary Condition :: West-South Line
    ! Set Boundary Condition :: West-North Line
    ! Set Boundary Condition :: West-Bottom Line
    ! Set Boundary Condition :: West-Top Line
    ! Set Boundary Condition :: East-South Line
    ! Set Boundary Condition :: East-North Line
    ! Set Boundary Condition :: East-Bottom Line
    ! Set Boundary Condition :: East-Top Line
    ! Set Boundary Condition :: South-Bottom Line
    ! Set Boundary Condition :: South-Top Line
    ! Set Boundary Condition :: North-Bottom Line
    ! Set Boundary Condition :: North-Top Line

    ! Set Boundary Condition :: West-South-Bottom Corner
    ! Set Boundary Condition :: West-North-Bottom Corner
    ! Set Boundary Condition :: East-South-Bottom Corner
    ! Set Boundary Condition :: East-North-Bottom Corner
    ! Set Boundary Condition :: West-South-Top Corner
    ! Set Boundary Condition :: West-North-Top Corner
    ! Set Boundary Condition :: East-South-Top Corner
    ! Set Boundary Condition :: East-North-Top Corner

  end if

  return

  end subroutine velocity3d_boundary