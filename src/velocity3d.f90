! velocity3d Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-16 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms, sets boundary values, and initializes the velocity grids.
!
! Boundary Conditions ::
!   0 - Wall Type Condition
!   1 - Fixed Value Type Condition
!   2 - Symmetry Type Condition

subroutine velocity3d

  ! Pull in standard variable header
  include "var3d.dec"

  ! Initalize Velocity Fields
  u = 0.
  v = 0.
  w = 0.

  ! ====================================== !
  ! Initialize coefficients :: u
  Ab_u = 0.
  As_u = 0.
  Aw_u = 0.
  Ae_u = 0.
  An_u = 0.
  At_u = 0.
  Ap_u = 1.
  b_u = 0.

  ! Initalize source terms :: u
  Sp_u = 0.
  Su_u = 0.

  ! Initialize coefficients :: v
  Ab_v = 0.
  As_v = 0.
  Aw_v = 0.
  Ae_v = 0.
  An_v = 0.
  At_v = 0.
  Ap_v = 1.
  b_v = 0.

  ! Initalize source terms :: v
  Sp_v = 0.
  Su_v = 0.

  ! Initialize coefficients :: w
  Ab_w = 0.
  As_w = 0.
  Aw_w = 0.
  Ae_w = 0.
  An_w = 0.
  At_w = 0.
  Ap_w = 1.
  b_w = 0.

  ! Initalize source terms :: w
  Sp_w = 0.
  Su_w = 0.

  ! ====================================== !
  call velocity3d_boundary

  return

end subroutine velocity3d

subroutine velocity3d_boundary

  ! Pull in standard variable header
  include "var3d.dec"

  ! Set Ap terms to 1.
  b_u = 0.
  b_v = 0.
  b_w = 0.
  
  Ap_u = 1.
  Ap_v = 1.
  Ap_w = 1.

  ! ====================================== !
  ! Set boundary conditions :: West Plane
  if (u_bc_wc .eq. 0) then
	Sp_u(1,2:n-2,2:l-2) = 0.
	Su_u(1,2:n-2,2:l-2) = 0.
  elseif (u_bc_wc .eq. 1) then
	Sp_u(1,2:n-2,2:l-2) = 0.
	Su_u(1,2:n-2,2:l-2) = u_bc_wv
  else
    Ae_u(1,2:n-2,2:l-2) = 1.
	Sp_u(1,2:n-2,2:l-2) = 0.
	Su_u(1,2:n-2,2:l-2) = 0.
  end if

  if (v_bc_wc .eq. 0) then
    Sp_v(1,2:n-1,2:l-2) = 2.*dx*dy/dz/Ra**(0.5)
	Su_v(1,2:n-1,2:l-2) = 0.
  elseif (v_bc_wc .eq. 1) then
    Sp_v(1,2:n-1,2:l-2) = 0.
	Su_v(1,2:n-1,2:l-2) = v_bc_wv
  else
    Ae_v(1,2:n-1,2:l-2) = 1.
    Sp_v(1,2:n-1,2:l-2) = 0.
	Su_v(1,2:n-1,2:l-2) = 0.
  end if

  if (w_bc_wc .eq. 0) then
    Sp_w(1,2:n-2,2:l-1) = 2.*dy*dz/dx/Ra**(0.5)
	Su_w(1,2:n-2,2:l-1) = 0.
  elseif (w_bc_wc .eq. 1) then
    Sp_w(1,2:n-2,2:l-1) = 0.
	Su_w(1,2:n-2,2:l-1) = w_bc_wv
  else
    Ae_w(1,2:n-2,2:l-1) = 1.
    Sp_w(1,2:n-2,2:l-1) = 0.
	Su_w(1,2:n-2,2:l-1) = 0.
  end if

  ! Set boundary conditions :: East Plane
  if (u_bc_ec .eq. 0) then
	Sp_u(m,2:n-2,2:l-2) = 0.
	Su_u(m,2:n-2,2:l-2) = 0.
  elseif (u_bc_ec .eq. 1) then
	Sp_u(m,2:n-2,2:l-2) = 0.
	Su_u(m,2:n-2,2:l-2) = u_bc_ev
  else
    Aw_u(m,2:n-2,2:l-2) = 1.
	Sp_u(m,2:n-2,2:l-2) = 0.
	Su_u(m,2:n-2,2:l-2) = 0.
  end if

  if (v_bc_ec .eq. 0) then
    Sp_v(m-1,2:n-1,2:l-2) = 2.*dx*dy/dz/Ra**(0.5)
	Su_v(m-1,2:n-1,2:l-2) = 0.
  elseif (v_bc_ec .eq. 1) then
    Sp_v(m-1,2:n-1,2:l-2) = 0.
	Su_v(m-1,2:n-1,2:l-2) = v_bc_ev
  else
    Aw_v(m-1,2:n-1,2:l-2) = 1.
    Sp_v(m-1,2:n-1,2:l-2) = 0.
	Su_v(m-1,2:n-1,2:l-2) = 0.
  end if

  if (w_bc_ec .eq. 0) then
    Sp_w(m-1,2:n-2,2:l-1) = 2.*dy*dz/dx/Ra**(0.5)
	Su_w(m-1,2:n-2,2:l-1) = 0.
  elseif (w_bc_ec .eq. 1) then
    Sp_w(m-1,2:n-2,2:l-1) = 0.
	Su_w(m-1,2:n-2,2:l-1) = w_bc_ev
  else
    Aw_w(m-1,2:n-2,2:l-1) = 1.
    Sp_w(m-1,2:n-2,2:l-1) = 0.
	Su_w(m-1,2:n-2,2:l-1) = 0.
  end if

  ! Set boundary conditions :: South Plane
  if (u_bc_sc .eq. 0) then
	Sp_u(2:m-1,1,2:l-2) = 2.*dz*dx/dy/Ra**(0.5)
	Su_u(2:m-1,1,2:l-2) = 0.
  elseif (u_bc_sc .eq. 1) then
	Sp_u(2:m-1,1,2:l-2) = 0.
	Su_u(2:m-1,1,2:l-2) = u_bc_sv
  else
    An_u(2:m-1,1,2:l-2) = 1.
	Sp_u(2:m-1,1,2:l-2) = 0.
	Su_u(2:m-1,1,2:l-2) = 0.
  end if

  if (v_bc_sc .eq. 0) then
    Sp_v(2:m-2,1,2:l-2) = 0.
	Su_v(2:m-2,1,2:l-2) = 0.
  elseif (v_bc_sc .eq. 1) then
    Sp_v(2:m-2,1,2:l-2) = 0.
	Su_v(2:m-2,1,2:l-2) = v_bc_sv
  else
    An_v(2:m-2,1,2:l-2) = 1.
    Sp_v(2:m-2,1,2:l-2) = 0.
	Su_v(2:m-2,1,2:l-2) = 0.
  end if

  if (w_bc_sc .eq. 0) then
    Sp_w(2:m-2,1,2:l-1) = 2.*dx*dy/dz/Ra**(0.5)
	Su_w(2:m-2,1,2:l-1) = 0.
  elseif (w_bc_sc .eq. 1) then
    Sp_w(2:m-2,1,2:l-1) = 0.
	Su_w(2:m-2,1,2:l-1) = w_bc_sv
  else
    An_w(2:m-2,1,2:l-1) = 1.
    Sp_w(2:m-2,1,2:l-1) = 0.
	Su_w(2:m-2,1,2:l-1) = 0.
  end if

  ! Set boundary conditions :: North Plane
  if (u_bc_nc .eq. 0) then
	Sp_u(2:m-1,n-1,2:l-2) = 2.*dz*dx/dy/Ra**(0.5)
	Su_u(2:m-1,n-1,2:l-2) = 0.
  elseif (u_bc_nc .eq. 1) then
	Sp_u(2:m-1,n-1,2:l-2) = 0.
	Su_u(2:m-1,n-1,2:l-2) = u_bc_nv
  else
    As_u(2:m-1,n-1,2:l-2) = 1.
	Sp_u(2:m-1,n-1,2:l-2) = 0.
	Su_u(2:m-1,n-1,2:l-2) = 0.
  end if

  if (v_bc_nc .eq. 0) then
    Sp_v(2:m-2,n-1,2:l-2) = 0.
	Su_v(2:m-2,n-1,2:l-2) = 0.
  elseif (v_bc_nc .eq. 1) then
    Sp_v(2:m-2,n-1,2:l-2) = 0.
	Su_v(2:m-2,n-1,2:l-2) = v_bc_nv
  else
    As_v(2:m-2,n-1,2:l-2) = 1.
    Sp_v(2:m-2,n-1,2:l-2) = 0.
	Su_v(2:m-2,n-1,2:l-2) = 0.
  end if

  if (w_bc_nc .eq. 0) then
    Sp_w(2:m-2,n-1,2:l-1) = 2.*dx*dy/dz/Ra**(0.5)
	Su_w(2:m-2,n-1,2:l-1) = 0.
  elseif (w_bc_nc .eq. 1) then
    Sp_w(2:m-2,n-1,2:l-1) = 0.
	Su_w(2:m-2,n-1,2:l-1) = w_bc_nv
  else
    As_w(2:m-2,n-1,2:l-1) = 1.
    Sp_w(2:m-2,n-1,2:l-1) = 0.
	Su_w(2:m-2,n-1,2:l-1) = 0.
  end if

  ! Set boundary conditions :: Bottom Plane
  if (u_bc_bc .eq. 0) then
	Sp_u(2:m-1,2:n-2,1) = 2.*dz*dx/dy/Ra**(0.5)
	Su_u(2:m-1,2:n-2,1) = 0.
  elseif (u_bc_bc .eq. 1) then
	Sp_u(2:m-1,2:n-2,1) = 0.
	Su_u(2:m-1,2:n-2,1) = u_bc_bv
  else
    At_u(2:m-1,2:n-2,1) = 1.
	  Sp_u(2:m-1,2:n-2,1) = 0.
	  Su_u(2:m-1,2:n-2,1) = 0.
  end if

  if (v_bc_bc .eq. 0) then
    Sp_v(2:m-2,2:n-1,1) = 2.*dx*dy/dz/Ra**(0.5)
	Su_v(2:m-2,2:n-1,1) = 0.
  elseif (v_bc_bc .eq. 1) then
    Sp_v(2:m-2,2:n-1,1) = 0.
	Su_v(2:m-2,2:n-1,1) = v_bc_bv
  else
    At_v(2:m-2,2:n-1,1) = 1.
    Sp_v(2:m-2,2:n-1,1) = 0.
	Su_v(2:m-2,2:n-1,1) = 0.
  end if

  if (w_bc_bc .eq. 0) then
    Sp_w(2:m-2,2:n-2,1) = 0.
	  Su_w(2:m-2,2:n-2,1) = 0.
  elseif (w_bc_bc .eq. 1) then
    Sp_w(2:m-2,2:n-2,1) = 0.
	  Su_w(2:m-2,2:n-2,1) = w_bc_bv
  else
    At_w(2:m-2,2:n-2,1) = 1.
    Sp_w(2:m-2,2:n-2,1) = 0.
	  Su_w(2:m-2,2:n-2,1) = 0.
  end if

  ! Set boundary conditions :: Top Plane
  if (u_bc_tc .eq. 0) then
	Sp_u(2:m-1,2:n-2,l-1) = 2.*dz*dx/dy/Ra**(0.5)
	Su_u(2:m-1,2:n-2,l-1) = 0.
  elseif (u_bc_tc .eq. 1) then
	Sp_u(2:m-1,2:n-2,l-1) = 0.
	Su_u(2:m-1,2:n-2,1) = u_bc_tv
  else
    Ab_u(2:m-1,2:n-2,l-1) = 1.
	Sp_u(2:m-1,2:n-2,l-1) = 0.
	Su_u(2:m-1,2:n-2,l-1) = 0.
  end if

  if (v_bc_tc .eq. 0) then
    Sp_v(2:m-2,2:n-1,l-1) = 2.*dx*dy/dz/Ra**(0.5)
	  Su_v(2:m-2,2:n-1,l-1) = 0.
  elseif (v_bc_tc .eq. 1) then
    Sp_v(2:m-2,2:n-1,l-1) = 0.
	  Su_v(2:m-2,2:n-1,l-1) = v_bc_tv
  else
    Ab_v(2:m-2,2:n-1,l-1) = 1.
    Sp_v(2:m-2,2:n-1,l-1) = 0.
	  Su_v(2:m-2,2:n-1,l-1) = 0.
  end if

  if (w_bc_tc .eq. 0) then
    Sp_w(2:m-2,2:n-2,l) = 0.
	  Su_w(2:m-2,2:n-2,l) = 0.
  elseif (w_bc_tc .eq. 1) then
    Sp_w(2:m-2,2:n-2,l) = 0.
	  Su_w(2:m-2,2:n-2,l) = w_bc_tv
  else
    Ab_w(2:m-2,2:n-2,l) = 1.
    Sp_w(2:m-2,2:n-2,l) = 0.
	  Su_w(2:m-2,2:n-2,l) = 0.
  end if

  ! ====================================== !

  ! ====================================== !
  ! Set boundary conditions :: West-South Line
  Ae_u(1,1,2:l-2) = 1.
  An_u(1,1,2:l-2) = 1.
  Ap_u(1,1,2:l-2) = 2.

  Ae_v(1,1,2:l-2) = 1.
  An_v(1,1,2:l-2) = 1.
  Ap_v(1,1,2:l-2) = 2.

  Ae_w(1,1,2:l-1) = 1.
  An_w(1,1,2:l-1) = 1.
  Ap_w(1,1,2:l-1) = 2.

  ! Set boundary conditions :: West-North Line
  As_u(1,n-1,2:l-2) = 1.
  Ae_u(1,n-1,2:l-2) = 1.
  Ap_u(1,n-1,2:l-2) = 2.

  As_v(1,n,2:l-2) = 1.
  Ae_v(1,n,2:l-2) = 1.
  Ap_v(1,n,2:l-2) = 2.

  As_w(1,n-1,2:l-1) = 1.
  Ae_w(1,n-1,2:l-1) = 1.
  Ap_w(1,n-1,2:l-1) = 2.

  ! Set boundary conditions :: West-Bottom Line
  Ae_u(1,2:n-2,1) = 1.
  At_u(1,2:n-2,1) = 1.
  Ap_u(1,2:n-2,1) = 2.

  Ae_v(1,2:n-1,1) = 1.
  At_v(1,2:n-1,1) = 1.
  Ap_v(1,2:n-1,1) = 2.

  Ae_w(1,2:n-2,1) = 1.
  At_w(1,2:n-2,1) = 1.
  Ap_w(1,2:n-2,1) = 2.

  ! Set boundary conditions :: West-Top Line
  Ab_u(1,2:n-2,l-1) = 1.
  Ae_u(1,2:n-2,l-1) = 1.
  Ap_u(1,2:n-2,l-1) = 2.

  Ab_v(1,2:n-1,l-1) = 1.
  Ae_v(1,2:n-1,l-1) = 1.
  Ap_v(1,2:n-1,l-1) = 2.

  Ab_w(1,2:n-2,l) = 1.
  Ae_w(1,2:n-2,l) = 1.
  Ap_w(1,2:n-2,l) = 2.

  ! Set boundary conditions :: East-South Line
  Aw_u(m,1,2:l-2) = 1.
  An_u(m,1,2:l-2) = 1.
  Ap_u(m,1,2:l-2) = 2.

  Aw_v(m-1,1,2:l-2) = 1.
  An_v(m-1,1,2:l-2) = 1.
  Ap_v(m-1,1,2:l-2) = 2.

  Aw_w(m-1,1,2:l-1) = 1.
  An_w(m-1,1,2:l-1) = 1.
  Ap_w(m-1,1,2:l-1) = 2.

  ! Set boundary conditions :: East-North Line
  As_u(m,n-1,2:l-2) = 1.
  Aw_u(m,n-1,2:l-2) = 1.
  Ap_u(m,n-1,2:l-2) = 2.

  As_v(m-1,n,2:l-2) = 1.
  Aw_v(m-1,n,2:l-2) = 1.
  Ap_v(m-1,n,2:l-2) = 2.

  As_w(m-1,n-1,2:l-1) = 1.
  Aw_w(m-1,n-1,2:l-1) = 1.
  Ap_w(m-1,n-1,2:l-1) = 2.

  ! Set boundary conditions :: East-Bottom Line
  Aw_u(m,2:n-2,1) = 1.
  At_u(m,2:n-2,1) = 1.
  Ap_u(m,2:n-2,1) = 2.

  Aw_v(m-1,2:n-1,1) = 1.
  At_v(m-1,2:n-1,1) = 1.
  Ap_v(m-1,2:n-1,1) = 2.

  Aw_w(m-1,2:n-2,1) = 1.
  At_w(m-1,2:n-2,1) = 1.
  Ap_w(m-1,2:n-2,1) = 2.

  ! Set boundary conditions :: East-Top Line
  Ab_u(m,2:n-2,l-1) = 1.
  Aw_u(m,2:n-2,l-1) = 1.
  Ap_u(m,2:n-2,l-1) = 2.

  Ab_v(m-1,2:n-1,l-1) = 1.
  Aw_v(m-1,2:n-1,l-1) = 1.
  Ap_v(m-1,2:n-1,l-1) = 2.

  Ab_w(m-1,2:n-2,l) = 1.
  Aw_w(m-1,2:n-2,l) = 1.
  Ap_w(m-1,2:n-2,l) = 2.

  ! Set boundary conditions :: South-Bottom Line
  An_u(2:m-1,1,1) = 1.
  At_u(2:m-1,1,1) = 1.
  Ap_u(2:m-1,1,1) = 2.

  An_v(2:m-2,1,1) = 1.
  At_v(2:m-2,1,1) = 1.
  Ap_v(2:m-2,1,1) = 2.

  An_w(2:m-2,1,1) = 1.
  At_w(2:m-2,1,1) = 1.
  Ap_w(2:m-2,1,1) = 2.

  ! Set boundary conditions :: South-Top Line
  Ab_u(2:m-1,1,l-1) = 1.
  An_u(2:m-1,1,l-1) = 1.
  Ap_u(2:m-1,1,l-1) = 2.

  Ab_v(2:m-2,1,l-1) = 1.
  An_v(2:m-2,1,l-1) = 1.
  Ap_v(2:m-2,1,l-1) = 2.

  Ab_w(2:m-2,1,l) = 1.
  An_w(2:m-2,1,l) = 1.
  Ap_w(2:m-2,1,l) = 2.

  ! Set boundary conditions :: North-Bottom Line
  As_u(2:m-1,n-1,1) = 1.
  At_u(2:m-1,n-1,1) = 1.
  Ap_u(2:m-1,n-1,1) = 2.

  As_v(2:m-2,n,1) = 1.
  At_v(2:m-2,n,1) = 1.
  Ap_v(2:m-2,n,1) = 2.

  As_w(2:m-2,n-1,1) = 1.
  At_w(2:m-2,n-1,1) = 1.
  Ap_w(2:m-2,n-1,1) = 2.

  ! Set boundary conditions :: North-Top Line
  Ab_u(2:m-1,n-1,l-1) = 1.
  As_u(2:m-1,n-1,l-1) = 1.
  Ap_u(2:m-1,n-1,l-1) = 2.

  Ab_v(2:m-2,n,l-1) = 1.
  As_v(2:m-2,n,l-1) = 1.
  Ap_v(2:m-2,n,l-1) = 2.

  Ab_w(2:m-2,n-1,l) = 1.
  As_w(2:m-2,n-1,l) = 1.
  Ap_w(2:m-2,n-1,l) = 2.

  ! ====================================== !

  ! ====================================== !
  ! Set boundary conditions :: West-South-Bottom Corner
  Ae_u(1,1,1) = 1.
  An_u(1,1,1) = 1.
  At_u(1,1,1) = 1.
  Ap_u(1,1,1) = 3.

  Ae_v(1,1,1) = 1.
  An_v(1,1,1) = 1.
  At_v(1,1,1) = 1.
  Ap_v(1,1,1) = 3.

  Ae_w(1,1,1) = 1.
  An_w(1,1,1) = 1.
  At_w(1,1,1) = 1.
  Ap_w(1,1,1) = 3.

  ! Set boundary conditions :: West-North-Bottom Corner
  As_u(1,n-1,1) = 1.
  Ae_u(1,n-1,1) = 1.
  At_u(1,n-1,1) = 1.
  Ap_u(1,n-1,1) = 3.

  As_v(1,n,1) = 1.
  Ae_v(1,n,1) = 1.
  At_v(1,n,1) = 1.
  Ap_v(1,n,1) = 3.

  As_w(1,n-1,1) = 1.
  Ae_w(1,n-1,1) = 1.
  At_w(1,n-1,1) = 1.
  Ap_w(1,n-1,1) = 3.

  ! Set boundary conditions :: East-South-Bottom Corner
  Aw_u(m,1,1) = 1.
  An_u(m,1,1) = 1.
  At_u(m,1,1) = 1.
  Ap_u(m,1,1) = 3.

  Aw_v(m-1,1,1) = 1.
  An_v(m-1,1,1) = 1.
  At_v(m-1,1,1) = 1.
  Ap_v(m-1,1,1) = 3.

  Aw_w(m-1,1,1) = 1.
  An_w(m-1,1,1) = 1.
  At_w(m-1,1,1) = 1.
  Ap_w(m-1,1,1) = 3.

  ! Set boundary conditions :: East-North-Bottom Corner
  As_u(m,n-1,1) = 1.
  Aw_u(m,n-1,1) = 1.
  At_u(m,n-1,1) = 1.
  Ap_u(m,n-1,1) = 3.

  As_v(m-1,n,1) = 1.
  Aw_v(m-1,n,1) = 1.
  At_v(m-1,n,1) = 1.
  Ap_v(m-1,n,1) = 3.

  As_w(m-1,n-1,1) = 1.
  Aw_w(m-1,n-1,1) = 1.
  At_w(m-1,n-1,1) = 1.
  Ap_w(m-1,n-1,1) = 3.

  ! Set boundary conditions :: West-South-Top Corner
  Ab_u(1,1,l-1) = 1.
  Ae_u(1,1,l-1) = 1.
  An_u(1,1,l-1) = 1.
  Ap_u(1,1,l-1) = 3.

  Ab_v(1,1,l-1) = 1.
  Ae_v(1,1,l-1) = 1.
  An_v(1,1,l-1) = 1.
  Ap_v(1,1,l-1) = 3.

  Ab_w(1,1,l) = 1.
  Ae_w(1,1,l) = 1.
  An_w(1,1,l) = 1.
  Ap_w(1,1,l) = 3.

  ! Set boundary conditions :: West-North-Top Corner
  Ab_u(1,n-1,l-1) = 1.
  As_u(1,n-1,l-1) = 1.
  Ae_u(1,n-1,l-1) = 1.
  Ap_u(1,n-1,l-1) = 3.

  Ab_v(1,n,l-1) = 1.
  As_v(1,n,l-1) = 1.
  Ae_v(1,n,l-1) = 1.
  Ap_v(1,n,l-1) = 3.

  Ab_w(1,n-1,l) = 1.
  As_w(1,n-1,l) = 1.
  Ae_w(1,n-1,l) = 1.
  Ap_w(1,n-1,l) = 3.

  ! Set boundary conditions :: East-South-Top Corner
  Ab_u(m,1,l-1) = 1.
  Aw_u(m,1,l-1) = 1.
  An_u(m,1,l-1) = 1.
  Ap_u(m,1,l-1) = 3.

  Ab_v(m-1,1,l-1) = 1.
  Aw_v(m-1,1,l-1) = 1.
  An_v(m-1,1,l-1) = 1.
  Ap_v(m-1,1,l-1) = 3.

  Ab_w(m-1,1,l) = 1.
  Aw_w(m-1,1,l) = 1.
  An_w(m-1,1,l) = 1.
  Ap_w(m-1,1,l) = 3.

  ! Set boundary conditions :: East-North-Top Corner
  Ab_u(m,n-1,l-1) = 1.
  As_u(m,n-1,l-1) = 1.
  Aw_u(m,n-1,l-1) = 1.
  Ap_u(m,n-1,l-1) = 3.

  Ab_v(m-1,n,l-1) = 1.
  As_v(m-1,n,l-1) = 1.
  Aw_v(m-1,n,l-1) = 1.
  Ap_v(m-1,n,l-1) = 3.

  Ab_w(m-1,n-1,l) = 1.
  As_w(m-1,n-1,l) = 1.
  Aw_w(m-1,n-1,l) = 1.
  Ap_w(m-1,n-1,l) = 3.

  ! ====================================== !

  return

end subroutine velocity3d_boundary
