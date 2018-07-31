! temperature3d Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms, sets boundary values, and initializes the temperature grid.
!
! Boundary Conditions ::
!   0 - Wall Type Condition
!   1 - Fixed Value Type Condition
!   2 - Symmetry Type Condition

subroutine temperature3d

  ! Pull in standard variable header
  include "var3d.dec"

  ! Initialize Temperature Field
  T = 0.

  ! Initialize coefficients
  Aw_T = 0.
  Ae_T = 0.
  As_T = 0.
  An_T = 0.
  Ab_T = 0.
  At_T = 0.
  Ap_T = 1.
  b_T = 0.

  ! Set temperature source terms
  Su_T = 0.
  Sp_T = 0.


  ! ====================================== !
  ! Set boundary conditions :: West Plane
  if (T_bc_wc .eq. 0) then
    Sp_T(1,2:n-2,2:l-2) = 0.
    Su_T(1,2:n-2,2:l-2) = 0.
  elseif (T_bc_wc .eq. 1) then
    Sp_T(1,2:n-2,2:l-2) = 0.
    Su_T(1,2:n-2,2:l-2) = (T_bc_wv-T_c)/(T_h-T_c)
  else
    Ae_T(1,2:n-2,2:l-2) = 1.
	Sp_T(1,2:n-2,2:l-2) = 0.
    Su_T(1,2:n-2,2:l-2) = 0.
  end if

  ! Set boundary conditions :: East Plane
  if (T_bc_ec .eq. 0) then
    Sp_T(m-1,2:n-2,2:l-2) = 0.
    Su_T(m-1,2:n-2,2:l-2) = 0.
  elseif (T_bc_ec .eq. 1) then
    Sp_T(m-1,2:n-2,2:l-2) = 0.
    Su_T(m-1,2:n-2,2:l-2) = (T_bc_ev-T_c)/(T_h-T_c)
  else
    Aw_T(m-1,2:n-2,2:l-2) = 1.
	Sp_T(m-1,2:n-2,2:l-2) = 0.
    Su_T(m-1,2:n-2,2:l-2) = 0.
  end if

  ! Set boundary conditions :: South Plane
  if (T_bc_sc .eq. 0) then
    Sp_T(2:m-2,1,2:l-2) = 0.
    Su_T(2:m-2,1,2:l-2) = 0.
  elseif (T_bc_sc .eq. 1) then
    Sp_T(2:m-2,1,2:l-2) = 0.
    Su_T(2:m-2,1,2:l-2) = (T_bc_sv-T_c)/(T_h-T_c)
  else
    An_T(2:m-2,1,2:l-2) = 1.
	Sp_T(2:m-2,1,2:l-2) = 0.
    Su_T(2:m-2,1,2:l-2) = 0.
  end if

  ! Set boundary conditions :: North Plane
  if (T_bc_nc .eq. 0) then
    Sp_T(2:m-2,n-1,2:l-2) = 0.
    Su_T(2:m-2,n-1,2:l-2) = 0.
  elseif (T_bc_nc .eq. 1) then
    Sp_T(2:m-2,n-1,2:l-2) = 0.
    Su_T(2:m-2,n-1,2:l-2) = (T_bc_nv-T_c)/(T_h-T_c)
  else
    As_T(2:m-2,n-1,2:l-2) = 1.
	Sp_T(2:m-2,n-1,2:l-2) = 0.
    Su_T(2:m-2,n-1,2:l-2) = 0.
  end if

  ! Set boundary conditions :: Bottom Plane
  if (T_bc_bc .eq. 0) then
    Sp_T(2:m-2,2:n-2,1) = 0.
    Su_T(2:m-2,2:n-2,1) = 0.
  elseif (T_bc_bc .eq. 1) then
    Sp_T(2:m-2,2:n-2,1) = 0.
    Su_T(2:m-2,2:n-2,1) = (T_bc_bv-T_c)/(T_h-T_c)
  else
    At_T(2:m-2,2:n-2,1) = 1.
	Sp_T(2:m-2,2:n-2,1) = 0.
    Su_T(2:m-2,2:n-2,1) = 0.
  end if

  ! Set boundary conditions :: Top Plane
  if (T_bc_tc .eq. 0) then
    Sp_T(2:m-2,2:n-2,l-1) = 0.
    Su_T(2:m-2,2:n-2,l-1) = 0.
  elseif (T_bc_tc .eq. 1) then
    Sp_T(2:m-2,2:n-2,l-1) = 0.
    Su_T(2:m-2,2:n-2,l-1) = (T_bc_tv-T_c)/(T_h-T_c)
  else
    Ab_T(2:m-2,2:n-2,l-1) = 1.
	Sp_T(2:m-2,2:n-2,l-1) = 0.
    Su_T(2:m-2,2:n-2,l-1) = 0.
  end if

  ! ====================================== !

  ! ====================================== !
  ! Set boundary conditions :: West-South Line
  Ae_T(1,1,2:l-2) = 1.
  An_T(1,1,2:l-2) = 1.
  Ap_T(1,1,2:l-2) = 2.
  
  ! Set boundary conditions :: West-North Line
  As_T(1,n-1,2:l-2) = 1.
  Ae_T(1,n-1,2:l-2) = 1.
  Ap_T(1,n-1,2:l-2) = 2.
  
  ! Set boundary conditions :: West-Bottom Line
  Ae_T(1,2:n-2,1) = 1.
  At_T(1,2:n-2,1) = 1.
  Ap_T(1,2:n-2,1) = 2.
  
  ! Set boundary conditions :: West-Top Line
  Ab_T(1,2:n-2,l-1) = 1.
  Ae_T(1,2:n-2,l-1) = 1.
  Ap_T(1,2:n-2,l-1) = 2.
  
  ! Set boundary conditions :: East-South Line
  Aw_T(m-1,1,2:l-2) = 1.
  An_T(m-1,1,2:l-2) = 1.
  Ap_T(m-1,1,2:l-2) = 2.
  
  ! Set boundary conditions :: East-North Line
  As_T(m-1,n-1,2:l-2) = 1.
  Aw_T(m-1,n-1,2:l-2) = 1.
  Ap_T(m-1,n-1,2:l-2) = 2.
  
  ! Set boundary conditions :: East-Bottom Line
  Aw_T(m-1,2:n-2,1) = 1.
  At_T(m-1,2:n-2,1) = 1.
  Ap_T(m-1,2:n-2,1) = 2.
  
  ! Set boundary conditions :: East-Top Line
  Ab_T(m-1,2:n-2,l-1) = 1.
  Aw_T(m-1,2:n-2,l-1) = 1.
  Ap_T(m-1,2:n-2,l-1) = 2.
  
  ! Set boundary conditions :: South-Bottom Line
  An_T(2:m-2,1,1) = 1.
  At_T(2:m-2,1,1) = 1.
  Ap_T(2:m-2,1,1) = 2.
  
  ! Set boundary conditions :: South-Top Line
  Ab_T(2:m-2,1,l-1) = 1.
  An_T(2:m-2,1,l-1) = 1.
  Ap_T(2:m-2,1,l-1) = 2.
  
  ! Set boundary conditions :: North-Bottom Line
  As_T(2:m-2,n-1,1) = 1.
  At_T(2:m-2,n-1,1) = 1.
  Ap_T(2:m-2,n-1,1) = 2.
  
  ! Set boundary conditions :: North-Top Line
  Ab_T(2:m-2,n-1,l-1) = 1.
  As_T(2:m-2,n-1,l-1) = 1.
  Ap_T(2:m-2,n-1,l-1) = 2.
  
  ! ====================================== !

  ! ====================================== !
  ! Set boundary conditions :: West-South-Bottom Corner
  Ae_T(1,1,1) = 1.
  An_T(1,1,1) = 1.
  At_T(1,1,1) = 1.
  Ap_T(1,1,1) = 3.
  
  ! Set boundary conditions :: West-North-Bottom Corner
  As_T(1,n-1,1) = 1.
  Ae_T(1,n-1,1) = 1.
  At_T(1,n-1,1) = 1.
  Ap_T(1,n-1,1) = 3.
  
  ! Set boundary conditions :: East-South-Bottom Corner
  Aw_T(m-1,1,1) = 1.
  An_T(m-1,1,1) = 1.
  At_T(m-1,1,1) = 1.
  Ap_T(m-1,1,1) = 3.
  
  ! Set boundary conditions :: East-North-Bottom Corner
  As_T(m-1,n-1,1) = 1.
  Aw_T(m-1,n-1,1) = 1.
  At_T(m-1,n-1,1) = 1.
  Ap_T(m-1,n-1,1) = 3.
  
  ! Set boundary conditions :: West-South-Top Corner
  Ab_T(1,1,l-1) = 1.
  Ae_T(1,1,l-1) = 1.
  An_T(1,1,l-1) = 1.
  Ap_T(1,1,l-1) = 3.
  
  ! Set boundary conditions :: West-North-Top Corner
  Ab_T(1,n-1,l-1) = 1.
  As_T(1,n-1,l-1) = 1.
  Ae_T(1,n-1,l-1) = 1.
  Ap_T(1,n-1,l-1) = 3.
  
  ! Set boundary conditions :: East-South-Top Corner
  Ab_T(m-1,1,l-1) = 1.
  Aw_T(m-1,1,l-1) = 1.
  An_T(m-1,1,l-1) = 1.
  Ap_T(m-1,1,l-1) = 3.
  
  ! Set boundary conditions :: East-North-Top Corner
  Ab_T(m-1,n-1,l-1) = 1.
  As_T(m-1,n-1,l-1) = 1.
  Aw_T(m-1,n-1,l-1) = 1.
  Ap_T(m-1,n-1,l-1) = 3.
  
  return

end subroutine temperature3d
