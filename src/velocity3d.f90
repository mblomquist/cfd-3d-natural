! velocity3d Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-16 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms, sets boundary values, and initializes the velocity grids.

subroutine velocity3d

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! ====================== U-Velocity ====================== !
  ! Initialize field :: Set to 0
  u = 0.
  u_hat = 0.
  u_star = 0.

  ! Initialize coefficients
  Aw_u = 0.
  Ae_u = 0.
  As_u = 0.
  An_u = 0.
  Ab_u = 0.
  At_u = 0.
  Ap_u = 1.
  b_u = 0.

  ! Initialize source terms
  Su_u = 0.
  Sp_u = 0.

  call velocity3d_boundary("u")

  ! ====================== V-Velocity ====================== !
  ! Initialize field :: Set to 0
  v = 0.
  v_hat = 0.
  v_star = 0.

  ! Initialize coefficients
  Aw_v = 0.
  Ae_v = 0.
  As_v = 0.
  An_v = 0.
  Ab_v = 0.
  At_v = 0.
  Ap_v = 1.
  b_v = 0.

  ! Initialize source terms
  Su_v = 0.
  Sp_v = 0.

  call velocity3d_boundary("v")

  ! ====================== W-Velocity ====================== !
  ! Initialize field :: Set to 0
  w = 0.
  w_hat = 0.
  w_star = 0.

  ! Initialize coefficients
  Aw_w = 0.
  Ae_w = 0.
  As_w = 0.
  An_w = 0.
  Ab_w = 0.
  At_w = 0.
  Ap_w = 1.
  b_w = 0.

  ! Initialize source terms
  Su_w = 0.
  Sp_w = 0.

  call velocity3d_boundary("w")

  return

end subroutine velocity3d
