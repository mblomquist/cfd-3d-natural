! pressure3d Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms to 0, sets boundary values, and initializes the pressure grid.
!

subroutine pressure3d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Initialize Pressure Field
  P = 0
  P_star = 0
  P_prime = 0
  
  ! Initialize coefficients
  Aw_p = 0.
  Ae_p = 0.
  As_p = 0.
  An_p = 0.
  Ab_p = 0.
  At_p = 0.
  Ap_p = 1.
  b_p = 0.

  ! Set pressure source terms
  Su_p = 0
  Sp_p = 0

  return

end subroutine pressure3d
