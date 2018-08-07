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

  return

end subroutine temperature3d
