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

  Sp_u(:,:,1) = 2.*dy*dz/dx*Pr*(Pr/Ra)**(0.5)
  Sp_u(:,:,l-1) = 2.*dy*dz/dx*Pr*(Pr/Ra)**(0.5)

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

  Sp_u(:,:,1) = 2.*dz*dx/dy*Pr*(Pr/Ra)**(0.5)
  Sp_u(:,:,l-1) = 2.*dz*dx/dy*Pr*(Pr/Ra)**(0.5)

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

  return

end subroutine velocity3d
