! initialize3d subroutine
!
! Written by Matt Blomquist
! Last Update: 2018-09-07 (YYYY-MM-DD)
!
! This subroutine reads the input file for a 3D natural convection problem
! and calculates the dimensionless quantities
!
subroutine initialize3d

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Read input file
  open(unit = 2, file = "input3d.txt")
  read(2,*)
  read(2,*)
  read(2,*) length, width, depth
  read(2,*)
  read(2,*) g, rho, mu, k_const, Cp, beta
  read(2,*)
  read(2,*) delta_T
  read(2,*)
  read(2,*) itrmax, maxit, solver_tol, simpler_tolc, simpler_tole, alpha_v, alpha_t, solver
  read(2,*)
  read(2,*) beta_1

  ! Calculate properties
  alpha = k_const / Cp / rho
  nu = mu / rho

  ! Calculate dimensionless numbers
  Ra = g*beta*delta_T*depth**3.0/alpha/nu
  Pr = nu/alpha

  ! Calculate characteristic velocity
  u0 = nu/depth

  ! Calculate dx, dy, dz
  dx = length/m/depth
  dy = width/n/depth
  dz = depth/l/depth

  ! Set inital error to 0
  R_c = 0.0
  R_e = 0.0

  return

 end subroutine initialize3d
