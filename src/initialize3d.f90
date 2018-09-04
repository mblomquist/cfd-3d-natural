! initialize3d Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-16 (YYYY-MM-DD)
!

subroutine initialize3d

  ! Pull in standard variable header
  include "var3d.dec"

  ! Read Input File ..........
  open(unit = 2, file = "input3d.txt")
  read(2,*)
  read(2,*)
  read(2,*) length, width, depth
  read(2,*)
  read(2,*) g, rho, mu, k_const, Cp, beta
  read(2,*)
  read(2,*) delta_T
  read(2,*)
  read(2,*) itrmax, maxit, solver_tol, simpler_tol, alpha_v, alpha_t, solver
  read(2,*)
  read(2,*) beta_u, beta_v
  read(2,*)
  close(2)

  ! Calculate parameters
  alpha = k_const / Cp / rho
  nu = mu / rho

  ! Calculate dimensionless numbers
  Ra = g*beta*delta_T*depth**3.0/alpha/nu
  Pr = nu/alpha

  ! Calculate Characteristic Velocity
  u0 = nu/depth
  Re = u0*length/nu

  ! Calculate geometry properties.
  call geometry3d

  ! Calculate pressure source terms and initialize grid
  call pressure3d_init

  ! Calculate velocity source terms and initialize grid
  call velocity3d_init

  ! Calculate temperature source terms and initialize grid
  call temperature3d_init

  return

end subroutine initialize3d
