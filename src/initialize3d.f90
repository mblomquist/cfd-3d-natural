! initialize3d Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-16 (YYYY-MM-DD)
!

subroutine initialize3d

  ! Pull in standard variable header
  include "var2d.dec"

  ! Read Input File ..........
  open(unit = 2, file = "input3d.txt")
  read(2,*)
  read(2,*)
  read(2,*) length, width, depth
  read(2,*)
  read(2,*) g, rho, mu, k_const, Cp, beta
  read(2,*)
  read(2,*) u_bc_wv, u_bc_ev, u_bc_nv, u_bc_sv, u_bc_bv, u_bc_tv
  read(2,*)
  read(2,*) u_bc_wc, u_bc_ec, u_bc_nc, u_bc_sc, u_bc_bc, u_bc_tc
  read(2,*)
  read(2,*) v_bc_wv, v_bc_ev, v_bc_nv, v_bc_sv, v_bc_bv, v_bc_tv
  read(2,*)
  read(2,*) v_bc_wc, v_bc_ec, v_bc_nc, v_bc_sc, v_bc_bc, v_bc_tc
  read(2,*)
  read(2,*) w_bc_wv, w_bc_ev, w_bc_nv, w_bc_sv, w_bc_bv, w_bc_tv
  read(2,*)
  read(2,*) w_bc_wc, w_bc_ec, w_bc_nc, w_bc_sc, w_bc_bc, w_bc_tc
  read(2,*)
  read(2,*) T_bc_wv, T_bc_ev, T_bc_nv, T_bc_sv, T_bc_bv, T_bc_tv
  read(2,*)
  read(2,*) T_bc_wc, T_bc_ec, T_bc_nc, T_bc_sc, T_bc_bc, T_bc_tc
  read(2,*)
  read(2,*) itrmax, maxit, solver_tol, simpler_tol, alpha_v, alpha_t, solver
  close(2)

  ! Calculate parameters
  alpha = k_const / Cp / rho
  nu = mu / rho
  delta_T = T_h - T_c

  ! Calculate Characteristic Velocity
  u0 = (g*beta*delta_T*L)**(0.5)

  ! Calculate dimensionless numbers
  Ra = g*beta*delta_T*length**3.0/alpha/nu
  Pr = nu/alpha

  ! Calculate geometry properties.
  call geometry3d

  ! Calculate pressure source terms and initialize grid
  call pressure3d

  ! Calculate velocity source terms and initialize grid
  call velocity3d

  ! Calculate temperature source terms and initialize grid
  call temperature3d

  return

end subroutine initialize2d
