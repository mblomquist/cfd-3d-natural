! solver2d_bicgstab
!
! Written by Matt Blomquist
! Last Update: 2018-02-02 (YYYY-MM-DD)
!
! This program solves a three-dimensional finite volume discretization problem
! using the bi-conjugate gradients stabilized algorithm.
!
! Definition of input arguments
! Inputs:
!   As, Aw, Ap, Ae, An :: These arrays represent the coefficients for adjacent nodes
!   b :: This array represents the right-hand side of the equation Ax=b
!   phi :: This value represents the appropriate solution array (pressure, velocity, temperature)
!   m, n :: These values represent the number of nodes for i and j for the phi value
!   tol :: represents the solution tolerance
!   maxit :: represents the maximum number of iterations of the BiCGStab Algorithm
!
! Outputs:
!   phi :: on exit, this value contains the updated solution
!   maxit :: on exit, this value contains the number of iterations of the BiCGStab algorithm
!   tol :: on exit, this value represents the normalized residual

subroutine solver3d_paradiso(Ab_u, As_u, Aw_u, Ap_u, Ae_u, An_u, At_u, b_u, u_star, m, n, l, solver_tol, maxit)

  return

end subroutine solver3d_paradiso
