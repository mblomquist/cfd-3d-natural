! velocity_solve2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!
! This subroutine updates the source terms for the solution of the momentum
! equations in the SIMPLER algorithm.
!
subroutine velocity3d_solve

  implicit none

  ! Include variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! ====================== U-Velocity ====================== !
  ! Update source terms :: u
  do i = 1,m
    do j = 1,n-1
	  do k = 1,l-1

        Ap_u(i,j,k) = Ap_u(i,j,k)/alpha_v

		if (i .eq. 1) then
		  b_u(i,j,k) = b_u(i,j,k)+dy*dz*(P_star(i-1,j,k)-P_star(i,j,k))+(1.0-alpha_v)*Ap_u(i,j,k)*u_hat(i,j,k)
		else
		  b_u(i,j,k) = b_u(i,j,k)+dy*dz*(P_star(i-1,j,k)-P_star(i,j,k))+(1.0-alpha_v)*Ap_u(i,j,k)*u_hat(i,j,k)
		end if

	  end do
    end do
  end do

  ! Solve u-velocity equation
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab_u, As_u, Aw_u, Ap_u, Ae_u, An_u, At_u, b_u, u_star, m, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab_u, As_u, Aw_u, Ap_u, Ae_u, An_u, At_u, b_u, u_star, m, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 2) then
    call solver3d_gmres(Ab_u, As_u, Aw_u, Ap_u, Ae_u, An_u, At_u, b_u, u_star, m, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 3) then
    call solver3d_paradiso(Ab_u, As_u, Aw_u, Ap_u, Ae_u, An_u, At_u, b_u, u_star, m, n-1, l-1, solver_tol, maxit)
  else
    call solver3d_tdma(Ab_u, As_u, Aw_u, Ap_u, Ae_u, An_u, At_u, b_u, u_star, m, n-1, l-1, solver_tol, maxit)
  end if

  ! ====================== V-Velocity ====================== !
  ! Update source terms :: v
  do i = 1, m-1
    do j = 1, n
	  do k = 1,l-1

        Ap_v(i,j,k) = Ap_v(i,j,k)/alpha_v

        if (j .eq. 1) then
		  b_v(i,j,k) = b_v(i,j,k)+(1.0-alpha_v)*Ap_v(i,j,k)*v_hat(i,j,k)
		else
		  b_v(i,j,k) = b_v(i,j,k)+dz*dx*(P_star(i,j-1,k)-P_star(i,j,k))+(1.0-alpha_v)*Ap_v(i,j,k)*v_hat(i,j,k)
		end if

	  end do
    end do
  end do

  ! Solve v-velocity equation
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab_v, As_v, Aw_v, Ap_v, Ae_v, An_v, At_v, b_v, v_star, m-1, n, l-1, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab_v, As_v, Aw_v, Ap_v, Ae_v, An_v, At_v, b_v, v_star, m-1, n, l-1, solver_tol, maxit)
  elseif (solver .eq. 2) then
    call solver3d_gmres(Ab_v, As_v, Aw_v, Ap_v, Ae_v, An_v, At_v, b_v, v_star, m-1, n, l-1, solver_tol, maxit)
  elseif (solver .eq. 3) then
    call solver3d_paradiso(Ab_v, As_v, Aw_v, Ap_v, Ae_v, An_v, At_v, b_v, v_star, m-1, n, l-1, solver_tol, maxit)
  else
    call solver3d_tdma(Ab_v, As_v, Aw_v, Ap_v, Ae_v, An_v, At_v, b_v, v_star, m-1, n, l-1, solver_tol, maxit)
  end if

  ! ====================== W-Velocity ====================== !
  ! Update source terms :: w
  do i = 1, m-1
    do j = 1, n-1
	  do k = 1,l

        Ap_w(i,j,k) = Ap_w(i,j,k)/alpha_v

		if (k .eq. 1) then
		  b_w(i,j,k) = b_w(i,j,k)+(1.0-alpha_v)*Ap_w(i,j,k)*w_hat(i,j,k)
		else
          b_w(i,j,k) = b_w(i,j,k)+dx*dy*(P_star(i,j,k-1)-P_star(i,j,k))+(1.0-alpha_v)*Ap_w(i,j,k)*w_hat(i,j,k)
		end if

	  end do
    end do
  end do

  ! Solve w-velocity equation
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab_w, As_w, Aw_w, Ap_w, Ae_w, An_w, At_w, b_w, w_star, m-1, n-1, l, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab_w, As_w, Aw_w, Ap_w, Ae_w, An_w, At_w, b_w, w_star, m-1, n-1, l, solver_tol, maxit)
  elseif (solver .eq. 2) then
    call solver3d_gmres(Ab_w, As_w, Aw_w, Ap_w, Ae_w, An_w, At_w, b_w, w_star, m-1, n-1, l, solver_tol, maxit)
  elseif (solver .eq. 3) then
    call solver3d_paradiso(Ab_w, As_w, Aw_w, Ap_w, Ae_w, An_w, At_w, b_w, w_star, m-1, n-1, l, solver_tol, maxit)
  else
    call solver3d_tdma(Ab_w, As_w, Aw_w, Ap_w, Ae_w, An_w, At_w, b_w, w_star, m-1, n-1, l, solver_tol, maxit)
  end if

  return

end subroutine velocity3d_solve
