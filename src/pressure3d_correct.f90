! pressure3d_correct Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine solves the pressure equation for the SIMPLER algorithm
! using u_star and v_star.

subroutine pressure3d_correct

  implicit none

  ! Include variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k, fault

  ! Update coefficients
  do i = 1,m-1
    do j = 1,n-1
	    do k = 1,l-1

	      ! Solve mass source term
		    b_p(i,j,k) = (u_star(i,j,k)-u_star(i+1,j,k))*dy*dz+(v_star(i,j,k)-v_star(i,j+1,k))*dz*dx+(w_star(i,j,k)-w_star(i,j,k+1))*dx*dy

	    end do
    end do
  end do

  do k = 1,l-1
    do j = 1,n-1
      do i = 1,m-1
        print *, i, j, Ap_p(i,j,k), b_p(i,j,k)
      end do
    end do
  end do

  ! Set reference pressure node (east-north corner)
  Aw_p(m-1,n-1,l-1) = 0.
  Ae_p(m-1,n-1,l-1) = 0.
  As_p(m-1,n-1,l-1) = 0.
  An_p(m-1,n-1,l-1) = 0.
  Ab_p(m-1,n-1,l-1) = 0.
  At_p(m-1,n-1,l-1) = 0.
  Ap_p(m-1,n-1,l-1) = 1.
  b_p(m-1,n-1,l-1) = 0.

  ! Initialize P_prime
  P_prime = 0.

  ! Solve pressure equation
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab_p, As_p, Aw_p, Ap_p, Ae_p, An_p, At_p, b_p, P_prime, m-1, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab_p, As_p, Aw_p, Ap_p, Ae_p, An_p, At_p, b_p, P_prime, m-1, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 2) then
    fault = 0
    do i = 3, maxit
      if (fault .eq. 0) then
        call solver3d_gmres(Ab_p, As_p, Aw_p, Ap_p, Ae_p, An_p, At_p, b_p, P_prime, m-1, n-1, l-1, solver_tol, maxit, fault)
      end if
    end do
  elseif (solver .eq. 3) then
    call solver3d_bicg(Ab_p, As_p, Aw_p, Ap_p, Ae_p, An_p, At_p, b_p, P_prime, m-1, n-1, l-1, solver_tol, maxit)
  else
    call solver3d_tdma(Ab_p, As_p, Aw_p, Ap_p, Ae_p, An_p, At_p, b_p, P_prime, m-1, n-1, l-1, solver_tol, maxit)
  end if

  return

end subroutine pressure3d_correct
