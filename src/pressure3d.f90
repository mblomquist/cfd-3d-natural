! pressure3d subroutines
!
! Written by Matt Blomquist
! Last Update: 2018-09-07 (YYYY-MM-DD)
!
! The subroutines below serve to calculate the coefficients and source terms
! for pressure in a 3D natural convection code. Additionally,
! the pressure solve and pressure correct subroutines are also
! included.
!
subroutine pressure3d_solve

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k, fault

  ! Initalize Coefficients
  Aw_p = 0.
  Ae_p = 0.
  As_p = 0.
  An_p = 0.
  Ab_p = 0.
  At_p = 0.
  Ap_p = 1.
  b_p = 0.

  ! Calculate Coefficients for Pressure
  do i = 1, m-1
    do j = 1, n-1
	  do k = 1, l-1

        Aw_p(i,j,k) = dy*dz*dy*dz/Ap_u(i,j,k)
        Ae_p(i,j,k) = dy*dz*dy*dz/Ap_u(i+1,j,k)
        As_p(i,j,k) = dx*dz*dz*dx/Ap_v(i,j,k)
        An_p(i,j,k) = dx*dz*dz*dx/Ap_v(i,j+1,k)
        Ab_p(i,j,k) = dy*dx*dx*dy/Ap_w(i,j,k)
        At_p(i,j,k) = dy*dx*dx*dy/Ap_w(i,j,k+1)

        ! Check Boundary :: West // Symmetry
        if (i .eq. 1) then
          Aw_p(i,j,k) = 0.
        end if

        ! Check Boundary :: East // Symmetry
        if (i .eq. m-1) then
          Ae_p(i,j,k) = 0.
        end if

        ! Check Boundary :: South // Symmetry
        if (j .eq. 1) then
          As_p(i,j,k) = 0.
        end if

        ! Check Boundary :: Noth // Symmetry
        if (j .eq. n-1) then
          An_p(i,j,k) = 0.
        end if

        ! Check Boundary :: Bottom // Wall
        if (k .eq. 1) then
          Ab_p(i,j,k) = 0.
        end if

        ! Check Boundary :: Top // Wall
        if (k .eq. l-1) then
          At_p(i,j,k) = 0.
        end if

        ! Update Ap_p
        Ap_p(i,j,k) = Aw_p(i,j,k)+Ae_p(i,j,k)+As_p(i,j,k)+An_p(i,j,k)+Ab_p(i,j,k)+At_p(i,j,k)

		    ! Calculate Mass Source Term from Pseudo-Velocities
		    b_p(i,j,k) = dz*dy*(u_hat(i,j,k)-u_hat(i+1,j,k)) + dx*dz*(v_hat(i,j,k)-v_hat(i,j+1,k)) + dx*dy*(w_hat(i,j,k)-w_hat(i,j,k+1))

	  end do
	end do
  end do

  ! Set reference pressure node
  Aw_p(1,1,1) = 0.
  Ae_p(1,1,1) = 0.
  As_p(1,1,1) = 0.
  An_p(1,1,1) = 0.
  Ab_p(1,1,1) = 0.
  At_p(1,1,1) = 0.
  Ap_p(1,1,1) = 1.
  b_p(1,1,1) = 0.

  ! Initialize Pressure Field to 0.
  P = 0.

  ! Solve Pressure Equation
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab_p, As_p, Aw_p, Ap_p, Ae_p, An_p, At_p, b_p, P, m-1, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab_p, As_p, Aw_p, Ap_p, Ae_p, An_p, At_p, b_p, P, m-1, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 2) then
    fault = 0
    do i = 3, maxit
      if (fault .eq. 0) then
        call solver3d_gmres(Ab_p, As_p, Aw_p, Ap_p, Ae_p, An_p, At_p, b_p, P, m-1, n-1, l-1, solver_tol, maxit, fault)
      end if
    end do
  elseif (solver .eq. 3) then
    call solver3d_bicg(Ab_p, As_p, Aw_p, Ap_p, Ae_p, An_p, At_p, b_p, P, m-1, n-1, l-1, solver_tol, maxit)
  else
    call solver3d_tdma(Ab_p, As_p, Aw_p, Ap_p, Ae_p, An_p, At_p, b_p, P, m-1, n-1, l-1, solver_tol, maxit)
  end if

  return
end subroutine pressure3d_solve

subroutine pressure3d_correct

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k, fault

  ! Calculate the new mass source term
  do i = 1, m-1
    do j = 1, n-1
	  do k = 1, l-1

	    b_p(i,j,k) = dz*dy*(u_star(i,j,k)-u_star(i+1,j,k)) + dx*dz*(v_star(i,j,k)-v_star(i,j+1,k)) + dx*dy*(w_star(i,j,k)-w_star(i,j,k+1))

	  end do
	end do
  end do

  ! Set reference pressure node
  Aw_p(1,1,1) = 0.
  Ae_p(1,1,1) = 0.
  As_p(1,1,1) = 0.
  An_p(1,1,1) = 0.
  Ab_p(1,1,1) = 0.
  At_p(1,1,1) = 0.
  Ap_p(1,1,1) = 1.
  b_p(1,1,1) = 0.

  ! Initialize Pressure Correction Field to 0.
  P_prime = 0.

  ! Solve Pressure Correction Equation
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
