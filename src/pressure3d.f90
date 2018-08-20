! pressure3d Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms to 0, sets boundary values, and initializes the pressure grid.
!

subroutine pressure3d_init

  ! Pull in standard variable header
  include "var3d.dec"

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

end subroutine pressure3d_init

subroutine pressure3d_solve

  implicit none

  ! Include variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k, fault

  ! Initialize coefficients
  Aw_p = 0.
  Ae_p = 0.
  As_p = 0.
  An_p = 0.
  Ab_p = 0.
  At_p = 0.
  Ap_p = 1.
  b_p = 0.

  ! Update coefficients
  do i = 1,m-1
    do j = 1,n-1
      do k = 1,l-1

        ! Update coefficients
        Aw_p(i,j,k) = dy*dz/Ap_u(i,j,k)*alpha_v
        Ae_p(i,j,k) = dy*dz/Ap_u(i+1,j,k)*alpha_v
        As_p(i,j,k) = dz*dx/Ap_v(i,j,k)*alpha_v
        An_p(i,j,k) = dz*dx/Ap_v(i,j+1,k)*alpha_v
        Ab_p(i,j,k) = dx*dy/Ap_w(i,j,k)*alpha_v
        At_p(i,j,k) = dx*dy/Ap_w(i,j,k+1)*alpha_v

        ! Check West Wall
        if (i .eq. 1) then
          Aw_p(i,j,k) = 0.
        end if

        ! Check East Wall
        if (i .eq. m-1) then
          Ae_p(i,j,k) = 0.
        end if

        ! Check South Wall
        if (j .eq. 1) then
          As_p(i,j,k) = 0.
        end if

        ! Check North Wall
        if (j .eq. n-1) then
          An_p(i,j,k) = 0.
        end if

        ! Check Bottom Wall
        if (k .eq. 1) then
          Ab_p(i,j,k) = 0.
        end if

        ! Check Top Wall
        if (k .eq. l-1) then
          At_p(i,j,k) = 0.
        end if

        ! Update Ap_p
        Ap_p(i,j,k) = Aw_p(i,j,k)+Ae_p(i,j,k)+As_p(i,j,k)+An_p(i,j,k)+Ab_p(i,j,k)+At_p(i,j,k)

        ! Solve mass source term
        b_p(i,j,k) = (u_hat(i,j,k)-u_hat(i+1,j,k))*dy*dz+(v_hat(i,j,k)-v_hat(i,j+1,k))*dz*dx+(w_hat(i,j,k)-w_hat(i,j,k+1))*dx*dy

      end do
    end do
  end do

  ! Set reference pressure node (east-north corner)
  Aw_p(:,:,1) = 0.
  Ae_p(:,:,1) = 0.
  As_p(:,:,1) = 0.
  An_p(:,:,1) = 0.
  Ab_p(:,:,1) = 0.
  At_p(:,:,1) = 0.
  Ap_p(:,:,1) = 1.
  b_p(:,:,1) = 0.

  ! Initialize P
  P = 0.

  ! Solve pressure equation
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

  !print *, ".............."
  !print *, "Ab_p", Ab_p
  !print *, "As_p", As_p
  !print *, "Aw_p", Aw_p
  !print *, "Ap_p", Ap_p
  !print *, "Ae_p", Ae_p
  !print *, "An_p", An_p
  !print *, "At_p", At_p
  !print *, "b_p", b_p
  !print *, "P", P
  !print *, ".............."

  return

end subroutine pressure3d_solve

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

  ! Set reference pressure node (east-north corner)
  Aw_p(:,:,1) = 0.
  Ae_p(:,:,1) = 0.
  As_p(:,:,1) = 0.
  An_p(:,:,1) = 0.
  Ab_p(:,:,1) = 0.
  At_p(:,:,1) = 0.
  Ap_p(:,:,1) = 1.
  b_p(:,:,1) = 0.

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
