! temperature3d_solve Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine solves the temperature equation of the SIMPLER algorithm
!
subroutine temperature3d_solve(start)

  ! Include variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k, fault, start

  ! Update source terms
  call temperature3d_source

  if (start .eq. 0) then
    alpha_temp = 1.0
  else
    alpha_temp = alpha_t
  end if

  ! Update source terms
  do i = 2, m-2
    do j = 2, n-2
      do k = 2, l-2
        Ap_T(i,j,k) = Ap_T(i,j,k)/alpha_temp
        b_T(i,j,k) = Su_T(i,j,k)+(1.0-alpha_temp)*Ap_T(i,j,k)*T(i,j,k)
      end do
    end do
  end do

  ! Solve velocity Equations
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab_T, As_T, Aw_T, Ap_T, Ae_T, An_T, At_T, b_T, T, m-1, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab_T, As_T, Aw_T, Ap_T, Ae_T, An_T, At_T, b_T, T, m-1, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 2) then
    fault = 0
    do i = 3, maxit
      if (fault .eq. 0) then
        call solver3d_gmres(Ab_T, As_T, Aw_T, Ap_T, Ae_T, An_T, At_T, b_T, T, m-1, n-1, l-1, solver_tol, maxit, fault)
      end if
    end do
  elseif (solver .eq. 3) then
    call solver3d_bicg(Ab_T, As_T, Aw_T, Ap_T, Ae_T, An_T, At_T, b_T, T, m-1, n-1, l-1, solver_tol, maxit)
  else
    call solver3d_tdma(Ab_T, As_T, Aw_T, Ap_T, Ae_T, An_T, At_T, b_T, T, m-1, n-1, l-1, solver_tol, maxit)
  end if

  return

end subroutine temperature3d_solve
