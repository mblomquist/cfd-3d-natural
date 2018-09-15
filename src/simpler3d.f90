! simpler3d subroutine
!
! Written by Matt Blomquist
! Last Update: 2018-09-07 (YYYY-MM-DD)
!
! This subroutine runs the SIMPLER algorithm for a 3D natural convection
! problem.
!
subroutine simpler3d

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i

  ! Step 0: Solve Inital Temperature Field (Assume Velocity = 0)
  call temperature3d_source
  call temperature3d_solve(0)
  !print *, "T:", T

  ! Step 1: Start with a guessed velocity field :: Set velocity to 0.
  u = 0.
  v = 0.
  w = 0.

  ! Start SIMPLER Algorithm
  do i = 1, itrmax

  	! Step 2: Calculate coefficients for velocity and
  	!         solve for u_hat, v_hat, w_hat
  	call velocity3d_source
  	call velocity3d_pseudo

  	! Step 3: Calculate coefficients for pressure and
  	!         solve for P
  	call pressure3d_solve

  	! Check for Convergence
  	call convergence3d(i)

  	if ((R_c(i,2) .le. simpler_tolc) .and. (R_e(i,1) .le. simpler_tole)) then

  	  ! Set u_hat, v_hat, w_hat to velocity solution
  	  u = u_hat
  	  v = v_hat
  	  w = w_hat

  	  ! Calculate the coefficients for temperature and
  	  ! solve for T
  	  call temperature3d_source
  	  call temperature3d_solve(1)

  	  print *, "SIMPLER Converged in: ", i
  	  print *, "Continuity Error: ", R_c(i,1), R_c(i,2)
  	  print *, "Energy Error: ", R_e(i,1), R_e(i,2)
  	  print *, ""

      return

  	else

  	  print *, "Iteration:", i
      print *, "Continuity Error: ", R_c(i,1), R_c(i,2)
  	  print *, "Energy Error: ", R_e(i,1), R_e(i,2)

  	end if

  	! Step 4: Solve the momentum equations
  	call velocity3d_solve

  	! Step 5: Calculate mass source and solve the
  	!         pressure correction equation, P_prime
  	call pressure3d_correct

  	! Step 6: Correct the velocity field
  	call velocity3d_correct

  	! Step 7: Calculate the coefficients for temperature and
  	!         solve for T
  	call temperature3d_source
  	call temperature3d_solve(1)

  end do

  return

end subroutine simpler3d
