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

  call cpu_time(t_step(1,1))

  ! Start SIMPLER Algorithm
  do i = 1, itrmax

    if (i .ge. 2) then
      call cpu_time(t_step(i,1))
    end if

  	! Step 2: Calculate coefficients for velocity and
  	!         solve for u_hat, v_hat, w_hat
  	call velocity3d_source
  	call velocity3d_pseudo

    call cpu_time(t_step(i,2))

  	! Step 3: Calculate coefficients for pressure and
  	!         solve for P
  	call pressure3d_solve

    call cpu_time(t_step(i,3))

  	! Check for Convergence
  	call convergence3d(i)

  	if ((R_c(i,3) .le. simpler_tolc) .and. (R_e(i,3) .le. simpler_tole)) then

  	  ! Set u_hat, v_hat, w_hat to velocity solution
  	  u = u_hat
  	  v = v_hat
  	  w = w_hat

  	  ! Calculate the coefficients for temperature and
  	  ! solve for T
  	  call temperature3d_source
  	  call temperature3d_solve(1)

  	  print *, "SIMPLER Converged in: ", i
  	  print *, "Continuity Error: ", R_c(i,1), R_c(i,3)
      print *, "U-Momentum Error: ", R_u(i,1), R_u(i,3)
      print *, "V-Momentum Error: ", R_v(i,1), R_v(i,3)
      print *, "W-Momentum Error: ", R_w(i,1), R_w(i,3)
  	  print *, "Energy Error: ", R_e(i,1), R_e(i,3)
  	  print *, ""

      return

  	else

  	  print *, "Iteration:", i
      print *, "Continuity Error: ", R_c(i,1), R_c(i,3)
      print *, "U-Momentum Error: ", R_u(i,1), R_u(i,3)
      print *, "V-Momentum Error: ", R_v(i,1), R_v(i,3)
      print *, "W-Momentum Error: ", R_w(i,1), R_w(i,3)
  	  print *, "Energy Error: ", R_e(i,1), R_e(i,3)
  	  print *, ""

  	end if

    call cpu_time(t_step(i,4))

  	! Step 4: Solve the momentum equations
  	call velocity3d_solve

    call cpu_time(t_step(i,5))

  	! Step 5: Calculate mass source and solve the
  	!         pressure correction equation, P_prime
  	call pressure3d_correct

    call cpu_time(t_step(i,6))

  	! Step 6: Correct the velocity field
  	call velocity3d_correct

    call cpu_time(t_step(i,7))

  	! Step 7: Calculate the coefficients for temperature and
  	!         solve for T
  	call temperature3d_source
  	call temperature3d_solve(1)

    call cpu_time(t_step(i,8))

  end do

  return

end subroutine simpler3d
