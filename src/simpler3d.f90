! simpler3d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine runs the SIMPLER algorithm for a 3D CFD problem.
!
subroutine simpler3d

  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define Internal Variables
  integer :: i, j, k

  print *, 'Start SIMPLER Algorithm.'

  ! Solve Temperature for Natural Convection First
  !print *, "Step 0: Solve Temperature Equation"
  call temperature3d_solve
  !print *, "................"
  !print *, "T:", T
  !print *, "................"

  do i = 1,itrmax

    ! Calculate velocity coefficients
    call velocity3d_boundary
    call velocity3d_source("u")
    call velocity3d_source("v")
    call velocity3d_source("w")

    ! Step 2: Calculate Pseudo-Velocities
    !print *, "Step 1: Solve Pseudo-Velocities"
    !call pseudo3d_solve
    u_hat = u
    v_hat = v
    w_hat = w

    !print *, "................"
    !print *, "u_hat:", u_hat
    !print *, "v_hat:", v_hat
    !print *, "w_hat:", w_hat
    !print *, "................"

    ! Step 3: Solve Pressure Equation
    !print *, "Step 2: Solve Pressure Equation"
    call pressure3d_solve

	  ! Set p_star := P
	  P_star = P

    !print *, "................"
    !print *, "P_star:", P_star
    !print *, "................"

    ! Step 4: Solve Momentum Equations
    !print *, "Step 4: Solve Momentum Equations"
    call velocity3d_solve
    !print *, "................"
    !print *, "u_star:", u_star
    !print *, "v_star:", v_star
    !print *, "w_star:", w_star
    !print *, "................"

    ! Step 5: Solve Pressure Equation
    !print *, "Step 5: Solve Pressure Correction"
    call pressure3d_correct

    ! Step 6: Correct Velocities
    !print *, "Step 6: Correct Velocities"
    call velocity3d_correct

    !print *, "................"
    !print *, "u:", u
    !print *, "v:", v
    !print *, "w:", w
    !print *, "................"

    ! Step 7: Solve Temperature Equation
    !print *, "Step 7: Solve Temperature Equation"
    call temperature3d_solve

    !print *, "................"
    !print *, "T:", T
    !print *, "................"

    ! Step 8: Check Convergence
    print *, "Check Convergence"
    call convergence3d(i)

    if (i .eq. 1) then

	  ! Print Current Information to Terminal
	  print *, "Iteration:", i
      print *, "Relative Momentum Error: ", R_e(i)
      print *, "Relative Energy Error:", R_t(i)

    else

	  ! Print Current Information to Terminal
	  print *, "Iteration:", i
      print *, "Relative Momentum Error: ", R_e(i)
      print *, "Relative Energy Error:", R_t(i)

	  ! Check for Convergence
	  if ((R_e(i) .le. simpler_tol) .and. (R_t(i) .le. simpler_tol)) then

        print *, "Simpler completed in: ", i
        exit

      end if

	end if

    ! Reset Initial  Initial Guesses
    P_star = P
	  u_star = u
	  v_star = v
    w_star = w

  end do

  return

end subroutine simpler3d
