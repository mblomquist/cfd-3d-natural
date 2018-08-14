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
  call temperature3d_solve(0)

  !print *, "................"
  !print *, "T:", T
  !print *, "................"

  do i = 1,itrmax

    ! Calculate velocity coefficients
    call velocity3d_source("u")
    call velocity3d_source("v")
    call velocity3d_source("w")

    print *, "................"
    print *, "Ab_v:", Ab_v
    print *, "As_v:", As_v
    print *, "Aw_v:", Aw_v
    print *, "Ap_v:", Ap_v
    print *, "Ae_v:", Ae_v
    print *, "An_v:", An_v
    print *, "At_v:", At_v
    print *, "b_v:", b_v
    print *, "................"

    ! Step 2: Calculate Pseudo-Velocities
    print *, "Step 1: Solve Pseudo-Velocities"
    call cpu_time(t_start)
    call pseudo3d_solve
    call cpu_time(t_end)
    t_1(i) = t_end - t_start

    print *, "................"
    !print *, "u_hat:", u_hat
    print *, "v_hat:", v_hat
    !print *, "w_hat:", w_hat
    print *, "................"

    ! Step 3: Solve Pressure Equation
    print *, "Step 2: Solve Pressure Equation"
    call cpu_time(t_start)
    call pressure3d_solve
    call cpu_time(t_end)
    t_2(i) = t_end - t_start

    ! Step 8: Check Convergence
    print *, "Check Convergence"
    call cpu_time(t_start)
    call convergence3d(i)
    call cpu_time(t_end)
    t_3(i) = t_end - t_start

    if (i .eq. 1) then

      ! Print Current Information to Terminal
      print *, "Iteration:", i
      print *, "Relative Momentum Error: ", R_e(i)
      print *, "Relative Energy Error:", R_t(i)

    else

	    !  Print Current Information to Terminal
	    print *, "Iteration:", i
      print *, "Relative Momentum Error: ", R_e(i)
      print *, "Relative Energy Error:", R_t(i)

      ! Check for Convergence
      if ((R_e(i) .le. simpler_tol) .and. (R_t(i) .le. simpler_tol)) then

        call temperature3d_solve(1)
        print *, "Simpler completed in: ", i
        exit

      end if
    end if

	  ! Set p_star := P
	  P_star = P

    !print *, "................"
    !print *, "P_star:", P_star
    !print *, "................"

    ! Step 4: Solve Momentum Equations
    print *, "Step 4: Solve Momentum Equations"
    call cpu_time(t_start)
    call velocity3d_solve
    call cpu_time(t_end)
    t_4(i) = t_end - t_start

    print *, "................"
    print *, "b_v:", b_v
    print *, "................"

    print *, "................"
    !print *, "u_star:", u_star
    print *, "v_star:", v_star
    !print *, "w_star:", w_star
    print *, "................"

    ! Step 5: Solve Pressure Equation
    print *, "Step 5: Solve Pressure Correction"
    call cpu_time(t_start)
    call pressure3d_correct
    call cpu_time(t_end)
    t_5(i) = t_end - t_start

    !print *, "................"
    !print *, "P_prime:", P_Prime
    !print *, "................"

    ! Step 6: Correct Velocities
    print *, "Step 6: Correct Velocities"
    call cpu_time(t_start)
    call velocity3d_correct
    call cpu_time(t_end)
    t_6(i) = t_end - t_start

    print *, "................"
    !print *, "u:", u
    print *, "v:", v
    !print *, "w:", w
    print *, "................"

    ! Step 7: Solve Temperature Equation
    print *, "Step 7: Solve Temperature Equation"
    call cpu_time(t_start)
    call temperature3d_solve(1)
    call cpu_time(t_end)
    t_7(i) = t_end - t_start

    !print *, "................"
    !print *, "T:", T
    !print *, "................"

    ! Reset Initial  Initial Guesses
	  u_star = u
	  v_star = v
    w_star = w

  end do

  return

end subroutine simpler3d
