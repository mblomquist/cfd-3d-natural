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

  do i = 1,itrmax

    ! Calculate velocity coefficients
    call velocity3d_source("u")
    call velocity3d_source("v")
    call velocity3d_source("w")

    ! Step 2: Calculate Pseudo-Velocities
    !print *, "Step 1: Solve Pseudo-Velocities"
    call cpu_time(t_start)
    call pseudo3d_solve
    call cpu_time(t_end)
    t_1(i) = t_end - t_start

    ! Step 3: Solve Pressure Equation
    !print *, "Step 2: Solve Pressure Equation"
    call cpu_time(t_start)
    call pressure3d_solve
    call cpu_time(t_end)
    t_2(i) = t_end - t_start

    ! Set p_star := P
    P_star = P

    ! Step 4: Solve Momentum Equations
    !print *, "Step 4: Solve Momentum Equations"
    call cpu_time(t_start)
    call velocity3d_solve
    call cpu_time(t_end)
    t_4(i) = t_end - t_start

    ! Step 8: Check Convergence
    !print *, "Check Convergence"
    call cpu_time(t_start)
    call convergence3d(i)
    call cpu_time(t_end)
    t_3(i) = t_end - t_start

    if (i .eq. 1) then

      ! Print Current Information to Terminal
      print *, ""
	    print *, "Iteration:", i
      print *, "Continuity Error: ", R_e(i,1), R_e(i,2)
      print *, "X Momentum Error: ", R_u(i,1), R_u(i,2)
      print *, "Y Momentum Error: ", R_v(i,1), R_v(i,2)
      print *, "Z Momentum Error: ", R_w(i,1), R_w(i,2)
      print *, "Temperature Error:", R_t(i,1), R_t(i,2)
      print *, ""

    else

	    !  Print Current Information to Terminal
      print *, ""
	    print *, "Iteration:", i
      print *, "Continuity Error: ", R_e(i,1), R_e(i,2)
      print *, "X Momentum Error: ", R_u(i,1), R_u(i,2)
      print *, "Y Momentum Error: ", R_v(i,1), R_v(i,2)
      print *, "Z Momentum Error: ", R_w(i,1), R_w(i,2)
      print *, "Temperature Error:", R_t(i,1), R_t(i,2)
      print *, ""

      ! Check for Convergence
      if ((R_e(i,2) .le. simpler_tol) .and. (R_t(i,2).le. simpler_tol)) then

        call temperature3d_solve(1)
        print *, "Simpler completed in: ", i
        exit

      end if
    end if 

    ! Step 5: Solve Pressure Equation
    !print *, "Step 5: Solve Pressure Correction"
    call cpu_time(t_start)
    call pressure3d_correct
    call cpu_time(t_end)
    t_5(i) = t_end - t_start

    ! Step 6: Correct Velocities
    !print *, "Step 6: Correct Velocities"
    call cpu_time(t_start)
    call velocity3d_correct
    call cpu_time(t_end)
    t_6(i) = t_end - t_start

    ! Step 7: Solve Temperature Equation
    !print *, "Step 7: Solve Temperature Equation"
    call cpu_time(t_start)
    call temperature3d_solve(1)
    call cpu_time(t_end)
    t_7(i) = t_end - t_start

    ! Reset Initial  Initial Guesses
	  u_star = u
	  v_star = v
    w_star = w

  end do

  return

end subroutine simpler3d
