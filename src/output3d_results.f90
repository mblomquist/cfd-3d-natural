! output_results3d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine generates data files for the pressure, temperature, and velocity fields.

subroutine output3d_results

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! Write field data to file
  open(unit=2, file="output/output3d_results.dat")

  write(2, *), 'TITLE = "3D CFD - Natural Convection - Field Scalar Data"'
  write(2, *), 'VARIABLES = "X", "Y", "Z", "Pressure", "U-Velocity", "V-Velocity", "W-Velocity", "Temperature"'
  write(2, *), 'ZONE I=16, J=16, K=8, DATAPACKING=POINT'

  do k = 1,l-1
    do j = 1,n-1
      do i = 1,m-1

        write(2, *), i*dx, j*dy, k*dz, P(i,j,k), (u(i,j,k)+u(i+1,j,k))/2, (v(i,j,k)+v(i,j+1,k))/2, &
                    (w(i,j,k)+w(i,j,k+1))/2, T(i,j,k)

      end do
    end do
  end do

  close(2)

  open(unit=5, file="output/terminal_data.dat")
  write (5, '("Grid size: ", 5I, 1X, 5I, 5I, /)'), m, n, l
  write (5, '("Rayleigh Number: ", E15.4, /)', advance="no"), Ra
  write (5, '("Prandtl Number: ", E15.4, /)', advance="no"), Pr
  write (5, '("delta_T: ", E15.4, /)', advance="no"), delta_T
  write (5, '("SIMPLER Algorithm Duration:", E15.4, /)', advance="no"), end_time-start_time
  write (5, *), "Solver #:", solver
  write (5, *), "R_e, R_t, t_1, t_2, t_3, t_4, t_5, t_6, t_7"
  do i = 1, itrmax
    write (5,'(E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5)', advance="no"), &
               R_e(i), R_t(i), t_1(i), t_2(i), t_3(i), t_4(i), t_5(i), t_6(i), t_7(i)
  end do

  close(5)

  return

end subroutine output3d_results
