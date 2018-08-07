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
  write(2, *), 'ZONE I=32, J=32, K=16, DATAPACKING=POINT'

  do k = 1,l-1
    do j = 1,n-1
      do i = 1,m-1

        write(2, *), length*i*dx, length*j*dy, length*k*dz, P(i,j,k), (u(i,j,k)+u(i+1,j,k))/2, (v(i,j,k)+v(i,j+1,k))/2, &
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
  write (5, *), "R_e, R_t"
  do i = 1, itrmax
    write (5,*), R_e(i), R_t(i)
  end do

  close(5)

  return

end subroutine output3d_results
