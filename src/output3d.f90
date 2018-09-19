! output3d subroutine
!
! Written by Matt Blomquist
! Last Update: 2018-09-07 (YYYY-MM-DD)
!
! The subroutines below outputs the pressure, velocity, and temperature
! fields. Additionally, the error at each iteration of the SIMPLER algorithm
! is plotted.
!
subroutine output3d

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! Write field data to file
  open(unit=2, file="output/output3d_results.dat")

  write(2, *), 'TITLE = "3D CFD - Natural Convection - Field Scalar Data"'
  write(2, *), 'VARIABLES = "X", "Y", "Z", "Pressure", "U-Velocity", "V-Velocity", "W-Velocity", "Temperature"'
  write(2, '("ZONE I=", I3, 1x, "J=", I3, 1x, "K=", I3, 1x, "DATAPACKING=POINT", /)'), m-1, n-1, l-1

  do k = 1,l-1
    do j = 1,n-1
      do i = 1,m-1

        write(2, *), depth*i*dx, depth*j*dy, depth*k*dz, P(i,j,k), (u(i,j,k)+u(i+1,j,k))/2, (v(i,j,k)+v(i,j+1,k))/2, &
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
  write (5, '("SIMPLER Algorithm Duration:", E15.4, /)', advance="no"), tf-t0
  write (5, *), "Solver #:", solver
  write (5, *), "R_c, R_e, R_u, R_v, R_w, dR_c, dR_e, dR_u, dR_v, dR_w"
  do i = 1, itrmax
    write (5,'(E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5, /)', advance="no"), &
           R_c(i,1), R_e(i,1), R_u(i,1), R_v(i,1), R_w(i,1), R_c(i,3), R_e(i,3), R_u(i,3), R_v(i,3), R_w(i,3)
  end do
  close(5)

  open(unit=8, file="output/time_data.dat")
  write (8, '("Grid size: ", 5I, 1X, 5I, 5I, /)'), m, n, l
  write (8, '("Rayleigh Number: ", E15.4, /)', advance="no"), Ra
  write (8, '("Prandtl Number: ", E15.4, /)', advance="no"), Pr
  write (8, '("delta_T: ", E15.4, /)', advance="no"), delta_T
  write (8, '("SIMPLER Algorithm Duration:", E15.4, /)', advance="no"), tf-t0
  write (8, '("SIMPLER Algorithm Start:", E15.4, /)', advance="no"), t0
  write (8, *), "Solver #:", solver
  write (8, *), "v_sol1, p_sol1, conv1, v_sol2, p_sol2, v_cor1, t_sol1"
  do i = 1, itrmax
    write (8,'(E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5, ",", E12.5, /)', advance="no"), &
           t_step(i,2)-t_step(i,1), t_step(i,3)-t_step(i,2), t_step(i,4)-t_step(i,3), t_step(i,5)-t_step(i,4), &
           t_step(i,6)-t_step(i,5), t_step(i,7)-t_step(i,6), t_step(i,8)-t_step(i,7)
  end do
  close(8)

  return
end subroutine output3d
