! convergence3d subroutine
!
! Written by Matt Blomquist
! Last Update: 2018-09-07 (YYYY-MM-DD)
!
! This subroutine calculates the maximum mass-source value from
! the pressure equation to determine the continuity error of the
! 3D natural convection problem. Additionally, the error from the
! energy equation is calculated.
!
subroutine convergence3d(itr)

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k, itr
  real(8) :: c_temp, e_temp, u_temp, v_temp, w_temp
  real(8) :: c_rms, e_rms, u_rms, v_rms, w_rms

  ! Find maximum mass source term
  R_c(itr,1) = 0.0
  R_e(itr,1) = 0.0
  R_u(itr,1) = 0.0
  R_v(itr,1) = 0.0
  R_w(itr,1) = 0.0

  c_temp = 0.0
  e_temp = 0.0
  u_temp = 0.0
  v_temp = 0.0
  w_temp = 0.0

  c_rms = 0.0
  e_rms = 0.0
  u_rms = 0.0
  v_rms = 0.0
  w_rms = 0.0

  do i = 2, m-2
    do j = 2, n-2
      do k = 2, l-2

        ! Check Continuity Equation
        c_temp = abs(b_p(i,j,k))

        !if (R_c(itr,1) .le. c_temp) then
        !  R_c(itr,1) = c_temp
        !end if

        ! Check U-Momentum Equation
        u_temp = abs(Ap_u(i,j,k)*u_hat(i,j,k)- &
                     Aw_u(i,j,k)*u_hat(i-1,j,k)- &
                     Ae_u(i,j,k)*u_hat(i+1,j,k)- &
                     As_u(i,j,k)*u_hat(i,j-1,k)- &
                     An_u(i,j,k)*u_hat(i,j+1,k)- &
                     Ab_u(i,j,k)*u_hat(i,j,k-1)- &
                     At_u(i,j,k)*u_hat(i,j,k+1)- &
                     b_u(i,j,k))

       ! Check V-Momentum Equation
       v_temp = abs(Ap_v(i,j,k)*v_hat(i,j,k)- &
                    Aw_v(i,j,k)*v_hat(i-1,j,k)- &
                    Ae_v(i,j,k)*v_hat(i+1,j,k)- &
                    As_v(i,j,k)*v_hat(i,j-1,k)- &
                    An_v(i,j,k)*v_hat(i,j+1,k)- &
                    Ab_v(i,j,k)*v_hat(i,j,k-1)- &
                    At_v(i,j,k)*v_hat(i,j,k+1)- &
                    b_v(i,j,k))

        ! Check W-Momentum Equation
        w_temp = abs(Ap_w(i,j,k)*w_hat(i,j,k)- &
                     Aw_w(i,j,k)*w_hat(i-1,j,k)- &
                     Ae_w(i,j,k)*w_hat(i+1,j,k)- &
                     As_w(i,j,k)*w_hat(i,j-1,k)- &
                     An_w(i,j,k)*w_hat(i,j+1,k)- &
                     Ab_w(i,j,k)*w_hat(i,j,k-1)- &
                     At_w(i,j,k)*w_hat(i,j,k+1)- &
                     b_w(i,j,k))

        ! Check Energy Equation
        e_temp = abs(Ap_T(i,j,k)*T(i,j,k)- &
                     Aw_T(i,j,k)*T(i-1,j,k)- &
                     Ae_T(i,j,k)*T(i+1,j,k)- &
                     As_T(i,j,k)*T(i,j-1,k)- &
                     An_T(i,j,k)*T(i,j+1,k)- &
                     Ab_T(i,j,k)*T(i,j,k-1)- &
                     At_T(i,j,k)*T(i,j,k+1)- &
                     b_T(i,j,k))

        R_c(itr, 1) = R_c(itr, 1) + c_temp**(2.0)
        R_u(itr, 1) = R_u(itr, 1) + u_temp**(2.0)
        R_v(itr, 1) = R_v(itr, 1) + v_temp**(2.0)
        R_w(itr, 1) = R_w(itr, 1) + w_temp**(2.0)
        R_e(itr, 1) = R_e(itr, 1) + e_temp**(2.0)

      end do
    end do
  end do

  R_c(itr, 1) = (R_c(itr, 1) / ((m-3)*(n-3)*(l-3)))**(0.5)
  R_e(itr, 1) = (R_e(itr, 1) / ((m-3)*(n-3)*(l-3)))**(0.5)
  R_u(itr, 1) = (R_u(itr, 1) / ((m-3)*(n-3)*(l-3)))**(0.5)
  R_v(itr, 1) = (R_v(itr, 1) / ((m-3)*(n-3)*(l-3)))**(0.5)
  R_w(itr, 1) = (R_w(itr, 1) / ((m-3)*(n-3)*(l-3)))**(0.5)

  if (itr .ge. 2) then
    R_c(itr,2) = abs(R_c(itr,1)-R_c(itr-1,1))
    R_e(itr,2) = abs(R_e(itr,1)-R_e(itr-1,1))
    R_u(itr,2) = abs(R_u(itr,1)-R_u(itr-1,1))
    R_v(itr,2) = abs(R_v(itr,1)-R_v(itr-1,1))
    R_w(itr,2) = abs(R_w(itr,1)-R_w(itr-1,1))
  else
    R_c(itr,2) = 1.0
    R_e(itr,2) = 1.0
    R_u(itr,2) = 1.0
    R_v(itr,2) = 1.0
    R_w(itr,2) = 1.0
  end if

  if (itr .ge. 10) then
    do i = itr, itr-10, -1
      c_rms = c_rms + R_c(i,2)**(2.0)
      e_rms = e_rms + R_e(i,2)**(2.0)
      u_rms = u_rms + R_u(i,2)**(2.0)
      v_rms = v_rms + R_v(i,2)**(2.0)
      w_rms = w_rms + R_w(i,2)**(2.0)
    end do

    R_c(itr,3) = (c_rms/10.0)**(0.5)
    R_e(itr,3) = (e_rms/10.0)**(0.5)
    R_u(itr,3) = (u_rms/10.0)**(0.5)
    R_v(itr,3) = (v_rms/10.0)**(0.5)
    R_w(itr,3) = (w_rms/10.0)**(0.5)

  else
    do i = itr, 1, -1
      c_rms = c_rms + R_c(i,2)**(2.0)
      e_rms = e_rms + R_e(i,2)**(2.0)
      u_rms = u_rms + R_u(i,2)**(2.0)
      v_rms = v_rms + R_v(i,2)**(2.0)
      w_rms = w_rms + R_w(i,2)**(2.0)
    end do

    R_c(itr,3) = (c_rms/(itr))**(0.5)
    R_e(itr,3) = (e_rms/(itr))**(0.5)
    R_u(itr,3) = (u_rms/(itr))**(0.5)
    R_v(itr,3) = (v_rms/(itr))**(0.5)
    R_w(itr,3) = (w_rms/(itr))**(0.5)

  end if


  return
end subroutine convergence3d
