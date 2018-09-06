! convergence2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!
! This subroutine computes the error of the SIMPLER solution based on the
! residuals for the u- and v-momentum equations.
!

subroutine convergence3d(itr)

  ! Include variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k, itr
  real(8) :: e2_sum, u2_sum, v2_sum, w2_sum, T2_sum

  R_e(itr, :) = 0.
  R_u(itr, :) = 0.
  R_v(itr, :) = 0.
  R_w(itr, :) = 0.
  R_t(itr, :) = 0.

  P_res = 0.
  T_res = 0.
  U_res = 0.
  V_res = 0.
  W_res = 0.

  do i = 2,m-2
    do j = 2,n-2
      do k = 2,l-2

        if (P_res .le. abs(b_p(i,j,k))) then
          P_res = abs(b_p(i,j,k))
        end if

        U_temp(i,j,k) = abs(Ap_u(i,j,k)*u_star(i,j,k)- &
                 Aw_u(i,j,k)*u_star(i-1,j,k)- &
                 Ae_u(i,j,k)*u_star(i+1,j,k)- &
                 As_u(i,j,k)*u_star(i,j-1,k)- &
                 An_u(i,j,k)*u_star(i,j+1,k)- &
                 Ab_u(i,j,k)*u_star(i,j,k-1)- &
                 At_u(i,j,k)*u_star(i,j,k+1)- &
                 b_u(i,j,k))

        if (U_res .le. U_temp(i,j,k)) then
          U_res = U_temp(i,j,k)
        end if

        V_temp(i,j,k) = abs(Ap_v(i,j,k)*v_star(i,j,k)- &
                 Aw_v(i,j,k)*v_star(i-1,j,k)- &
                 Ae_v(i,j,k)*v_star(i+1,j,k)- &
                 As_v(i,j,k)*v_star(i,j-1,k)- &
                 An_v(i,j,k)*v_star(i,j+1,k)- &
                 Ab_v(i,j,k)*v_star(i,j,k-1)- &
                 At_v(i,j,k)*v_star(i,j,k+1)- &
                 b_v(i,j,k))

        if (V_res .le. V_temp(i,j,k)) then
          V_res = V_temp(i,j,k)
        end if

        W_temp(i,j,k) = abs(Ap_w(i,j,k)*w_star(i,j,k)- &
                 Aw_w(i,j,k)*w_star(i-1,j,k)- &
                 Ae_w(i,j,k)*w_star(i+1,j,k)- &
                 As_w(i,j,k)*w_star(i,j-1,k)- &
                 An_w(i,j,k)*w_star(i,j+1,k)- &
                 Ab_w(i,j,k)*w_star(i,j,k-1)- &
                 At_w(i,j,k)*w_star(i,j,k+1)- &
                 b_w(i,j,k))

        if (W_res .le. W_temp(i,j,k)) then
          W_res = W_temp(i,j,k)
        end if

        T_temp(i,j,k) = abs(Ap_T(i,j,k)*T(i,j,k)- &
		             Aw_T(i,j,k)*T(i-1,j,k)- &
                 Ae_T(i,j,k)*T(i+1,j,k)- &
                 As_T(i,j,k)*T(i,j-1,k)- &
                 An_T(i,j,k)*T(i,j+1,k)- &
                 Ab_T(i,j,k)*T(i,j,k-1)- &
                 At_T(i,j,k)*T(i,j,k+1)- &
                 b_T(i,j,k))

        if (T_res .le. T_temp(i,j,k)) then
          T_res = T_temp(i,j,k)
        end if

      end do
    end do
  end do

  if (itr .eq. 1) then

    R_e(itr, :) = 1.
    R_u(itr, :) = 1.
    R_v(itr, :) = 1.
    R_w(itr, :) = 1.
    R_t(itr, :) = 1.

  else

    R_e(itr, 1) = P_res
    R_u(itr, 1) = u_res
    R_v(itr, 1) = v_res
    R_w(itr, 1) = w_res
    R_t(itr, 1) = T_res

    e2_sum = 0.0
    u2_sum = 0.0
    v2_sum = 0.0
    w2_sum = 0.0
    T2_sum = 0.0

    do i = itr, 1, -1

      e2_sum = R_e(i,1)**2.0 + e2_sum
      u2_sum = R_u(i,1)**2.0 + u2_sum
      v2_sum = R_v(i,1)**2.0 + v2_sum
      w2_sum = R_w(i,1)**2.0 + w2_sum
      T2_sum = R_T(i,1)**2.0 + T2_sum

    end do

    R_e(itr, 2) = (e2_sum/itr)**(0.5)
    R_u(itr, 2) = (u2_sum/itr)**(0.5)
    R_v(itr, 2) = (v2_sum/itr)**(0.5)
    R_w(itr, 2) = (w2_sum/itr)**(0.5)
    R_t(itr, 2) = (T2_sum/itr)**(0.5)

    R_e(itr, 3) = abs(abs(R_e(itr-1, 2))-abs(R_e(itr, 2)))/abs(R_e(itr, 2))
    R_u(itr, 3) = abs(abs(R_u(itr-1, 2))-abs(R_u(itr, 2)))/abs(R_u(itr, 2))
    R_v(itr, 3) = abs(abs(R_v(itr-1, 2))-abs(R_v(itr, 2)))/abs(R_v(itr, 2))
    R_w(itr, 3) = abs(abs(R_w(itr-1, 2))-abs(R_w(itr, 2)))/abs(R_w(itr, 2))
    R_t(itr, 3) = abs(abs(R_t(itr-1, 2))-abs(R_t(itr, 2)))/abs(R_t(itr, 2))

  end if

  return

end subroutine convergence3d
