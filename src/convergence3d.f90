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
  real(8) :: P_res, T_res, T_temp

  R_e(itr) = 0.0
  R_t(itr) = 0.0

  P_res = 0.
  T_res = 0.

  do i = 2,m-2
    do j = 2,n-2
	  do k = 2,l-2
	    
		if (P_res .le. abs(b_p(i,j,k))) then
		  P_res = abs(b_p(i,j,k))
		end if

		T_temp = abs(Ap_T(i,j,k)*T(i,j,k)- &
		             Aw_T(i,j,k)*T(i-1,j,k)- &
					 Ae_T(i,j,k)*T(i+1,j,k)- &
					 As_T(i,j,k)*T(i,j-1,k)- &
					 An_T(i,j,k)*T(i,j+1,k)- &
					 Ab_T(i,j,k)*T(i,j,k-1)- &
					 At_T(i,j,k)*T(i,j,k+1)- &
					 b_T(i,j,k))

		if (T_res .le. T_temp) then
		  T_res = T_temp
		end if

	  end do
	end do
  end do

  R_e(itr) = P_res
  R_t(itr) = T_res

  return

end subroutine convergence3d
