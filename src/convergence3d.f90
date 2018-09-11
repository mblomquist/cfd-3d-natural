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
  real(8) :: c_temp, e_temp

  ! Find maximum mass source term
  R_c(itr) = 0.0
  R_e(itr) = 0.0

  do i = 1, m-1
    do j = 1, n-1
	  do k = 1, l-1

	    ! Check Continuity Equation
	    c_temp = abs(b_p(i,j,k))

		if (R_c(itr) .le. c_temp) then
		  R_c(itr) = c_temp
		end if

		! Check Energy Equation
		if ((i .ne. 1) .and. (i .ne. m-1) .and. (j .ne. 1) .and. (j .ne. n-1) .and. (k .ne. 1) .and. (k .ne. l-1)) then
		  e_temp = abs(Ap_T(i,j,k)*T(i,j,k)- &
		               Aw_T(i,j,k)*T(i-1,j,k)- &
			  		   Ae_T(i,j,k)*T(i+1,j,k)- &
					   As_T(i,j,k)*T(i,j-1,k)- &
					   An_T(i,j,k)*T(i,j+1,k)- &
					   Ab_T(i,j,k)*T(i,j,k-1)- &
					   At_T(i,j,k)*T(i,j,k+1)- &
					   b_T(i,j,k))
		end if

		if (R_e(itr) .le. e_temp) then
		  R_e(itr) = e_temp
		end if

	  end do
	end do
  end do

  return
end subroutine convergence3d
