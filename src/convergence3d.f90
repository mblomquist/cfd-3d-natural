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
  real(8) :: c_temp, e_temp, c_rms, e_rms

  ! Find maximum mass source term
  R_c(itr,1) = 0.0
  R_e(itr,1) = 0.0

  c_rms = 0.0
  e_rms = 0.0

  do i = 2, m-2
    do j = 2, n-2
      do k = 2, l-2
        ! Check Continuity Equation
        c_temp = abs(b_p(i,j,k))

        if (R_c(itr,1) .le. c_temp) then
          R_c(itr,1) = c_temp
        end if

        ! Check Energy Equation
        e_temp = abs(Ap_T(i,j,k)*T(i,j,k)- &
                     Aw_T(i,j,k)*T(i-1,j,k)- &
                     Ae_T(i,j,k)*T(i+1,j,k)- &
                     As_T(i,j,k)*T(i,j-1,k)- &
                     An_T(i,j,k)*T(i,j+1,k)- &
                     Ab_T(i,j,k)*T(i,j,k-1)- &
                     At_T(i,j,k)*T(i,j,k+1)- &
                     b_T(i,j,k))

        if (R_e(itr,1) .le. e_temp) then
          R_e(itr,1) = e_temp
        end if

      end do
    end do
  end do

  if (itr .eq. 1) then

    R_c(itr,2) = 1.0
    R_e(itr,2) = 1.0

  else

    R_c(itr,2) = abs(R_c(itr,1)-R_c(itr-1,1))/R_c(itr,1)
    R_e(itr,2) = abs(R_e(itr,1)-R_e(itr-1,1))/R_e(itr,1)

  end if 

  return
end subroutine convergence3d
