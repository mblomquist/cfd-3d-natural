! velocity3d_correct Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine corrects the velocity values for a 2D
! CFD problem
!
subroutine velocity3d_correct

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! ====================== U-Velocity ====================== !
  ! Correct velocity values
  do i = 2,m-1
    do j = 1,n-1
      do k = 1,l-1

        u(i,j,k) = u_star(i,j,k)+dy/Ap_u(i,j,k)*(P_prime(i-1,j,k)-P_prime(i,j,k))*alpha_v

      end do
    end do
  end do

  ! ====================== V-Velocity ====================== !
  ! Correct velocity values
  do i = 1,m-1
    do j = 2,n-1
      do k = 1,l-1
        v(i,j,k) = v_star(i,j,k)+dx/Ap_v(i,j,k)*(P_prime(i,j-1,k)-P_prime(i,j,k))*alpha_v
      end do
    end do
  end do

  ! ====================== W-Velocity ====================== !
  ! Correct velocity values
  do i = 1,m-1
    do j = 1,n-1
      do k = 2,l-1
        w(i,j,k) = w_star(i,j,k)+dx/Ap_w(i,j,k)*(P_prime(i,j,k-1)-P_prime(i,j,k))*alpha_v
      end do
    end do
  end do

  return

end subroutine velocity3d_correct
