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
    do j = 2,n-2
      do k = 1,l-1

        u(i,j,k) = u_star(i,j,k)+dy/Ap_u(i,j,k)*(P_prime(i-1,j,k)-P_prime(i,j,k))*alpha_v

      end do
    end do
  end do

  u(1,:,:) = u(2,:,:)
  u(m,:,:) = u(m-1,:,:)
  u(:,1,:) = u(:,2,:)
  u(:,n-1,:) = u(:,n-2,:)

  ! ====================== V-Velocity ====================== !
  ! Correct velocity values
  do i = 2,m-2
    do j = 2,n-1
      do k = 1,l-1
        v(i,j,k) = v_star(i,j,k)+dx/Ap_v(i,j,k)*(P_prime(i,j-1,k)-P_prime(i,j,k))*alpha_v
      end do
    end do
  end do

  v(1,:,:) = v(2,:,:)
  v(m-1,:,:) = v(m-2,:,:)
  v(:,1,:) = v(:,2,:)
  v(:,n,:) = v(:,n-1,:)

  ! ====================== W-Velocity ====================== !
  ! Correct velocity values
  do i = 2,m-2
    do j = 2,n-2
      do k = 2,l-1
        w(i,j,k) = w_star(i,j,k)+dx/Ap_w(i,j,k)*(P_prime(i,j,k-1)-P_prime(i,j,k))*alpha_v
      end do
    end do
  end do

  w(1,:,:) = w(2,:,:)
  w(m-1,:,:) = w(m-2,:,:)
  w(:,1,:) = w(:,2,:)
  w(:,n-1,:) = w(:,n-2,:)
  w(:,:,1) = 0.
  w(:,:,l) = 0.

  return

end subroutine velocity3d_correct
