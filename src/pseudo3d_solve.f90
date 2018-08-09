! pseudo3d_solve Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine computes the pseudo velocities (u_hat, v_hat, and w_hat for a 3D
! CFD problem.
!

subroutine pseudo3d_solve

  implicit none

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! ========================== u_hat ========================== !

  do i = 2, m-1
    do j = 2, n-2
      do k = 1, l-1

        if (k .eq. 1) then
          u_hat(i,j,k) = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
                        Ae_u(i,j,k)*u_star(i+1,j,k)+ &
                        As_u(i,j,k)*u_star(i,j-1,k)+ &
                        An_u(i,j,k)*u_star(i,j+1,k)+ &
                        At_u(i,j,k)*u_star(i,j,k+1)+ &
                        b_u(i,j,k))/Ap_u(i,j,k)
        elseif (k .eq. l-1) then
          u_hat(i,j,k) = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
                        Ae_u(i,j,k)*u_star(i+1,j,k)+ &
                        As_u(i,j,k)*u_star(i,j-1,k)+ &
                        An_u(i,j,k)*u_star(i,j+1,k)+ &
                        Ab_u(i,j,k)*u_star(i,j,k-1)+ &
                        b_u(i,j,k))/Ap_u(i,j,k)
        else
          u_hat(i,j,k) = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
                        Ae_u(i,j,k)*u_star(i+1,j,k)+ &
                        As_u(i,j,k)*u_star(i,j-1,k)+ &
                        An_u(i,j,k)*u_star(i,j+1,k)+ &
                        Ab_u(i,j,k)*u_star(i,j,k-1)+ &
                        At_u(i,j,k)*u_star(i,j,k+1)+ &
                        b_u(i,j,k))/Ap_u(i,j,k)
        end if

      end do
    end do
  end do

  u_hat(1,:,:) = u_hat(2,:,:)
  u_hat(m,:,:) = u_hat(m-1,:,:)
  u_hat(:,1,:) = u_hat(:,2,:)
  u_hat(:,n-1,:) = u_hat(:,n-2,:)

  ! ========================== v_hat ========================== !

  do i = 2, m-2
    do j = 2, n-1
      do k = 1, l-1

        if (k .eq. 1) then
          v_hat(i,j,k) = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
                          Ae_v(i,j,k)*v_star(i+1,j,k)+ &
                          As_v(i,j,k)*v_star(i,j-1,k)+ &
                          An_v(i,j,k)*v_star(i,j+1,k)+ &
                          At_v(i,j,k)*v_star(i,j,k+1)+ &
                          b_v(i,j,k))/Ap_v(i,j,k)
        elseif (k .eq. l-1) then
          v_hat(i,j,k) = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
                          Ae_v(i,j,k)*v_star(i+1,j,k)+ &
                          As_v(i,j,k)*v_star(i,j-1,k)+ &
                          An_v(i,j,k)*v_star(i,j+1,k)+ &
                          Ab_v(i,j,k)*v_star(i,j,k-1)+ &
                          b_v(i,j,k))/Ap_v(i,j,k)
        else
          v_hat(i,j,k) = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
                          Ae_v(i,j,k)*v_star(i+1,j,k)+ &
                          As_v(i,j,k)*v_star(i,j-1,k)+ &
                          An_v(i,j,k)*v_star(i,j+1,k)+ &
                          Ab_v(i,j,k)*v_star(i,j,k-1)+ &
                          At_v(i,j,k)*v_star(i,j,k+1)+ &
                          b_v(i,j,k))/Ap_v(i,j,k)
        end if

        v_hat(1,:,:) = v_hat(2,:,:)
        v_hat(m-1,:,:) = v_hat(m-2,:,:)
        v_hat(:,1,:) = v_hat(:,2,:)
        v_hat(:,n,:) = v_hat(:,n-1,:)

      end do
    end do
  end do

  ! ========================== w_hat ========================== !

  do i = 2, m-2
    do j = 2, n-2
      do k = 2, l-1

			  w_hat(i,j,k) = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
                        Ae_w(i,j,k)*w_star(i+1,j,k)+ &
                        As_w(i,j,k)*w_star(i,j-1,k)+ &
                        An_w(i,j,k)*w_star(i,j+1,k)+ &
                        Ab_w(i,j,k)*w_star(i,j,k-1)+ &
                        At_w(i,j,k)*w_star(i,j,k+1)+ &
                        b_w(i,j,k))/Ap_w(i,j,k)

	    end do
    end do
  end do

  w_hat(1,:,:) = w_hat(2,:,:)
  w_hat(m-1,:,:) = w_hat(m-2,:,:)
  w_hat(:,1,:) = w_hat(:,2,:)
  w_hat(:,n-1,:) = w_hat(:,n-2,:)
  w_hat(:,:,1) = 0.
  w_hat(:,:,l) = 0.

  return

end subroutine pseudo3d_solve
