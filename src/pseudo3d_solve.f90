! pseudo3d_solve Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine computes the pseudo velocities (u_hat, v_hat, and w_hat) for a 3D
! CFD problem.
!

subroutine pseudo3d_solve

  implicit none

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! ========================== u_hat ========================== !

  do i = 1, m
    do j = 1, n-1
	  do k = 1, l-1

	    if (i .eq. 1) then
		  if (j .eq. 1) then
		    if (k .eq. 1) then
			  u_hat = (Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else if (k .eq. l-1) then
			  u_hat = (Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else
			  u_hat = (Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			end if
		  elseif (j .eq. n-1) then
		    if (k .eq. 1) then
			  u_hat = (Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else if (k .eq. l-1)
			  u_hat = (Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else
			  u_hat = (Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			end if
		  else
		    if (k .eq. 1) then
			  u_hat = (Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else if (k .eq. l-1) then
			  u_hat = (Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else
			  u_hat = (Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			end if
		  end if
		elseif (i .eq. m) then
		  if (j .eq. 1) then
		    if (k .eq. 1) then
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           An_u(i,j,k)*u_star(i,j+1,k)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else if (k .eq. l-1)
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			end if
		  elseif (j .eq. n-1) then
		    if (k .eq. 1) then
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           As_u(i,j,k)*u_star(i,j-1,k)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else if (k .eq. l-1) then
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           As_u(i,j,k)*u_star(i,j-1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           As_u(i,j,k)*u_star(i,j-1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			end if
		  else
		    if (k .eq. 1) then
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           As_u(i,j,k)*u_star(i,j-1,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else if (k .eq. l-1) then
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           As_u(i,j,k)*u_star(i,j-1,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           As_u(i,j,k)*u_star(i,j-1,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			end if
		  end if
		else
		  if (j .eq. 1) then
		    if (k .eq. 1) then
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else if (k .eq. l-1) then
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			end if
		  elseif (j .eq. n-1) then
		    if (k .eq. 1) then
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else if (k .eq. l-1) then
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			end if
		  else
		    if (k .eq. 1) then
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else if (k .eq. l-1) then
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			else
			  u_hat = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
			           Ae_u(i,j,k)*u_star(i+1,j,k)+ &
					   As_u(i,j,k)*u_star(i,j-1,k)+ &
					   An_u(i,j,k)*u_star(i,j+1,k)+ &
					   Ab_u(i,j,k)*u_star(i,j,k-1)+ &
					   At_u(i,j,k)*u_star(i,j,k+1)+ &
					   b_u(i,j,k))/Ap_u(i,j,k)
			end if
		  end if
		end if

	  end do
	end do
  end do

  ! ========================== v_hat ========================== !

  do i = 1, m-1
    do j = 1, n
	  do k = 1, l-1

	    if (i .eq. 1) then
		  if (j .eq. 1) then
		    if (k .eq. 1) then
			  v_hat = (Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else if (k .eq. l-1) then
			  v_hat = (Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else
			  v_hat = (Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			end if
		  elseif (j .eq. n) then
		    if (k .eq. 1) then
			  v_hat = (Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else if (k .eq. l-1)
			  v_hat = (Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else
			  v_hat = (Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			end if
		  else
		    if (k .eq. 1) then
			  v_hat = (Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else if (k .eq. l-1) then
			  v_hat = (Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else
			  v_hat = (Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			end if
		  end if
		elseif (i .eq. m-1) then
		  if (j .eq. 1) then
		    if (k .eq. 1) then
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           An_v(i,j,k)*v_star(i,j+1,k)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else if (k .eq. l-1)
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			end if
		  elseif (j .eq. n) then
		    if (k .eq. 1) then
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           As_v(i,j,k)*v_star(i,j-1,k)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else if (k .eq. l-1) then
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           As_v(i,j,k)*v_star(i,j-1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           As_v(i,j,k)*v_star(i,j-1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			end if
		  else
		    if (k .eq. 1) then
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           As_v(i,j,k)*v_star(i,j-1,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else if (k .eq. l-1) then
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           As_v(i,j,k)*v_star(i,j-1,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           As_v(i,j,k)*v_star(i,j-1,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			end if
		  end if
		else
		  if (j .eq. 1) then
		    if (k .eq. 1) then
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else if (k .eq. l-1) then
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			end if
		  elseif (j .eq. n) then
		    if (k .eq. 1) then
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else if (k .eq. l-1) then
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			end if
		  else
		    if (k .eq. 1) then
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else if (k .eq. l-1) then
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			else
			  v_hat = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
			           Ae_v(i,j,k)*v_star(i+1,j,k)+ &
					   As_v(i,j,k)*v_star(i,j-1,k)+ &
					   An_v(i,j,k)*v_star(i,j+1,k)+ &
					   Ab_v(i,j,k)*v_star(i,j,k-1)+ &
					   At_v(i,j,k)*v_star(i,j,k+1)+ &
					   b_v(i,j,k))/Ap_v(i,j,k)
			end if
		  end if
		end if

	  end do
	end do
  end do

  ! ========================== w_hat ========================== !

  do i = 1, m-1
    do j = 1, n-1
	  do k = 1, l

	    if (i .eq. 1) then
		  if (j .eq. 1) then
		    if (k .eq. 1) then
			  w_hat = (Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else if (k .eq. l) then
			  w_hat = (Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else
			  w_hat = (Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			end if
		  elseif (j .eq. n-1) then
		    if (k .eq. 1) then
			  w_hat = (Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else if (k .eq. l)
			  w_hat = (Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else
			  w_hat = (Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			end if
		  else
		    if (k .eq. 1) then
			  w_hat = (Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else if (k .eq. l) then
			  w_hat = (Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else
			  w_hat = (Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			end if
		  end if
		elseif (i .eq. m-1) then
		  if (j .eq. 1) then
		    if (k .eq. 1) then
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           An_w(i,j,k)*w_star(i,j+1,k)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else if (k .eq. l)
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			end if
		  elseif (j .eq. n-1) then
		    if (k .eq. 1) then
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           As_w(i,j,k)*w_star(i,j-1,k)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else if (k .eq. l) then
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           As_w(i,j,k)*w_star(i,j-1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           As_w(i,j,k)*w_star(i,j-1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			end if
		  else
		    if (k .eq. 1) then
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           As_w(i,j,k)*w_star(i,j-1,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else if (k .eq. l) then
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           As_w(i,j,k)*w_star(i,j-1,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           As_w(i,j,k)*w_star(i,j-1,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			end if
		  end if
		else
		  if (j .eq. 1) then
		    if (k .eq. 1) then
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else if (k .eq. l) then
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			end if
		  elseif (j .eq. n-1) then
		    if (k .eq. 1) then
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else if (k .eq. l) then
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			end if
		  else
		    if (k .eq. 1) then
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else if (k .eq. l) then
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			else
			  w_hat = (Aw_w(i,j,k)*w_star(i-1,j,k)+ &
			           Ae_w(i,j,k)*w_star(i+1,j,k)+ &
					   As_w(i,j,k)*w_star(i,j-1,k)+ &
					   An_w(i,j,k)*w_star(i,j+1,k)+ &
					   Ab_w(i,j,k)*w_star(i,j,k-1)+ &
					   At_w(i,j,k)*w_star(i,j,k+1)+ &
					   b_w(i,j,k))/Ap_w(i,j,k)
			end if
		  end if
		end if

	  end do
	end do
  end do

  return

end subroutine pseudo3d_solve
