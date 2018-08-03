! velocity3d_source Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine updates the source terms for the solution of the momentum
! equations in the SIMPLER algorithm.
!

subroutine velocity3d_source(direction)

  ! Include standard variable header
  include "var3d.dec"

  ! Define input variables
  character :: direction

  ! Define internal variables
  integer :: i, j, k
  real(8) :: Fw, Fe, Fs, Fn, Fb, Ft, Dw, De, Ds, Dn, Db, Dt

  ! ====================== U-Velocity ====================== !
  ! u-velocity update loop
  if (direction .eq. "u") then

    ! Reset velocity b values
    !b_u = 0.

    ! Calculate interior coefficients
    do i = 2,m-1
      do j = 2,n-2
        do k = 2,l-2

          ! Update convection terms
          Fw = dy*dz*(u_star(i-1,j,k)+u_star(i,j,k))/2
		      Fe = dy*dz*(u_star(i,j,k)+u_star(i+1,j,k))/2
		      Fs = dz*dx*(v_star(i-1,j,k)+v_star(i,j,k))/2
		      Fn = dz*dx*(v_star(i-1,j+1,k)+v_star(i,j+1,k))/2
		      Fb = dx*dy*(w_star(i-1,j,k)+w_star(i,j,k))/2
		      Ft = dx*dy*(w_star(i-1,j,k+1)+w_star(i,j,k+1))/2

          ! Update diffusion terms
          Dw = dy*dz/dx*Pr*(Pr/Ra)**(0.5)
          De = dy*dz/dx*Pr*(Pr/Ra)**(0.5)
          Ds = dz*dx/dy*Pr*(Pr/Ra)**(0.5)
          Dn = dz*dx/dy*Pr*(Pr/Ra)**(0.5)
          Db = dx*dy/dz*Pr*(Pr/Ra)**(0.5)
          Dt = dx*dy/dz*Pr*(Pr/Ra)**(0.5)

		      ! Compute Coefficients - Power Law Differening Scheme
		      Aw_u(i,j,k) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
		      Ae_u(i,j,k) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
		      As_u(i,j,k) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
		      An_u(i,j,k) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)
		      Ab_u(i,j,k) = Db*max(0.0,(1-0.1*abs(Fb/Db))**5)+max(Fb,0.0)
		      At_u(i,j,k) = Dt*max(0.0,(1-0.1*abs(Ft/Dt))**5)+max(-Ft,0.0)

		      ! Update Ap coefficient
		      Ap_u(i,j,k) = Aw_u(i,j,k)+Ae_u(i,j,k)+As_u(i,j,k)+An_u(i,j,k)+Ab_u(i,j,k)+At_u(i,j,k)+Sp_u(i,j,k)*dx*dy*dz

		      ! Check False Diffusion
		      if (Ap_u(i,j,k) .eq. 0.) then
            print *, "False Diffusion @:", i, j, k
            Ap_u(i,j,k) = 1.
          end if

		      ! Update b values
		      b_u(i,j,k) = Su_u(i,j,k)*dx*dy*dz

	     end do
      end do
    end do

  end if

  ! ====================== V-Velocity ====================== !
  ! v-velocity update loop
  if (direction .eq. "v") then

    ! Reset velocity b values
    !b_v = 0.

    ! Calculate interior coefficients
    do i = 2,m-2
      do j = 2,n-1
        do k = 2,l-2

          ! Update convection terms
		      Fw = dy*dz*(u_star(i,j-1,k)+u_star(i,j,k))/2
		      Fe = dy*dz*(u_star(i+1,j,k)+u_star(i+1,j+1,k))/2
		      Fs = dz*dx*(v_star(i,j-1,k)+v_star(i,j,k))/2
		      Fn = dz*dx*(v_star(i,j,k)+v_star(i,j+1,k))/2
		      Fb = dx*dy*(w_star(i,j-1,k)+w_star(i,j,k))/2
		      Ft = dx*dy*(w_star(i,j,k+1)+w_star(i,j+1,k+1))/2

          ! Update diffusion terms
          Dw = dy*dz/dx*Pr*(Pr/Ra)**(0.5)
          De = dy*dz/dx*Pr*(Pr/Ra)**(0.5)
          Ds = dz*dx/dy*Pr*(Pr/Ra)**(0.5)
          Dn = dz*dx/dy*Pr*(Pr/Ra)**(0.5)
		      Db = dx*dy/dz*Pr*(Pr/Ra)**(0.5)
		      Dt = dx*dy/dz*Pr*(Pr/Ra)**(0.5)

		      ! Compute Coefficients - Power Law Differening Scheme
		      Aw_v(i,j,k) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
		      Ae_v(i,j,k) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
		      As_v(i,j,k) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
		      An_v(i,j,k) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)
		      Ab_v(i,j,k) = Db*max(0.0,(1-0.1*abs(Fb/Db))**5)+max(Fb,0.0)
		      At_v(i,j,k) = Dt*max(0.0,(1-0.1*abs(Ft/Dt))**5)+max(-Ft,0.0)

		      ! Update Ap coefficient
		      Ap_v(i,j,k) = Aw_v(i,j,k)+Ae_v(i,j,k)+As_v(i,j,k)+An_v(i,j,k)+Ab_v(i,j,k)+At_v(i,j,k)+Sp_v(i,j,k)*dx*dy*dz

		      ! Check False Diffusion
		      if (Ap_v(i,j,k) .eq. 0.) then
            print *, "False Diffusion @:", i, j, k
            Ap_v(i,j,k) = 1.
          end if

          ! Update b values
          b_v(i,j,k) = Su_v(i,j,k)*dx*dy*dz

        end do
      end do
    end do


  end if

  ! ====================== W-Velocity ====================== !
  ! w-velocity update loop
  if (direction .eq. "w") then

    ! Reset velocity b values
    !b_w = 0.

    ! Calculate interior coefficients
    do i = 2,m-2
      do j = 2,n-2
        do k = 2,l-1

          ! Update convection terms
		      Fw = dy*dz*(u_star(i,j,k-1)+u_star(i,j,k))/2
		      Fe = dy*dz*(u_star(i+1,j,k)+u_star(i+1,j,k+1))/2
		      Fs = dz*dx*(v_star(i,j,k-1)+v_star(i,j,k))/2
		      Fn = dz*dx*(v_star(i,j+1,k)+v_star(i,j+1,k+1))/2
		      Fb = dx*dy*(w_star(i,j,k-1)+w_star(i,j,k))/2
		      Ft = dx*dy*(w_star(i,j,k)+w_star(i,j,k+1))/2

          ! Update diffusion terms
          Dw = dy*dz/dx*Pr*(Pr/Ra)**(0.5)
          De = dy*dz/dx*Pr*(Pr/Ra)**(0.5)
          Ds = dz*dx/dy*Pr*(Pr/Ra)**(0.5)
          Dn = dz*dx/dy*Pr*(Pr/Ra)**(0.5)
		      Db = dx*dy/dz*Pr*(Pr/Ra)**(0.5)
		      Dt = dx*dy/dz*Pr*(Pr/Ra)**(0.5)

		      ! Compute Coefficients - Power Law Differening Scheme
		      Aw_w(i,j,k) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
		      Ae_w(i,j,k) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
		      As_w(i,j,k) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
		      An_w(i,j,k) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)
		      Ab_w(i,j,k) = Db*max(0.0,(1-0.1*abs(Fb/Db))**5)+max(Fb,0.0)
		      At_w(i,j,k) = Dt*max(0.0,(1-0.1*abs(Ft/Dt))**5)+max(-Ft,0.0)

		      ! Update Ap coefficient
		      Ap_w(i,j,k) = Aw_w(i,j,k)+Ae_w(i,j,k)+As_w(i,j,k)+An_w(i,j,k)+Ab_w(i,j,k)+At_w(i,j,k)+Sp_w(i,j,k)*dx*dy*dz

		      ! Check False Diffusion
		      if (Ap_w(i,j,k) .eq. 0.) then
            print *, "False Diffusion @:", i, j, k
			      Ap_w(i,j,k) = 1.
          end if

		      ! Update b values
		      b_w(i,j,k) = Su_w(i,j,k)*dx*dy*dz+Pr*(((T(i,j,k)+T(i,j,k-1))/2.0))*dx*dy*dz-(Pr/2.0)*dx*dy*dz

        end do
      end do
    end do


  end if

  return

end subroutine velocity3d_source
