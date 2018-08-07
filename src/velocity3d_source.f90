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

    ! Calculate interior coefficients
    do i = 2,m-1
      do j = 2,n-2
        do k = 1,l-1

          ! Update convection terms
          Fw = dy*dz*(u_star(i-1,j,k)+u_star(i,j,k))/2
		      Fe = dy*dz*(u_star(i,j,k)+u_star(i+1,j,k))/2
		      Fs = dz*dx*(v_star(i-1,j,k)+v_star(i,j,k))/2
		      Fn = dz*dx*(v_star(i-1,j+1,k)+v_star(i,j+1,k))/2
		      Fb = dx*dy*(w_star(i-1,j,k)+w_star(i,j,k))/2
		      Ft = dx*dy*(w_star(i-1,j,k+1)+w_star(i,j,k+1))/2

          ! Update diffusion terms
          Dw = dy*dz/dx*Pr
          De = dy*dz/dx*Pr
          Ds = dz*dx/dy*Pr
          Dn = dz*dx/dy*Pr
          Db = dx*dy/dz*Pr
          Dt = dx*dy/dz*Pr

		      ! Compute Coefficients - Power Law Differening Scheme
		      Aw_u(i,j,k) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
		      Ae_u(i,j,k) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
		      As_u(i,j,k) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
		      An_u(i,j,k) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)
		      Ab_u(i,j,k) = Db*max(0.0,(1-0.1*abs(Fb/Db))**5)+max(Fb,0.0)
		      At_u(i,j,k) = Dt*max(0.0,(1-0.1*abs(Ft/Dt))**5)+max(-Ft,0.0)

          ! Check bottom / top
          if (k .eq. 1) then
            As_u(i,j,k) = 0.
          elseif (k .eq. l-1) then
            An_u(i,j,k) = 0.
          end if

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

    ! South
    Aw_u(:,1,:) = 0.
    Ae_u(:,1,:) = 0.
    As_u(:,1,:) = 0.
    An_u(:,1,:) = 1.
    Ab_u(:,1,:) = 0.
    At_u(:,1,:) = 0.

    Ap_u(:,1,:) = 1.
    b_u(:,1,:) = 0.

    ! North
    Aw_u(:,n-1,:) = 0.
    Ae_u(:,n-1,:) = 0.
    As_u(:,n-1,:) = 1.
    An_u(:,n-1,:) = 0.
    Ab_u(:,n-1,:) = 0.
    At_u(:,n-1,:) = 0.

    Ap_u(:,n-1,:) = 1.
    b_u(:,n-1,:) = 0.

    ! West-South
    Aw_u(1,1,:) = 0.
    Ae_u(1,1,:) = 1.
    As_u(1,1,:) = 0.
    An_u(1,1,:) = 1.
    Ab_u(1,1,:) = 0.
    At_u(1,1,:) = 0.

    Ap_u(1,1,:) = 2.
    b_u(1,1,:) = 0.

    ! West-North
    Aw_u(1,n-1,:) = 0.
    Ae_u(1,n-1,:) = 1.
    As_u(1,n-1,:) = 1.
    An_u(1,n-1,:) = 0.
    Ab_u(1,n-1,:) = 0.
    At_u(1,n-1,:) = 0.

    Ap_u(1,n-1,:) = 2.
    b_u(1,n-1,:) = 0.

    ! East-South
    Aw_u(m,1,:) = 1.
    Ae_u(m,1,:) = 0.
    As_u(m,1,:) = 0.
    An_u(m,1,:) = 1.
    Ab_u(m,1,:) = 0.
    At_u(m,1,:) = 0.

    Ap_u(m,1,:) = 2.
    b_u(m,1,:) = 0.

    ! East-North
    Aw_u(m,n-1,:) = 1.
    Ae_u(m,n-1,:) = 0.
    As_u(m,n-1,:) = 1.
    An_u(m,n-1,:) = 0.
    Ab_u(m,n-1,:) = 0.
    At_u(m,n-1,:) = 0.

    Ap_u(m,n-1,:) = 2.
    b_u(m,n-1,:) = 0.

    ! West
    Aw_u(1,:,:) = 0.
    Ae_u(1,:,:) = 0.
    As_u(1,:,:) = 0.
    An_u(1,:,:) = 0.
    Ab_u(1,:,:) = 0.
    At_u(1,:,:) = 0.

    Ap_u(1,:,:) = 1.
    b_u(1,:,:) = 0.

    ! East
    Aw_u(m,:,:) = 0.
    Ae_u(m,:,:) = 0.
    As_u(m,:,:) = 0.
    An_u(m,:,:) = 0.
    Ab_u(m,:,:) = 0.
    At_u(m,:,:) = 0.

    Ap_u(m,:,:) = 1.
    b_u(m,:,:) = 0.

  end if

  ! ====================== V-Velocity ====================== !
  ! v-velocity update loop
  if (direction .eq. "v") then

    ! Calculate interior coefficients
    do i = 2,m-2
      do j = 2,n-1
        do k = 1,l-1

          ! Update convection terms
		      Fw = dy*dz*(u_star(i,j-1,k)+u_star(i,j,k))/2
		      Fe = dy*dz*(u_star(i+1,j-1,k)+u_star(i+1,j,k))/2
		      Fs = dz*dx*(v_star(i,j-1,k)+v_star(i,j,k))/2
		      Fn = dz*dx*(v_star(i,j,k)+v_star(i,j+1,k))/2
		      Fb = dx*dy*(w_star(i,j-1,k)+w_star(i,j,k))/2
		      Ft = dx*dy*(w_star(i,j-1,k+1)+w_star(i,j,k+1))/2

          ! Update diffusion terms
          Dw = dy*dz/dx*Pr
          De = dy*dz/dx*Pr
          Ds = dz*dx/dy*Pr
          Dn = dz*dx/dy*Pr
		      Db = dx*dy/dz*Pr
		      Dt = dx*dy/dz*Pr

		      ! Compute Coefficients - Power Law Differening Scheme
		      Aw_v(i,j,k) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
		      Ae_v(i,j,k) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
		      As_v(i,j,k) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
		      An_v(i,j,k) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)
		      Ab_v(i,j,k) = Db*max(0.0,(1-0.1*abs(Fb/Db))**5)+max(Fb,0.0)
		      At_v(i,j,k) = Dt*max(0.0,(1-0.1*abs(Ft/Dt))**5)+max(-Ft,0.0)

          ! Check bottom / top
          if (k .eq. 1) then
            As_v(i,j,k) = 0.
          elseif (k .eq. l-1) then
            An_v(i,j,k) = 0.
          end if

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

    ! West
    Aw_v(1,:,:) = 0.
    Ae_v(1,:,:) = 1.
    As_v(1,:,:) = 0.
    An_v(1,:,:) = 0.
    Ab_v(1,:,:) = 0.
    At_v(1,:,:) = 0.

    Ap_v(1,:,:) = 1.
    b_v(1,:,:) = 0.

    ! East
    Aw_v(m-1,:,:) = 1.
    Ae_v(m-1,:,:) = 0.
    As_v(m-1,:,:) = 0.
    An_v(m-1,:,:) = 0.
    Ab_v(m-1,:,:) = 0.
    At_v(m-1,:,:) = 0.

    Ap_v(m-1,:,:) = 1.
    b_v(m-1,:,:) = 0.

    ! West-South
    Aw_v(1,1,:) = 0.
    Ae_v(1,1,:) = 1.
    As_v(1,1,:) = 0.
    An_v(1,1,:) = 1.
    Ab_v(1,1,:) = 0.
    At_v(1,1,:) = 0.

    Ap_v(1,1,:) = 2.
    b_v(1,1,:) = 0.

    ! West-North
    Aw_v(1,n,:) = 0.
    Ae_v(1,n,:) = 1.
    As_v(1,n,:) = 1.
    An_v(1,n,:) = 0.
    Ab_v(1,n,:) = 0.
    At_v(1,n,:) = 0.

    Ap_v(1,n,:) = 2.
    b_v(1,n,:) = 0.

    ! East-South
    Aw_v(1,n,:) = 1.
    Ae_v(1,n,:) = 0.
    As_v(1,n,:) = 0.
    An_v(1,n,:) = 1.
    Ab_v(1,n,:) = 0.
    At_v(1,n,:) = 0.

    Ap_v(1,n,:) = 2.
    b_v(1,n,:) = 0.

    ! East-North
    Aw_v(m-1,n,:) = 1.
    Ae_v(m-1,n,:) = 0.
    As_v(m-1,n,:) = 1.
    An_v(m-1,n,:) = 0.
    Ab_v(m-1,n,:) = 0.
    At_v(m-1,n,:) = 0.

    Ap_v(m-1,n,:) = 2.
    b_v(m-1,n,:) = 0.

    ! South
    Aw_v(:,1,:) = 0.
    Ae_v(:,1,:) = 0.
    As_v(:,1,:) = 0.
    An_v(:,1,:) = 0.
    Ab_v(:,1,:) = 0.
    At_v(:,1,:) = 0.

    Ap_v(:,1,:) = 1.
    b_v(:,1,:) = 0.

    ! North
    Aw_v(:,n,:) = 0.
    Ae_v(:,n,:) = 0.
    As_v(:,n,:) = 0.
    An_v(:,n,:) = 0.
    Ab_v(:,n,:) = 0.
    At_v(:,n,:) = 0.

    Ap_v(:,n,:) = 1.
    b_v(:,n,:) = 0.

  end if

  ! ====================== W-Velocity ====================== !
  ! w-velocity update loop
  if (direction .eq. "w") then

    ! Calculate interior coefficients
    do i = 2,m-2
      do j = 2,n-2
        do k = 2,l-1

          ! Update convection terms
		      Fw = dy*dz*(u_star(i,j,k-1)+u_star(i,j,k))/2
		      Fe = dy*dz*(u_star(i+1,j,k-1)+u_star(i+1,j,k))/2
		      Fs = dz*dx*(v_star(i,j,k-1)+v_star(i,j,k))/2
		      Fn = dz*dx*(v_star(i,j+1,k-1)+v_star(i,j+1,k))/2
		      Fb = dx*dy*(w_star(i,j,k-1)+w_star(i,j,k))/2
		      Ft = dx*dy*(w_star(i,j,k)+w_star(i,j,k+1))/2

          ! Update diffusion terms
          Dw = dy*dz/dx*Pr
          De = dy*dz/dx*Pr
          Ds = dz*dx/dy*Pr
          Dn = dz*dx/dy*Pr
		      Db = dx*dy/dz*Pr
		      Dt = dx*dy/dz*Pr

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
		      b_w(i,j,k) = Su_w(i,j,k)*dx*dy*dz+Ra*Pr*(((T(i,j,k)+T(i,j,k-1))/2.0)-0.5)*dx*dy*dz

        end do
      end do
    end do

    ! West
    Aw_w(1,:,:) = 0.
    Ae_w(1,:,:) = 1.
    As_w(1,:,:) = 0.
    An_w(1,:,:) = 0.
    Ab_w(1,:,:) = 0.
    At_w(1,:,:) = 0.

    Ap_w(1,:,:) = 1.
    b_w(1,:,:) = 0.

    ! East
    Aw_w(m-1,:,:) = 1.
    Ae_w(m-1,:,:) = 0.
    As_w(m-1,:,:) = 0.
    An_w(m-1,:,:) = 0.
    Ab_w(m-1,:,:) = 0.
    At_w(m-1,:,:) = 0.

    Ap_w(m-1,:,:) = 1.
    b_w(m-1,:,:) = 0.

    ! South
    Aw_w(:,1,:) = 0.
    Ae_w(:,1,:) = 0.
    As_w(:,1,:) = 0.
    An_w(:,1,:) = 1.
    Ab_w(:,1,:) = 0.
    At_w(:,1,:) = 0.

    Ap_w(:,1,:) = 1.
    b_w(:,1,:) = 0.

    ! North
    Aw_w(:,n-1,:) = 0.
    Ae_w(:,n-1,:) = 0.
    As_w(:,n-1,:) = 1.
    An_w(:,n-1,:) = 0.
    Ab_w(:,n-1,:) = 0.
    At_w(:,n-1,:) = 0.

    Ap_w(:,n-1,:) = 1.
    b_w(:,n-1,:) = 0.

    ! West-South
    Aw_w(1,1,:) = 0.
    Ae_w(1,1,:) = 1.
    As_w(1,1,:) = 0.
    An_w(1,1,:) = 1.
    Ab_w(1,1,:) = 0.
    At_w(1,1,:) = 0.

    Ap_w(1,1,:) = 2.
    b_w(1,1,:) = 0.

    ! West-North
    Aw_w(1,n-1,:) = 0.
    Ae_w(1,n-1,:) = 1.
    As_w(1,n-1,:) = 1.
    An_w(1,n-1,:) = 0.
    Ab_w(1,n-1,:) = 0.
    At_w(1,n-1,:) = 0.

    Ap_w(1,n-1,:) = 2.
    b_w(1,n-1,:) = 0.

    ! East-South
    Aw_w(m-1,1,:) = 1.
    Ae_w(m-1,1,:) = 0.
    As_w(m-1,1,:) = 0.
    An_w(m-1,1,:) = 1.
    Ab_w(m-1,1,:) = 0.
    At_w(m-1,1,:) = 0.

    Ap_w(m-1,1,:) = 2.
    b_w(m-1,1,:) = 0.

    ! East-North
    Aw_w(m-1,n-1,:) = 1.
    Ae_w(m-1,n-1,:) = 0.
    As_w(m-1,n-1,:) = 0.
    An_w(m-1,n-1,:) = 1.
    Ab_w(m-1,n-1,:) = 0.
    At_w(m-1,n-1,:) = 0.

    Ap_w(m-1,n-1,:) = 2.
    b_w(m-1,n-1,:) = 0.

    ! Bottom
    Aw_w(:,:,1) = 0.
    Ae_w(:,:,1) = 0.
    As_w(:,:,1) = 0.
    An_w(:,:,1) = 0.
    Ab_w(:,:,1) = 0.
    At_w(:,:,1) = 0.

    Ap_w(:,:,1) = 1.
    b_w(:,:,1) = 0.

    ! Top
    Aw_w(:,:,l) = 0.
    Ae_w(:,:,l) = 0.
    As_w(:,:,l) = 0.
    An_w(:,:,l) = 0.
    Ab_w(:,:,l) = 0.
    At_w(:,:,l) = 0.

    Ap_w(:,:,l) = 1.
    b_w(:,:,l) = 0.

  end if

  return

end subroutine velocity3d_source
