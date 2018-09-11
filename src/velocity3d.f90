! velocity3d subroutines
!
! Written by Matt Blomquist
! Last Update: 2018-09-07 (YYYY-MM-DD)
!
! The subroutines below serve to calculate the coefficients and source terms
! for the u, v, and w velocities in a 3D natural convection code. Additionally,
! the velocity solve, pseudo solve, and velocity correct subroutines are also
! included.
!
subroutine velocity3d_source

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k
  real(8) :: Dw, De, Ds, Dn, Db, Dt
  real(8) :: Fw, Fe, Fs, Fn, Fb, Ft

  ! Initlaize Coefficients
  Aw_u = 0.
  Ae_u = 0.
  As_u = 0.
  An_u = 0.
  Ab_u = 0.
  At_u = 0.
  Ap_u = 0.
  b_u = 0.

  Aw_v = 0.
  Ae_v = 0.
  As_v = 0.
  An_v = 0.
  Ab_v = 0.
  At_v = 0.
  Ap_v = 0.
  b_v = 0.

  Aw_w = 0.
  Ae_w = 0.
  As_w = 0.
  An_w = 0.
  Ab_w = 0.
  At_w = 0.
  Ap_w = 0.
  b_w = 0.

  ! ====================================== !
  ! ====================================== !
  ! Calculate Coefficients for U-Velocity
  do i = 2, m-1
    do j = 2, n-2
	     do k = 2, l-2

	    ! Calculate diffusion terms
		Dw = dy*dz/dx
		De = dy*dz/dx
		Ds = dz*dx/dy
		Dn = dz*dx/dy
		Db = dx*dy/dz
		Dt = dx*dy/dz

		! Calculate adconvection terms
		Fw = dy*dz*(u(i-1,j,k)+u(i,j,k))/2.
		Fe = dy*dz*(u(i,j,k)+u(i+1,j,k))/2.
		Fs = dz*dx*(v(i-1,j,k)+v(i,j,k))/2.
		Fn = dz*dx*(v(i-1,j+1,k)+v(i,j+1,k))/2.
		Fb = dx*dy*(w(i-1,j,k)+w(i,j,k))/2.
		Ft = dx*dy*(w(i-1,j,k+1)+w(i,j,k+1))/2.

		! Calculate Coefficients :: Power Law Differencing Scheme
		Aw_u(i,j,k) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
		Ae_u(i,j,k) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
		As_u(i,j,k) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
		An_u(i,j,k) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)
		Ab_u(i,j,k) = Db*max(0.0,(1-0.1*abs(Fb/Db))**5)+max(Fb,0.0)
		At_u(i,j,k) = Dt*max(0.0,(1-0.1*abs(Ft/Dt))**5)+max(-Ft,0.0)

		! Update Ap_u Coefficient
		Ap_u(i,j,k) = Aw_u(i,j,k)+Ae_u(i,j,k)+As_u(i,j,k)+An_u(i,j,k)+Ab_u(i,j,k)+At_u(i,j,k)

		! Calculate Source Term
		b_u(i,j,k) = 0.

	  end do
	end do
  end do

  ! Update Boundary Conditions :: West // Symmetry Condition
  Aw_u(1,:,2:l-2) = 0.
  Ae_u(1,:,2:l-2) = 1.
  As_u(1,:,2:l-2) = 0.
  An_u(1,:,2:l-2) = 0.
  Ab_u(1,:,2:l-2) = 0.
  At_u(1,:,2:l-2) = 0.
  Ap_u(1,:,2:l-2) = 1.
  b_u(1,:,2:l-2) = 0.

  ! Update Boundary Conditions :: East // Symmetry Condition
  Aw_u(m,:,2:l-2) = 1.
  Ae_u(m,:,2:l-2) = 0.
  As_u(m,:,2:l-2) = 0.
  An_u(m,:,2:l-2) = 0.
  Ab_u(m,:,2:l-2) = 0.
  At_u(m,:,2:l-2) = 0.
  Ap_u(m,:,2:l-2) = 1.
  b_u(m,:,2:l-2) = 0.

  ! Update Boundary Conditions :: South // Symmetry Condition
  Aw_u(2:m-1,1,2:l-2) = 0.
  Ae_u(2:m-1,1,2:l-2) = 0.
  As_u(2:m-1,1,2:l-2) = 0.
  An_u(2:m-1,1,2:l-2) = 1.
  Ab_u(2:m-1,1,2:l-2) = 0.
  At_u(2:m-1,1,2:l-2) = 0.
  Ap_u(2:m-1,1,2:l-2) = 1.
  b_u(2:m-1,1,2:l-2) = 0.

  ! Update Boundary Conditions :: North // Symmetry Condition
  Aw_u(2:m-1,n-1,2:l-2) = 0.
  Ae_u(2:m-1,n-1,2:l-2) = 0.
  As_u(2:m-1,n-1,2:l-2) = 1.
  An_u(2:m-1,n-1,2:l-2) = 0.
  Ab_u(2:m-1,n-1,2:l-2) = 0.
  At_u(2:m-1,n-1,2:l-2) = 0.
  Ap_u(2:m-1,n-1,2:l-2) = 1.
  b_u(2:m-1,n-1,2:l-2) = 0.

  ! Update Boundary Conditions :: Bottom // Wall Condition Test with 0 Velocity
  Aw_u(:,:,1) = 0.
  Ae_u(:,:,1) = 0.
  As_u(:,:,1) = 0.
  An_u(:,:,1) = 0.
  Ab_u(:,:,1) = 0.
  At_u(:,:,1) = 0.
  Ap_u(:,:,1) = 1.
  b_u(:,:,1) = 0.

  ! Update Boundary Conditions :: Top // Wall Condition Test with 0 Velocity
  Aw_u(:,:,l-1) = 0.
  Ae_u(:,:,l-1) = 0.
  As_u(:,:,l-1) = 0.
  An_u(:,:,l-1) = 0.
  Ab_u(:,:,l-1) = 0.
  At_u(:,:,l-1) = 0.
  Ap_u(:,:,l-1) = 1.
  b_u(:,:,l-1) = 0.

  ! ====================================== !
  ! ====================================== !

  ! ====================================== !
  ! ====================================== !
  ! Calculate Coefficients for V-Velocity
  do i = 2, m-2
    do j = 2, n-1
	  do k = 2, l-2

	    ! Calculate diffusion terms
		Dw = dy*dz/dx
		De = dy*dz/dx
		Ds = dz*dx/dy
		Dn = dz*dx/dy
		Db = dx*dy/dz
		Dt = dx*dy/dz

		! Calculate adconvection terms
		Fw = dy*dz*(u(i,j-1,k)+u(i,j,k))/2.
		Fe = dy*dz*(u(i+1,j-1,k)+u(i+1,j,k))/2.
		Fs = dz*dx*(v(i,j-1,k)+v(i,j,k))/2.
		Fn = dz*dx*(v(i,j,k)+v(i,j+1,k))/2.
		Fb = dx*dy*(w(i,j-1,k)+w(i,j,k))/2.
		Ft = dx*dy*(w(i,j-1,k+1)+w(i,j,k+1))/2.

		! Calculate Coefficients :: Power Law Differencing Scheme
		Aw_v(i,j,k) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
		Ae_v(i,j,k) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
		As_v(i,j,k) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
		An_v(i,j,k) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)
		Ab_v(i,j,k) = Db*max(0.0,(1-0.1*abs(Fb/Db))**5)+max(Fb,0.0)
		At_v(i,j,k) = Dt*max(0.0,(1-0.1*abs(Ft/Dt))**5)+max(-Ft,0.0)

		! Update Ap_u Coefficient
		Ap_v(i,j,k) = Aw_v(i,j,k)+Ae_v(i,j,k)+As_v(i,j,k)+An_v(i,j,k)+Ab_v(i,j,k)+At_v(i,j,k)

		! Calculate Source Term
		b_v(i,j,k) = 0.

	  end do
	end do
  end do

  ! Update Boundary Conditions :: West // Symmetry Condition
  Aw_v(1,:,2:l-2) = 0.
  Ae_v(1,:,2:l-2) = 1.
  As_v(1,:,2:l-2) = 0.
  An_v(1,:,2:l-2) = 0.
  Ab_v(1,:,2:l-2) = 0.
  At_v(1,:,2:l-2) = 0.
  Ap_v(1,:,2:l-2) = 1.
  b_v(1,:,2:l-2) = 0.

  ! Update Boundary Conditions :: East // Symmetry Condition
  Aw_v(m-1,:,2:l-2) = 1.
  Ae_v(m-1,:,2:l-2) = 0.
  As_v(m-1,:,2:l-2) = 0.
  An_v(m-1,:,2:l-2) = 0.
  Ab_v(m-1,:,2:l-2) = 0.
  At_v(m-1,:,2:l-2) = 0.
  Ap_v(m-1,:,2:l-2) = 1.
  b_v(m-1,:,2:l-2) = 0.

  ! Update Boundary Conditions :: South // Symmetry Condition
  Aw_v(2:m-2,1,2:l-2) = 0.
  Ae_v(2:m-2,1,2:l-2) = 0.
  As_v(2:m-2,1,2:l-2) = 0.
  An_v(2:m-2,1,2:l-2) = 1.
  Ab_v(2:m-2,1,2:l-2) = 0.
  At_v(2:m-2,1,2:l-2) = 0.
  Ap_v(2:m-2,1,2:l-2) = 1.
  b_v(2:m-2,1,2:l-2) = 0.

  ! Update Boundary Conditions :: North // Symmetry Condition
  Aw_v(2:m-2,n,2:l-2) = 0.
  Ae_v(2:m-2,n,2:l-2) = 0.
  As_v(2:m-2,n,2:l-2) = 1.
  An_v(2:m-2,n,2:l-2) = 0.
  Ab_v(2:m-2,n,2:l-2) = 0.
  At_v(2:m-2,n,2:l-2) = 0.
  Ap_v(2:m-2,n,2:l-2) = 1.
  b_v(2:m-2,n,2:l-2) = 0.

  ! Update Boundary Conditions :: Bottom // Wall Condition Test with 0 Velocity
  Aw_v(:,:,1) = 0.
  Ae_v(:,:,1) = 0.
  As_v(:,:,1) = 0.
  An_v(:,:,1) = 0.
  Ab_v(:,:,1) = 0.
  At_v(:,:,1) = 0.
  Ap_v(:,:,1) = 1.
  b_v(:,:,1) = 0.

  ! Update Boundary Conditions :: Top // Wall Condition Test with 0 Velocity
  Aw_v(:,:,l-1) = 0.
  Ae_v(:,:,l-1) = 0.
  As_v(:,:,l-1) = 0.
  An_v(:,:,l-1) = 0.
  Ab_v(:,:,l-1) = 0.
  At_v(:,:,l-1) = 0.
  Ap_v(:,:,l-1) = 1.
  b_v(:,:,l-1) = 0.

  ! ====================================== !
  ! ====================================== !

  ! ====================================== !
  ! ====================================== !
  ! Calculate Coefficients for W-Velocity
  do i = 2, m-2
    do j = 2, n-2
	  do k = 2, l-1

	    ! Calculate diffusion terms
		Dw = dy*dz/dx
		De = dy*dz/dx
		Ds = dz*dx/dy
		Dn = dz*dx/dy
		Db = dx*dy/dz
		Dt = dx*dy/dz

		! Calculate adconvection terms
		Fw = dy*dz*(u(i,j,k-1)+u(i,j,k))/2.
		Fe = dy*dz*(u(i+1,j,k-1)+u(i+1,j,k))/2.
		Fs = dz*dx*(v(i,j,k-1)+v(i,j,k))/2.
		Fn = dz*dx*(v(i,j+1,k-1)+v(i,j+1,k))/2.
		Fb = dx*dy*(w(i,j,k-1)+w(i,j,k))/2.
		Ft = dx*dy*(w(i,j,k)+w(i,j,k+1))/2.

		! Calculate Coefficients :: Power Law Differencing Scheme
		Aw_w(i,j,k) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
		Ae_w(i,j,k) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
		As_w(i,j,k) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
		An_w(i,j,k) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)
		Ab_w(i,j,k) = Db*max(0.0,(1-0.1*abs(Fb/Db))**5)+max(Fb,0.0)
		At_w(i,j,k) = Dt*max(0.0,(1-0.1*abs(Ft/Dt))**5)+max(-Ft,0.0)

		! Update Ap_u Coefficient
		Ap_w(i,j,k) = Aw_w(i,j,k)+Ae_w(i,j,k)+As_w(i,j,k)+An_w(i,j,k)+Ab_w(i,j,k)+At_w(i,j,k)

		! Calculate Source Term
		b_w(i,j,k) = (Ra/Pr*((T(i,j,k)+T(i,j,k+1))/2.-0.5))*dx*dy*dz

	  end do
	end do
  end do

  ! Update Boundary Conditions :: West // Symmetry Condition
  Aw_w(1,:,2:l-1) = 0.
  Ae_w(1,:,2:l-1) = 1.
  As_w(1,:,2:l-1) = 0.
  An_w(1,:,2:l-1) = 0.
  Ab_w(1,:,2:l-1) = 0.
  At_w(1,:,2:l-1) = 0.
  Ap_w(1,:,2:l-1) = 1.
  b_w(1,:,2:l-1) = 0.

  ! Update Boundary Conditions :: East // Symmetry Condition
  Aw_w(m-1,:,2:l-1) = 1.
  Ae_w(m-1,:,2:l-1) = 0.
  As_w(m-1,:,2:l-1) = 0.
  An_w(m-1,:,2:l-1) = 0.
  Ab_w(m-1,:,2:l-1) = 0.
  At_w(m-1,:,2:l-1) = 0.
  Ap_w(m-1,:,2:l-1) = 1.
  b_w(m-1,:,2:l-1) = 0.

  ! Update Boundary Conditions :: South // Symmetry Condition
  Aw_w(2:m-2,1,2:l-1) = 0.
  Ae_w(2:m-2,1,2:l-1) = 0.
  As_w(2:m-2,1,2:l-1) = 0.
  An_w(2:m-2,1,2:l-1) = 1.
  Ab_w(2:m-2,1,2:l-1) = 0.
  At_w(2:m-2,1,2:l-1) = 0.
  Ap_w(2:m-2,1,2:l-1) = 1.
  b_w(2:m-2,1,2:l-1) = 0.

  ! Update Boundary Conditions :: North // Symmetry Condition
  Aw_w(2:m-2,n-1,2:l-1) = 0.
  Ae_w(2:m-2,n-1,2:l-1) = 0.
  As_w(2:m-2,n-1,2:l-1) = 1.
  An_w(2:m-2,n-1,2:l-1) = 0.
  Ab_w(2:m-2,n-1,2:l-1) = 0.
  At_w(2:m-2,n-1,2:l-1) = 0.
  Ap_w(2:m-2,n-1,2:l-1) = 1.
  b_w(2:m-2,n-1,2:l-1) = 0.

  ! Update Boundary Conditions :: Bottom // Wall Condition Test with 0 Velocity
  Aw_w(:,:,1) = 0.
  Ae_w(:,:,1) = 0.
  As_w(:,:,1) = 0.
  An_w(:,:,1) = 0.
  Ab_w(:,:,1) = 0.
  At_w(:,:,1) = 0.
  Ap_w(:,:,1) = 1.
  b_w(:,:,1) = 0.

  ! Update Boundary Conditions :: Top // Wall Condition Test with 0 Velocity
  Aw_w(:,:,l) = 0.
  Ae_w(:,:,l) = 0.
  As_w(:,:,l) = 0.
  An_w(:,:,l) = 0.
  Ab_w(:,:,l) = 0.
  At_w(:,:,l) = 0.
  Ap_w(:,:,l) = 1.
  b_w(:,:,l) = 0.

  ! ====================================== !
  ! ====================================== !

  return
end subroutine velocity3d_source

subroutine velocity3d_pseudo

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! ====================================== !
  ! Calculate Interior Psudeo U-Velocities
  do i = 2, m-1
    do j = 2, n-2
	  do k = 2, l-2

	    u_hat(i,j,k) = (Aw_u(i,j,k)*u(i-1,j,k)+ &
                        Ae_u(i,j,k)*u(i+1,j,k)+ &
                        As_u(i,j,k)*u(i,j-1,k)+ &
                        An_u(i,j,k)*u(i,j+1,k)+ &
                        Ab_u(i,j,k)*u(i,j,k-1)+ &
                        At_u(i,j,k)*u(i,j,k+1)+ &
                        b_u(i,j,k))/Ap_u(i,j,k)

	  end do
	end do
  end do

  ! Update Boundary U-Velocity :: West // Symmetry Condition
  u_hat(1,:,2:l-2) = u_hat(2,:,2:l-2)

  ! Update Boundary U-Velocity :: East // Symmetry Condition
  u_hat(m,:,2:l-2) = u_hat(m-1,:,2:l-2)

  ! Update Boundary U-Velocity :: South // Symmetry Condition
  u_hat(2:m-1,1,2:l-2) = u_hat(2:m-1,2,2:l-2)

  ! Update Boundary U-Velocity :: North // Symmetry Condition
  u_hat(2:m-1,n-1,2:l-2) = u_hat(2:m-1,n-2,2:l-2)

  ! Update Boundary U-Velocity :: Bottom // Wall Condition - Test with 0 Velocity
  u_hat(:,:,1) = 0.

  ! Update Boundary U-Velocity :: Top // Wall Condition - Test with 0 Velocity
  u_hat(:,:,l-1) = 0.

  ! ====================================== !
  ! ====================================== !
  ! Calculate Interior Psudeo V-Velocities
  do i = 2, m-2
    do j = 2, n-1
	  do k = 2, l-2

	    v_hat(i,j,k) = (Aw_v(i,j,k)*v(i-1,j,k)+ &
                        Ae_v(i,j,k)*v(i+1,j,k)+ &
                        As_v(i,j,k)*v(i,j-1,k)+ &
                        An_v(i,j,k)*v(i,j+1,k)+ &
                        Ab_v(i,j,k)*v(i,j,k-1)+ &
                        At_v(i,j,k)*v(i,j,k+1)+ &
                        b_v(i,j,k))/Ap_v(i,j,k)

	  end do
	end do
  end do

  ! Update Boundary V-Velocity :: West // Symmetry Condition
  v_hat(1,:,2:l-2) = v_hat(2,:,2:l-2)

  ! Update Boundary V-Velocity :: East // Symmetry Condition
  v_hat(m-1,:,2:l-2) = v_hat(m-2,:,2:l-2)

  ! Update Boundary V-Velocity :: South // Symmetry Condition
  v_hat(2:m-2,1,2:l-2) = v_hat(2:m-2,2,2:l-2)

  ! Update Boundary V-Velocity :: North // Symmetry Condition
  v_hat(2:m-2,n,2:l-2) = v_hat(2:m-2,n-1,2:l-2)

  ! Update Boundary V-Velocity :: Bottom // Wall Condition - Test with 0 Velocity
  v_hat(:,:,1) = 0.

  ! Update Boundary V-Velocity :: Top // Wall Condition - Test with 0 Velocity
  v_hat(:,:,l-1) = 0.

  ! ====================================== !
  ! ====================================== !
  ! Calculate Interior Psudeo W-Velocities
  do i = 2, m-2
    do j = 2, n-2
	  do k = 2, l-1

	    w_hat(i,j,k) = (Aw_w(i,j,k)*w(i-1,j,k)+ &
                        Ae_w(i,j,k)*w(i+1,j,k)+ &
                        As_w(i,j,k)*w(i,j-1,k)+ &
                        An_w(i,j,k)*w(i,j+1,k)+ &
                        Ab_w(i,j,k)*w(i,j,k-1)+ &
                        At_w(i,j,k)*w(i,j,k+1)+ &
                        b_w(i,j,k))/Ap_w(i,j,k)

	  end do
	end do
  end do

  ! Update Boundary W-Velocity :: West // Symmetry Condition
  w_hat(1,:,2:l-1) = w_hat(2,:,2:l-1)

  ! Update Boundary W-Velocity :: East // Symmetry Condition
  w_hat(m-1,:,2:l-1) = w_hat(m-2,:,2:l-1)

  ! Update Boundary W-Velocity :: South // Symmetry Condition
  w_hat(2:m-2,1,2:l-1) = w_hat(2:m-2,2,2:l-1)

  ! Update Boundary W-Velocity :: North // Symmetry Condition
  w_hat(2:m-2,n-1,2:l-1) = w_hat(2:m-2,n-2,2:l-1)

  ! Update Boundary W-Velocity :: Bottom // Wall Condition
  w_hat(:,:,1) = 0.

  ! Update Boundary W-Velocity :: Top // Wall Condition
  w_hat(:,:,l) = 0.

  ! ====================================== !
  ! ====================================== !

  return
end subroutine velocity3d_pseudo

subroutine velocity3d_solve

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k, fault

  ! ====================================== !
  ! ====================================== !
  ! Update U-Velocity Source Term and Add Relaxation
  do i = 2, m-1
    do j = 2, n-2
	  do k = 2, l-2

	    Ap_u(i,j,k) = Ap_u(i,j,k)/alpha_v
		  b_u(i,j,k) = b_u(i,j,k) + dy*dz*(P(i-1,j,k)-P(i,j,k))*beta_1 + (1.0-alpha_v)*Ap_u(i,j,k)*u_hat(i,j,k)

	  end do
	end do
  end do

  ! Solve U-Velocity Equation
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab_u, As_u, Aw_u, Ap_u, Ae_u, An_u, At_u, b_u, u_star, m, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab_u, As_u, Aw_u, Ap_u, Ae_u, An_u, At_u, b_u, u_star, m, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 2) then
    fault = 0
    do i = 3, maxit
      if (fault .eq. 0) then
        call solver3d_gmres(Ab_u, As_u, Aw_u, Ap_u, Ae_u, An_u, At_u, b_u, u_star, m, n-1, l-1, solver_tol, maxit, fault)
      end if
    end do
  elseif (solver .eq. 3) then
    call solver3d_bicg(Ab_u, As_u, Aw_u, Ap_u, Ae_u, An_u, At_u, b_u, u_star, m, n-1, l-1, solver_tol, maxit)
  else
    call solver3d_tdma(Ab_u, As_u, Aw_u, Ap_u, Ae_u, An_u, At_u, b_u, u_star, m, n-1, l-1, solver_tol, maxit)
  end if
  ! ====================================== !
  ! ====================================== !

  ! ====================================== !
  ! ====================================== !
  ! Update V-Velocity Source Term and Add Relaxation
  do i = 2, m-2
    do j = 2, n-1
	  do k = 2, l-2

	    Ap_v(i,j,k) = Ap_v(i,j,k)/alpha_v
		  b_v(i,j,k) = b_v(i,j,k) + dz*dx*(P(i,j-1,k)-P(i,j,k))*beta_1 + (1.0-alpha_v)*Ap_v(i,j,k)*v_hat(i,j,k)

	  end do
	end do
  end do

  ! Solve V-Velocity Equation
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab_v, As_v, Aw_v, Ap_v, Ae_v, An_v, At_v, b_v, v_star, m-1, n, l-1, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab_v, As_v, Aw_v, Ap_v, Ae_v, An_v, At_v, b_v, v_star, m-1, n, l-1, solver_tol, maxit)
  elseif (solver .eq. 2) then
    fault = 0
    do i = 3, maxit
      if (fault .eq. 0) then
        call solver3d_gmres(Ab_v, As_v, Aw_v, Ap_v, Ae_v, An_v, At_v, b_v, v_star, m-1, n, l-1, solver_tol, maxit, fault)
      end if
    end do
  elseif (solver .eq. 3) then
    call solver3d_bicg(Ab_v, As_v, Aw_v, Ap_v, Ae_v, An_v, At_v, b_v, v_star, m-1, n, l-1, solver_tol, maxit)
  else
    call solver3d_tdma(Ab_v, As_v, Aw_v, Ap_v, Ae_v, An_v, At_v, b_v, v_star, m-1, n, l-1, solver_tol, maxit)
  end if
  ! ====================================== !
  ! ====================================== !

  ! ====================================== !
  ! ====================================== !
  ! Update V-Velocity Source Term and Add Relaxation
  do i = 2, m-2
    do j = 2, n-2
	  do k = 2, l-1

	    Ap_w(i,j,k) = Ap_w(i,j,k)/alpha_v
		  b_w(i,j,k) = b_w(i,j,k) + dx*dy*(P(i,j,k-1)-P(i,j,k))*beta_1 + (1.0-alpha_v)*Ap_w(i,j,k)*w_hat(i,j,k)

	  end do
	end do
  end do

  ! Solve W-Velocity Equation
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab_w, As_w, Aw_w, Ap_w, Ae_w, An_w, At_w, b_w, w_star, m-1, n-1, l, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab_w, As_w, Aw_w, Ap_w, Ae_w, An_w, At_w, b_w, w_star, m-1, n-1, l, solver_tol, maxit)
  elseif (solver .eq. 2) then
    fault = 0
    do i = 3, maxit
      if (fault .eq. 0) then
        call solver3d_gmres(Ab_w, As_w, Aw_w, Ap_w, Ae_w, An_w, At_w, b_w, w_star, m-1, n-1, l, solver_tol, maxit, fault)
      end if
    end do
  elseif (solver .eq. 3) then
    call solver3d_bicg(Ab_w, As_w, Aw_w, Ap_w, Ae_w, An_w, At_w, b_w, w_star, m-1, n-1, l, solver_tol, maxit)
  else
    call solver3d_tdma(Ab_w, As_w, Aw_w, Ap_w, Ae_w, An_w, At_w, b_w, w_star, m-1, n-1, l, solver_tol, maxit)
  end if
  ! ====================================== !
  ! ====================================== !

  return
end subroutine velocity3d_solve

subroutine velocity3d_correct

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! ====================================== !
  ! ====================================== !
  ! Correct U-Velocity
  do i = 2, m-1
    do j = 2, n-2
	  do k = 2, l-2

	    u(i,j,k) = u_star(i,j,k) + dy*dz/Ap_u(i,j,k)*(P_prime(i-1,j,k)-P_prime(i,j,k))

	  end do
	end do
  end do

  ! Update Boundary U-Velocity :: West // Symmetry Condition
  u(1,:,2:l-2) = u(2,:,2:l-2)

  ! Update Boundary U-Velocity :: East // Symmetry Condition
  u(m,:,2:l-2) = u(m-1,:,2:l-2)

  ! Update Boundary U-Velocity :: South // Symmetry Condition
  u(2:m-1,1,2:l-2) = u(2:m-1,2,2:l-2)

  ! Update Boundary U-Velocity :: North // Symmetry Condition
  u(2:m-1,n-1,2:l-2) = u(2:m-1,n-2,2:l-2)

  ! Update Boundary U-Velocity :: Bottom // Wall Condition
  u(:,:,1) = 0.

  ! Update Boundary U-Velocity :: Top // Wall Condition
  u(:,:,l-1) = 0.
  ! ====================================== !
  ! ====================================== !

  ! ====================================== !
  ! ====================================== !
  ! Correct V-Velocity
  do i = 2, m-2
    do j = 2, n-1
	  do k = 2, l-2

	    v(i,j,k) = v_star(i,j,k) + dz*dx/Ap_v(i,j,k)*(P_prime(i,j-1,k)-P_prime(i,j,k))

	  end do
	end do
  end do

  ! Update Boundary V-Velocity :: West // Symmetry Condition
  v(1,:,2:l-2) = v(2,:,2:l-2)

  ! Update Boundary V-Velocity :: East // Symmetry Condition
  v(m-1,:,2:l-2) = v(m-2,:,2:l-2)

  ! Update Boundary V-Velocity :: South // Symmetry Condition
  v(2:m-2,1,2:l-2) = v(2:m-2,2,2:l-2)

  ! Update Boundary V-Velocity :: North // Symmetry Condition
  v(2:m-2,n,2:l-2) = v(2:m-2,n-1,2:l-2)

  ! Update Boundary V-Velocity :: Bottom // Wall Condition
  v(:,:,1) = 0.

  ! Update Boundary V-Velocity :: Top // Wall Condition
  v(:,:,l-1) = 0.
  ! ====================================== !
  ! ====================================== !

  ! ====================================== !
  ! ====================================== !
  ! Correct W-Velocity
  do i = 2, m-2
    do j = 2, n-2
	  do k = 2, l-1

	    w(i,j,k) = w_star(i,j,k) + dx*dy/Ap_w(i,j,k)*(P_prime(i,j,k-1)-P_prime(i,j,k))

	  end do
	end do
  end do

  ! Update Boundary V-Velocity :: West // Symmetry Condition
  w(1,:,2:l-1) = w(2,:,2:l-1)

  ! Update Boundary V-Velocity :: East // Symmetry Condition
  w(m-1,:,2:l-1) = w(m-2,:,2:l-1)

  ! Update Boundary V-Velocity :: South // Symmetry Condition
  w(2:m-2,1,2:l-1) = w(2:m-2,2,2:l-1)

  ! Update Boundary V-Velocity :: North // Symmetry Condition
  w(2:m-2,n-1,2:l-2) = w(2:m-2,n-2,2:l-2)

  ! Update Boundary V-Velocity :: Bottom // Wall Condition
  w(:,:,1) = 0.

  ! Update Boundary V-Velocity :: Top // Wall Condition
  w(:,:,l) = 0.
  ! ====================================== !
  ! ====================================== !

  return

end subroutine velocity3d_correct
