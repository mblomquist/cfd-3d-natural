! velocity3d subroutines for 3d cfd problems

subroutine velocity3d_init
  ! velocity3d_init Subroutine for 3D CFD Problems
  !
  ! Written by Matt Blomquist
  ! Last Update: 2018-07-16 (YYYY-MM-DD)
  !
  ! This subroutine calculates the boundary source terms, sets interior source
  ! terms, sets boundary values, and initializes the velocity grids.
  !
  ! Boundary Conditions ::
  !   0 - Wall Type Condition
  !   1 - Fixed Value Type Condition
  !   2 - Symmetry Type Condition
  ! Pull in standard variable header

  include "var3d.dec"

  ! Initalize Velocity Fields
  u = 0.
  v = 0.
  w = 0.

  ! ====================================== !
  ! Initialize coefficients :: u
  Ab_u = 0.
  As_u = 0.
  Aw_u = 0.
  Ae_u = 0.
  An_u = 0.
  At_u = 0.
  Ap_u = 1.
  b_u = 0.

  ! Initalize source terms :: u
  Sp_u = 0.
  Su_u = 0.

  ! ====================================== !

  ! Initialize coefficients :: v
  Ab_v = 0.
  As_v = 0.
  Aw_v = 0.
  Ae_v = 0.
  An_v = 0.
  At_v = 0.
  Ap_v = 1.
  b_v = 0.

  ! Initalize source terms :: v
  Sp_v = 0.
  Su_v = 0.

  ! ====================================== !

  ! Initialize coefficients :: w
  Ab_w = 0.
  As_w = 0.
  Aw_w = 0.
  Ae_w = 0.
  An_w = 0.
  At_w = 0.
  Ap_w = 1.
  b_w = 0.

  ! Initalize source terms :: w
  Sp_w = 0.
  Su_w = 0.

  ! ====================================== !

  return

end subroutine velocity3d_init

subroutine velocity3d_boundary(direction)

  ! Boundary Conditions ::
  !   0 - Wall Type Condition
  !   1 - Fixed Value Type Condition
  !   2 - Symmetry Type Condition
  !
  ! Note:
  ! 1. For wall conditions the principle direction (i.e. west/east for u-velocity),
  !    the velocity value is set to zero. Non-principle walls are set to no slip.
  ! 2. For symmetry conditions, the source loop needs to be adjusted not to overwrite
  !    the assignment. This is the same for fixed value conditions.



  ! Include standard variable header
  include "var3d.dec"

  ! Define input variables
  character :: direction

  if (direction .eq. "u") then

    ! Initialize coefficients :: u
    Ab_u = 0.
    As_u = 0.
    Aw_u = 0.
    Ae_u = 0.
    An_u = 0.
    At_u = 0.
    Ap_u = 1.
    b_u = 0.

    ! Initalize source terms :: u
    Sp_u = 0.
    Su_u = 0.

    ! West
    b_u(1,2:n-2,2:l-2) = u0

    ! East
    b_u(m,2:n-2,2:l-2) = -u0

    ! South
    An_u(2:m-1,1,2:l-2) = 1.

    ! North
    As_u(2:m-1,n-1,2:l-2) = 1.

    ! Bottom
    b_u(:,:,1) = 0. !-2.0*mu*u0*dx*dy/dz

    ! Top
    b_u(:,:,l-1) = 0. !-2.0*mu*u0*dx*dy/dz

    ! West-South
    Ae_u(1,1,:) = 1.
    An_u(1,1,:) = 1.
    Ap_u(1,1,:) = 2.

    ! West-North
    Ae_u(1,n-1,:) = 1.
    As_u(1,n-1,:) = 1.
    Ap_u(1,n-1,:) = 2.

    ! East-South
    Aw_u(m,1,:) = 1.
    An_u(m,1,:) = 1.
    Ap_u(m,1,:) = 2.

    ! East-North
    Aw_u(m,n-1,:) = 1.
    As_u(m,n-1,:) = 1.
    Ap_u(m,n-1,:) = 2.

  end if

  if (direction .eq. "v") then

    ! Initialize coefficients :: u
    Ab_v = 0.
    As_v = 0.
    Aw_v = 0.
    Ae_v = 0.
    An_v = 0.
    At_v = 0.
    Ap_v = 1.
    b_v = 0.

    ! Initalize source terms :: u
    Sp_v = 0.
    Su_v = 0.

    ! West
    Ae_v(1,2:n-1,2:l-2) = 1.

    ! East
    Aw_v(m-1,2:n-1,2:l-2) = 1.

    ! South
    An_v(2:m-2,1,2:l-2) = 1.

    ! North
    As_v(2:m-2,n,2:l-2) = 1.

    ! Bottom
    b_v(:,:,1) = 0. !-2.0*mu*u0*dx*dy/dz

    ! Top
    b_v(:,:,l-1) = 0. !-2.0*mu*u0*dx*dy/dz

    ! West-South
    Ae_v(1,1,:) = 1.
    An_v(1,1,:) = 1.
    Ap_v(1,1,:) = 2.

    ! West-North
    Ae_v(1,n,:) = 1.
    As_v(1,n,:) = 1.
    Ap_v(1,n,:) = 2.

    ! East-South
    Aw_v(m-1,1,:) = 1.
    An_v(m-1,1,:) = 1.
    Ap_v(m-1,1,:) = 2.

    ! East-North
    Aw_v(m-1,n,:) = 1.
    As_v(m-1,n,:) = 1.
    Ap_v(m-1,n,:) = 2.

  end if

  if (direction .eq. "w") then

    ! Initialize coefficients :: u
    Ab_w = 0.
    As_w = 0.
    Aw_w = 0.
    Ae_w = 0.
    An_w = 0.
    At_w = 0.
    Ap_w = 1.
    b_w = 0.

    ! Initalize source terms :: u
    Sp_w = 0.
    Su_w = 0.

    ! West
    Ae_w(1,2:n-2,2:l-1) = 1.

    ! East
    Aw_w(m-1,2:n-2,2:l-1) = 1.

    ! South
    An_w(2:m-2,1,2:l-1) = 1.

    ! North
    As_w(2:m-2,n-1,2:l-1) = 1.

    ! Bottom
    b_w(:,:,1) = 0.

    ! Top
    b_w(:,:,l) = 0.

    ! West-South
    Ae_w(1,1,2:l-1) = 1.
    An_w(1,1,2:l-1) = 1.
    Ap_w(1,1,2:l-1) = 2.

    ! West-North
    Ae_w(1,n-1,2:l-1) = 1.
    As_w(1,n-1,2:l-1) = 1.
    Ap_w(1,n-1,2:l-1) = 2.

    ! East-South
    Aw_w(m-1,1,2:l-1) = 1.
    An_w(m-1,1,2:l-1) = 1.
    Ap_w(m-1,1,2:l-1) = 2.

    ! East-North
    Aw_w(m-1,n-1,2:l-1) = 1.
    As_w(m-1,n-1,2:l-1) = 1.
    Ap_w(m-1,n-1,2:l-1) = 2.

  end if

  return

end subroutine velocity3d_boundary

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

    ! Update coefficients at boundaries
    call velocity3d_boundary("u")

    ! Calculate interior coefficients
    do i = 2, m-1
      do j = 2, n-2
        do k = 2, l-2

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

          ! Check bottom / top
          if (k .eq. 1) then
            Ab_u(i,j,k) = 0.
          elseif (k .eq. l-1) then
            At_u(i,j,k) = 0.
          end if

		      ! Update Ap coefficient
		      Ap_u(i,j,k) = Aw_u(i,j,k)+Ae_u(i,j,k)+As_u(i,j,k)+An_u(i,j,k)+Ab_u(i,j,k)+At_u(i,j,k)-Sp_u(i,j,k)

		      ! Update b values
		      b_u(i,j,k) = beta_u*dx*dy*dz

	     end do
      end do
    end do

  end if

  ! ====================== V-Velocity ====================== !
  ! v-velocity update loop
  if (direction .eq. "v") then

    ! Update coefficients at boundaries
    call velocity3d_boundary("v")

    ! Calculate interior coefficients
    do i = 2,m-2
      do j = 2,n-1
        do k = 2,l-2

          ! Update convection terms
		      Fw = dy*dz*(u_star(i,j-1,k)+u_star(i,j,k))/2
		      Fe = dy*dz*(u_star(i+1,j-1,k)+u_star(i+1,j,k))/2
		      Fs = dz*dx*(v_star(i,j-1,k)+v_star(i,j,k))/2
		      Fn = dz*dx*(v_star(i,j,k)+v_star(i,j+1,k))/2
		      Fb = dx*dy*(w_star(i,j-1,k)+w_star(i,j,k))/2
		      Ft = dx*dy*(w_star(i,j-1,k+1)+w_star(i,j,k+1))/2

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

          ! Check bottom / top
          if (k .eq. 1) then
            Ab_v(i,j,k) = 0.
          elseif (k .eq. l-1) then
            At_v(i,j,k) = 0.
          end if

		      ! Update Ap coefficient
		      Ap_v(i,j,k) = Aw_v(i,j,k)+Ae_v(i,j,k)+As_v(i,j,k)+An_v(i,j,k)+Ab_v(i,j,k)+At_v(i,j,k)-Sp_v(i,j,k)+beta_v

          ! Update b values
          b_v(i,j,k) = 0.

        end do
      end do
    end do

  end if

  ! ====================== W-Velocity ====================== !
  ! w-velocity update loop
  if (direction .eq. "w") then

    ! Update coefficients at boundaries
    call velocity3d_boundary("w")

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
		      Ap_w(i,j,k) = Aw_w(i,j,k)+Ae_w(i,j,k)+As_w(i,j,k)+An_w(i,j,k)+Ab_w(i,j,k)+At_w(i,j,k)-Sp_w(i,j,k)

		      ! Update b values
		      b_w(i,j,k) = Pr*(((T(i,j,k)+T(i,j,k-1))/2.0)-0.5)*dx*dy*dz

        end do
      end do
    end do

  end if

  return

end subroutine velocity3d_source

subroutine velocity3d_solve

  implicit none

  ! Include variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k, fault

  ! ====================== U-Velocity ====================== !
  ! Update source terms :: u
  do i = 2,m-1
    do j = 2,n-2
      do k = 2,l-2

        Ap_u(i,j,k) = Ap_u(i,j,k)/alpha_v
        b_u(i,j,k) = b_u(i,j,k)+(P_star(i-1,j,k)-P_star(i,j,k))+(1.0-alpha_v)*Ap_u(i,j,k)*u_hat(i,j,k)

	     end do
    end do
  end do

  ! Solve u-velocity equation
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

  ! ====================== V-Velocity ====================== !
  ! Update source terms :: v
  do i = 2, m-2
    do j = 2, n-1
      do k = 2,l-2

        Ap_v(i,j,k) = Ap_v(i,j,k)/alpha_v
        b_v(i,j,k) = b_v(i,j,k)+(P_star(i,j-1,k)-P_star(i,j,k))+(1.0-alpha_v)*Ap_v(i,j,k)*v_hat(i,j,k)

      end do
    end do
  end do

  ! Solve v-velocity equation
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

  ! ====================== W-Velocity ====================== !
  ! Update source terms :: w
  do i = 2, m-2
    do j = 2, n-2
      do k = 2,l-1

        Ap_w(i,j,k) = Ap_w(i,j,k)/alpha_v
        b_w(i,j,k) = b_w(i,j,k)+(P_star(i,j,k-1)-P_star(i,j,k))+(1.0-alpha_v)*Ap_w(i,j,k)*w_hat(i,j,k)

      end do
    end do
  end do

  ! Solve w-velocity equation
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

  return

end subroutine velocity3d_solve

subroutine velocity3d_correct

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! ====================== U-Velocity ====================== !
  u = u_star
  v = v_star
  w = w_star

  ! Correct velocity values
  ! U-Velocity Update
  do i = 2,m-1
    do j = 2,n-2
      do k = 2,l-2
        u(i,j,k) = u_star(i,j,k)+dy*dz/Ap_u(i,j,k)*(P_prime(i-1,j,k)-P_prime(i,j,k))*alpha_v
      end do
    end do
  end do

  ! V-Velocity Update
  do i = 2,m-2
    do j = 2,n-1
      do k = 2,l-2
        v(i,j,k) = v_star(i,j,k)+dx*dz/Ap_v(i,j,k)*(P_prime(i,j-1,k)-P_prime(i,j,k))*alpha_v
      end do
    end do
  end do

  ! W-Velocity Update
  do i = 2,m-2
    do j = 2,n-2
      do k = 2,l-1
        w(i,j,k) = w_star(i,j,k)+dx*dy/Ap_w(i,j,k)*(P_prime(i,j,k-1)-P_prime(i,j,k))*alpha_v
      end do
    end do
  end do

  return

end subroutine velocity3d_correct

subroutine pseudo3d_solve

  implicit none

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! ========================== u_hat ========================== !
  u_hat = u_star
  v_hat = v_star
  w_hat = w_star

  do i = 2, m-1
    do j = 2, n-2
      do k = 2, l-2

          u_hat(i,j,k) = (Aw_u(i,j,k)*u_star(i-1,j,k)+ &
                        Ae_u(i,j,k)*u_star(i+1,j,k)+ &
                        As_u(i,j,k)*u_star(i,j-1,k)+ &
                        An_u(i,j,k)*u_star(i,j+1,k)+ &
                        Ab_u(i,j,k)*u_star(i,j,k-1)+ &
                        At_u(i,j,k)*u_star(i,j,k+1)+ &
                        b_u(i,j,k))/Ap_u(i,j,k)

      end do
    end do
  end do

  ! ========================== v_hat ========================== !

  do i = 2, m-2
    do j = 2, n-1
      do k = 2, l-2

          v_hat(i,j,k) = (Aw_v(i,j,k)*v_star(i-1,j,k)+ &
                          Ae_v(i,j,k)*v_star(i+1,j,k)+ &
                          As_v(i,j,k)*v_star(i,j-1,k)+ &
                          An_v(i,j,k)*v_star(i,j+1,k)+ &
                          Ab_v(i,j,k)*v_star(i,j,k-1)+ &
                          At_v(i,j,k)*v_star(i,j,k+1)+ &
                          b_v(i,j,k))/Ap_v(i,j,k)

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

  return

end subroutine pseudo3d_solve
