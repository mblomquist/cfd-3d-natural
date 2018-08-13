! temperature3d Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This subroutine calculates the boundary source terms, sets interior source
! terms, sets boundary values, and initializes the temperature grid.
!
! Boundary Conditions ::
!   0 - Wall Type Condition
!   1 - Fixed Value Type Condition
!   2 - Symmetry Type Condition

subroutine temperature3d_init

  ! Pull in standard variable header
  include "var3d.dec"

  ! Initialize Temperature Field
  T = 0.

  ! Initialize coefficients
  Aw_T = 0.
  Ae_T = 0.
  As_T = 0.
  An_T = 0.
  Ab_T = 0.
  At_T = 0.
  Ap_T = 1.
  b_T = 0.

  ! Set temperature source terms
  Su_T = 0.
  Sp_T = 0.

  ! Define source update bounds
  if (T_bc_wc .eq. 0) then
    istart_T = 1
  else
    istart_T = 2
  end if

  if (T_bc_ec .eq. 0) then
    iend_T = m-1
  else
    iend_T = m-2
  end if

  if (T_bc_sc .eq. 0) then
    jstart_T = 1
  else
    jstart_T = 2
  end if

  if (T_bc_nc .eq. 0) then
    jend_T = n-1
  else
    jend_T = n-2
  end if

  if (T_bc_bc .eq. 0) then
    kstart_T = 1
  else
    kstart_T = 2
  end if

  if (T_bc_tc .eq. 0) then
    kend_T = l-1
  else
    kend_T = l-2
  end if

  return

end subroutine temperature3d_init

subroutine temperature3d_boundary

  ! Pull in standard variable header
  include "var3d.dec"

  ! West
  if (T_bc_wc .eq. 1) then
    Ap_T(1,2:n-2,2:l-2) = 1.
    b_T(1,2:n-2,2:l-2) = (T_bc_wv-T_c)/(delta_T)
  elseif (T_bc_Tc .eq. 2) then
    Ap_T(1,2:n-2,2:l-2) = 1.
    Ae_T(1,2:n-2,2:l-2) = 1.
  else
    Su_T(1,2:n-2,2:l-2) = 0.
  end if

  ! East
  if (T_bc_ec .eq. 1) then
    Ap_T(m-1,2:n-2,2:l-2) = 1.
    b_T(m-1,2:n-2,2:l-2) = (T_bc_ev-T_c)/(delta_T)
  elseif (T_bc_ec .eq. 2) then
    Ap_T(m-1,2:n-2,2:l-2) = 1.
    Aw_T(m-1,2:n-2,2:l-2) = 1.
  else
    Su_T(m-1,2:n-2,2:l-2) = 0.
  end if

  ! South
  if (T_bc_sc .eq. 1) then
    Ap_T(2:m-2,1,2:l-2) = 1.
    b_T(2:m-2,1,2:l-2) = (T_bc_sv-T_c)/(delta_T)
  elseif (T_bc_sc .eq. 2) then
    Ap_T(2:m-2,1,2:l-2) = 1.
    An_T(2:m-2,1,2:l-2) = 1.
  else
    Su_T(2:m-2,1,2:l-2) = 0.
  end if

  ! North
  if (T_bc_nc .eq. 1) then
    Ap_T(2:m-2,n-1,2:l-2) = 1.
    b_T(2:m-2,n-1,2:l-2) = (T_bc_nv-T_c)/(delta_T)
  elseif (T_bc_nc .eq. 2) then
    Ap_T(2:m-2,n-1,2:l-2) = 1.
    As_T(2:m-2,n-1,2:l-2) = 1.
  else
    Su_T(2:m-2,n-1,2:l-2) = 0.
  end if

  ! Bottom
  if (T_bc_bc .eq. 1) then
    Ap_T(2:m-2,2:n-2,1) = 1.
    b_T(2:m-2,2:n-2,1) = (T_bc_bv-T_c)/(delta_T)
  elseif (T_bc_bc .eq. 2) then
    Ap_T(2:m-2,2:n-2,1) = 1.
    At_T(2:m-2,2:n-2,1) = 1.
  else
    Su_T(2:m-2,2:n-2,1) = 0.
  end if

  ! Top
  if (T_bc_tc .eq. 1) then
    Ap_T(2:m-2,2:n-2,l-1) = 1.
    b_T(2:m-2,2:n-2,l-1) = (T_bc_tv-T_c)/(delta_T)
  elseif (T_bc_tc .eq. 2) then
    Ap_T(2:m-2,2:n-2,l-1) = 1.
    Ab_T(2:m-2,2:n-2,l-1) = 1.
  else
    Su_T(2:m-2,2:n-2,l-1) = 0.
  end if

  ! West-South
  Aw_T(1,1,2:l-2) = 0.
  Ae_T(1,1,2:l-2) = 1.
  As_T(1,1,2:l-2) = 0.
  An_T(1,1,2:l-2) = 1.
  Ab_T(1,1,2:l-2) = 0.
  At_T(1,1,2:l-2) = 0.

  Ap_T(1,1,2:l-2) = 2.
  b_T(1,1,2:l-2) = 0.

  ! West-North
  Aw_T(1,1,2:l-2) = 0.
  Ae_T(1,1,2:l-2) = 1.
  As_T(1,1,2:l-2) = 1.
  An_T(1,1,2:l-2) = 0.
  Ab_T(1,1,2:l-2) = 0.
  At_T(1,1,2:l-2) = 0.

  Ap_T(1,1,2:l-2) = 2.
  b_T(1,1,2:l-2) = 0.

  ! West-Bottom
  Aw_T(1,2:n-2,1) = 0.
  Ae_T(1,2:n-2,1) = 1.
  As_T(1,2:n-2,1) = 0.
  An_T(1,2:n-2,1) = 0.
  Ab_T(1,2:n-2,1) = 0.
  At_T(1,2:n-2,1) = 1.

  Ap_T(1,2:n-2,1) = 2.
  b_T(1,2:n-2,1) = 0.

  ! West-Top
  Aw_T(1,2:n-2,l-1) = 0.
  Ae_T(1,2:n-2,l-1) = 1.
  As_T(1,2:n-2,l-1) = 0.
  An_T(1,2:n-2,l-1) = 0.
  Ab_T(1,2:n-2,l-1) = 1.
  At_T(1,2:n-2,l-1) = 0.

  Ap_T(1,2:n-2,l-1) = 2.
  b_T(1,2:n-2,l-1) = 0.

  ! East-South
  Aw_T(m-1,1,2:l-2) = 1.
  Ae_T(m-1,1,2:l-2) = 0.
  As_T(m-1,1,2:l-2) = 0.
  An_T(m-1,1,2:l-2) = 1.
  Ab_T(m-1,1,2:l-2) = 0.
  At_T(m-1,1,2:l-2) = 0.

  Ap_T(m-1,1,2:l-2) = 2.
  b_T(m-1,1,2:l-2) = 0.

  ! East-North
  Aw_T(m-1,n-1,2:l-2) = 1.
  Ae_T(m-1,n-1,2:l-2) = 0.
  As_T(m-1,n-1,2:l-2) = 1.
  An_T(m-1,n-1,2:l-2) = 0.
  Ab_T(m-1,n-1,2:l-2) = 0.
  At_T(m-1,n-1,2:l-2) = 0.

  Ap_T(m-1,n-1,2:l-2) = 2.
  b_T(m-1,n-1,2:l-2) = 0.

  ! East Bottom
  Aw_T(m-1,2:n-2,1) = 1.
  Ae_T(m-1,2:n-2,1) = 0.
  As_T(m-1,2:n-2,1) = 0.
  An_T(m-1,2:n-2,1) = 0.
  Ab_T(m-1,2:n-2,1) = 0.
  At_T(m-1,2:n-2,1) = 1.

  Ap_T(m-1,2:n-2,1) = 2.
  b_T(m-1,2:n-2,1) = 0.

  ! East Top
  Aw_T(m-1,2:n-2,l-1) = 1.
  Ae_T(m-1,2:n-2,l-1) = 0.
  As_T(m-1,2:n-2,l-1) = 0.
  An_T(m-1,2:n-2,l-1) = 0.
  Ab_T(m-1,2:n-2,l-1) = 1.
  At_T(m-1,2:n-2,l-1) = 0.

  Ap_T(m-1,2:n-2,l-1) = 2.
  b_T(m-1,2:n-2,l-1) = 0.

  ! South-Bottom
  Aw_T(2:m-2,1,1) = 0.
  Ae_T(2:m-2,1,1) = 0.
  As_T(2:m-2,1,1) = 0.
  An_T(2:m-2,1,1) = 1.
  Ab_T(2:m-2,1,1) = 0.
  At_T(2:m-2,1,1) = 1.

  Ap_T(2:m-2,1,1) = 2.
  b_T(2:m-2,1,1) = 0.

  ! South-Top
  Aw_T(2:m-2,1,l-1) = 0.
  Ae_T(2:m-2,1,l-1) = 0.
  As_T(2:m-2,1,l-1) = 0.
  An_T(2:m-2,1,l-1) = 1.
  Ab_T(2:m-2,1,l-1) = 1.
  At_T(2:m-2,1,l-1) = 0.

  Ap_T(2:m-2,1,l-1) = 2.
  b_T(2:m-2,1,l-1) = 0.

  ! North-Bottom
  Aw_T(2:m-2,n-1,1) = 0.
  Ae_T(2:m-2,n-1,1) = 0.
  As_T(2:m-2,n-1,1) = 1.
  An_T(2:m-2,n-1,1) = 0.
  Ab_T(2:m-2,n-1,1) = 0.
  At_T(2:m-2,n-1,1) = 1.

  Ap_T(2:m-2,n-1,1) = 2.
  b_T(2:m-2,n-1,1) = 0.

  ! North-Top
  Aw_T(2:m-2,n-1,l-1) = 0.
  Ae_T(2:m-2,n-1,l-1) = 0.
  As_T(2:m-2,n-1,l-1) = 1.
  An_T(2:m-2,n-1,l-1) = 0.
  Ab_T(2:m-2,n-1,l-1) = 1.
  At_T(2:m-2,n-1,l-1) = 0.

  Ap_T(2:m-2,n-1,l-1) = 2.
  b_T(2:m-2,n-1,l-1) = 0.

  ! West-South-Bottom
  Aw_T(1,1,1) = 0.
  Ae_T(1,1,1) = 1.
  As_T(1,1,1) = 0.
  An_T(1,1,1) = 1.
  Ab_T(1,1,1) = 0.
  At_T(1,1,1) = 1.

  Ap_T(1,1,1) = 3.
  b_T(1,1,1) = 0.

  ! West-North-Bottom
  Aw_T(1,n-1,1) = 0.
  Ae_T(1,n-1,1) = 1.
  As_T(1,n-1,1) = 1.
  An_T(1,n-1,1) = 0.
  Ab_T(1,n-1,1) = 0.
  At_T(1,n-1,1) = 1.

  Ap_T(1,n-1,1) = 3.
  b_T(1,n-1,1) = 0.

  ! West-South-Top
  Aw_T(1,1,l-1) = 0.
  Ae_T(1,1,l-1) = 1.
  As_T(1,1,l-1) = 0.
  An_T(1,1,l-1) = 1.
  Ab_T(1,1,l-1) = 1.
  At_T(1,1,l-1) = 0.

  Ap_T(1,1,l-1) = 3.
  b_T(1,1,l-1) = 0.

  ! West-North-Top
  Aw_T(1,n-1,l-1) = 0.
  Ae_T(1,n-1,l-1) = 1.
  As_T(1,n-1,l-1) = 1.
  An_T(1,n-1,l-1) = 0.
  Ab_T(1,n-1,l-1) = 1.
  At_T(1,n-1,l-1) = 0.

  Ap_T(1,n-1,l-1) = 3.
  b_T(1,n-1,l-1) = 0.

  ! East-South-Bottom
  Aw_T(m-1,1,1) = 1.
  Ae_T(m-1,1,1) = 0.
  As_T(m-1,1,1) = 0.
  An_T(m-1,1,1) = 1.
  Ab_T(m-1,1,1) = 0.
  At_T(m-1,1,1) = 1.

  Ap_T(m-1,1,1) = 3.
  b_T(m-1,1,1) = 0.

  ! East-North-Bottom
  Aw_T(m-1,n-1,1) = 1.
  Ae_T(m-1,n-1,1) = 0.
  As_T(m-1,n-1,1) = 1.
  An_T(m-1,n-1,1) = 0.
  Ab_T(m-1,n-1,1) = 0.
  At_T(m-1,n-1,1) = 1.

  Ap_T(m-1,n-1,1) = 3.
  b_T(m-1,n-1,1) = 0.

  ! East-South-Top
  Aw_T(m-1,1,l-1) = 1.
  Ae_T(m-1,1,l-1) = 0.
  As_T(m-1,1,l-1) = 0.
  An_T(m-1,1,l-1) = 1.
  Ab_T(m-1,1,l-1) = 1.
  At_T(m-1,1,l-1) = 0.

  Ap_T(m-1,1,l-1) = 3.
  b_T(m-1,1,l-1) = 0.

  ! East-North-Top
  Aw_T(m-1,n-1,l-1) = 1.
  Ae_T(m-1,n-1,l-1) = 0.
  As_T(m-1,n-1,l-1) = 1.
  An_T(m-1,n-1,l-1) = 0.
  Ab_T(m-1,n-1,l-1) = 1.
  At_T(m-1,n-1,l-1) = 0.

  Ap_T(m-1,n-1,l-1) = 3.
  b_T(m-1,n-1,l-1) = 0.

  return

end subroutine temperature3d_boundary

subroutine temperature3d_source

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k
  real(8) :: Fw, Fe, Fs, Fn, Fb, Ft, Dw, De, Dn, Ds, Db, Dt

  ! Initialize coefficients
  Ab_T = 0.
  As_T = 0.
  Aw_T = 0.
  Ap_T = 0.
  Ae_T = 0.
  An_T = 0.
  At_T = 0.
  b_T = 0.

  ! Update boundary conditions
  call temperature3d_boundary

  ! Solve for source coefficients
  do i = istart_T, iend_T
    do j = jstart_T, jend_T
      do k = kstart_T, kend_T

        ! Update convective terms
        Fw = dy*dz*u(i,j,k)
        Fe = dy*dz*u(i+1,j,k)
        Fs = dz*dx*v(i,j,k)
        Fn = dz*dx*v(i,j+1,k)
        Fb = dx*dy*w(i,j,k)
        Ft = dx*dy*w(i,j,k+1)

        ! Update diffusion terms
        Dw = dy*dz/dx/(Ra*Pr)**(0.5)
        De = dy*dz/dx/(Ra*Pr)**(0.5)
        Dn = dz*dx/dy/(Ra*Pr)**(0.5)
        Ds = dz*dx/dy/(Ra*Pr)**(0.5)
        Db = dx*dy/dz/(Ra*Pr)**(0.5)
        Dt = dx*dy/dz/(Ra*Pr)**(0.5)

        ! Compute Coefficients - Power Law Differening Scheme
        Aw_T(i,j,k) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
        Ae_T(i,j,k) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
        As_T(i,j,k) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
        An_T(i,j,k) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)
        Ab_T(i,j,k) = Db*max(0.0,(1-0.1*abs(Fb/Db))**5)+max(Fb,0.0)
        At_T(i,j,k) = Dt*max(0.0,(1-0.1*abs(Ft/Dt))**5)+max(-Ft,0.0)

  	    ! Update Ap coefficient
  	    Ap_T(i,j,k) = Aw_T(i,j,k)+Ae_T(i,j,k)+As_T(i,j,k)+An_T(i,j,k)+Ab_T(i,j,k)+At_T(i,j,k)

        if (Ap_T(i,j,k) .eq. 0.) then

          print *, "False Diffusion (temperature)@:", i, j, k
          Ap_T(i,j,k) = 1.0

        end if

  	    ! Update b values
  	    b_T(i,j,k) = 0.

	     end do
    end do
  end do

  return

end subroutine temperature3d_source

subroutine temperature3d_solve(start)

  ! Include variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k, fault, start

  ! Update source terms
  call temperature3d_source

  if (start .eq. 0) then
    alpha_temp = 1.0
  else
    alpha_temp = alpha_t
  end if

  ! Update source terms
  do i = istart_T, iend_T
    do j = jstart_T, jend_T
      do k = kstart_T, kend_T
        Ap_T(i,j,k) = Ap_T(i,j,k)/alpha_temp
        b_T(i,j,k) = Su_T(i,j,k)+(1.0-alpha_temp)*Ap_T(i,j,k)*T(i,j,k)
      end do
    end do
  end do

  ! Solve velocity Equations
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab_T, As_T, Aw_T, Ap_T, Ae_T, An_T, At_T, b_T, T, m-1, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab_T, As_T, Aw_T, Ap_T, Ae_T, An_T, At_T, b_T, T, m-1, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 2) then
    fault = 0
    do i = 3, maxit
      if (fault .eq. 0) then
        call solver3d_gmres(Ab_T, As_T, Aw_T, Ap_T, Ae_T, An_T, At_T, b_T, T, m-1, n-1, l-1, solver_tol, maxit, fault)
      end if
    end do
  elseif (solver .eq. 3) then
    call solver3d_bicg(Ab_T, As_T, Aw_T, Ap_T, Ae_T, An_T, At_T, b_T, T, m-1, n-1, l-1, solver_tol, maxit)
  else
    call solver3d_tdma(Ab_T, As_T, Aw_T, Ap_T, Ae_T, An_T, At_T, b_T, T, m-1, n-1, l-1, solver_tol, maxit)
  end if

  return

end subroutine temperature3d_solve
