! temperature3d subroutines
!
! Written by Matt Blomquist
! Last Update: 2018-09-07 (YYYY-MM-DD)
!
! The subroutines below serve to calculate the coefficients and source terms
! for temperature in a 3D natural convection code. Additionally,
! the temperature subroutines is also included.
!
subroutine temperature3d_source

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k
  real(8) :: Dw, De, Ds, Dn, Db, Dt
  real(8) :: Fw, Fe, Fs, Fn, Fb, Ft

  ! Initalize Coefficients
  Aw_T = 0.
  Ae_T = 0.
  As_T = 0.
  An_T = 0.
  Ab_T = 0.
  At_T = 0.
  Ap_T = 1.
  b_T = 0.

  ! Calculate Interior Coefficients for Temperature
  do i = 2, m-2
    do j = 2, n-2
	  do k = 2, l-2

	    ! Update diffusion terms
      Dw = dy*dz/dx/Pr
      De = dy*dz/dx/Pr
      Ds = dz*dx/dy/Pr
      Dn = dz*dx/dy/Pr
      Db = dx*dy/dz/Pr
      Dt = dx*dy/dz/Pr

        ! Update advection terms
        Fw = dy*dz*u(i,j,k)
        Fe = dy*dz*u(i+1,j,k)
        Fs = dz*dx*v(i,j,k)
        Fn = dz*dx*v(i,j+1,k)
        Fb = dx*dy*w(i,j,k)
        Ft = dx*dy*w(i,j,k+1)

        ! Compute Coefficients - Power Law Differening Scheme
        Aw_T(i,j,k) = Dw*max(0.0,(1-0.1*abs(Fw/Dw))**5)+max(Fw,0.0)
        Ae_T(i,j,k) = De*max(0.0,(1-0.1*abs(Fe/De))**5)+max(-Fe,0.0)
        As_T(i,j,k) = Ds*max(0.0,(1-0.1*abs(Fs/Ds))**5)+max(Fs,0.0)
        An_T(i,j,k) = Dn*max(0.0,(1-0.1*abs(Fn/Dn))**5)+max(-Fn,0.0)
        Ab_T(i,j,k) = Db*max(0.0,(1-0.1*abs(Fb/Db))**5)+max(Fb,0.0)
        At_T(i,j,k) = Dt*max(0.0,(1-0.1*abs(Ft/Dt))**5)+max(-Ft,0.0)

  	    ! Update Ap coefficient
  	    Ap_T(i,j,k) = Aw_T(i,j,k)+Ae_T(i,j,k)+As_T(i,j,k)+An_T(i,j,k)+Ab_T(i,j,k)+At_T(i,j,k)

  	    ! Update b values
  	    b_T(i,j,k) = 0.

	  end do
	end do
  end do

  ! Update Boundary Conditions :: West // Symmetry Condition
  Aw_T(1,:,2:l-2) = 0.
  Ae_T(1,:,2:l-2) = 1.
  As_T(1,:,2:l-2) = 0.
  An_T(1,:,2:l-2) = 0.
  Ab_T(1,:,2:l-2) = 0.
  At_T(1,:,2:l-2) = 0.
  Ap_T(1,:,2:l-2) = 1.
  b_T(1,:,2:l-2) = 0.

  ! Update Boundary Conditions :: East // Symmetry Condition
  Aw_T(m-1,:,2:l-2) = 1.
  Ae_T(m-1,:,2:l-2) = 0.
  As_T(m-1,:,2:l-2) = 0.
  An_T(m-1,:,2:l-2) = 0.
  Ab_T(m-1,:,2:l-2) = 0.
  At_T(m-1,:,2:l-2) = 0.
  Ap_T(m-1,:,2:l-2) = 1.
  b_T(m-1,:,2:l-2) = 0.

  ! Update Boundary Conditions :: South // Symmetry Condition
  Aw_T(2:m-2,1,2:l-2) = 0.
  Ae_T(2:m-2,1,2:l-2) = 0.
  As_T(2:m-2,1,2:l-2) = 0.
  An_T(2:m-2,1,2:l-2) = 1.
  Ab_T(2:m-2,1,2:l-2) = 0.
  At_T(2:m-2,1,2:l-2) = 0.
  Ap_T(2:m-2,1,2:l-2) = 1.
  b_T(2:m-2,1,2:l-2) = 0.

  ! Update Boundary Conditions :: North // Symmetry Condition
  Aw_T(2:m-2,n-1,2:l-2) = 0.
  Ae_T(2:m-2,n-1,2:l-2) = 0.
  As_T(2:m-2,n-1,2:l-2) = 1.
  An_T(2:m-2,n-1,2:l-2) = 0.
  Ab_T(2:m-2,n-1,2:l-2) = 0.
  At_T(2:m-2,n-1,2:l-2) = 0.
  Ap_T(2:m-2,n-1,2:l-2) = 1.
  b_T(2:m-2,n-1,2:l-2) = 0.

  ! Update Boundary Conditions :: Bottom // Wall Condition - Hot
  Aw_T(:,:,1) = 0.
  Ae_T(:,:,1) = 0.
  As_T(:,:,1) = 0.
  An_T(:,:,1) = 0.
  Ab_T(:,:,1) = 0.
  At_T(:,:,1) = 0.
  Ap_T(:,:,1) = 1.
  b_T(:,:,1) = 1.

  ! Update Boundary Conditions :: Top // Wall Condition - Cold
  Aw_T(:,:,l-1) = 0.
  Ae_T(:,:,l-1) = 0.
  As_T(:,:,l-1) = 0.
  An_T(:,:,l-1) = 0.
  Ab_T(:,:,l-1) = 0.
  At_T(:,:,l-1) = 0.
  Ap_T(:,:,l-1) = 1.
  b_T(:,:,l-1) = 0.

  return

end subroutine temperature3d_source

subroutine temperature3d_solve(start)

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k, start, fault

  if (start .eq. 1) then

    ! Set Values for Boundaries
    Ap_Tr = Ap_T
    b_Tr = b_T

    ! Update Values with Relaxation
    do i = 2, m-2
      do j = 2, n-2
        do k = 2, l-2

          Ap_Tr(i,j,k) = Ap_T(i,j,k)/alpha_t
          b_Tr(i,j,k) = b_T(i,j,k) + (1.0-alpha_t)*Ap_Tr(i,j,k)*T(i,j,k)

        end do
      end do
    end do

  else

    ! Set Values
    Ap_Tr = Ap_T
    b_Tr = b_T

  end if

  ! Solve Temperature Equation
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab_T, As_T, Aw_T, Ap_Tr, Ae_T, An_T, At_T, b_Tr, T, m-1, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab_T, As_T, Aw_T, Ap_Tr, Ae_T, An_T, At_T, b_Tr, T, m-1, n-1, l-1, solver_tol, maxit)
  elseif (solver .eq. 2) then
    fault = 0
    do i = 3, maxit
      if (fault .eq. 0) then
        call solver3d_gmres(Ab_T, As_T, Aw_T, Ap_Tr, Ae_T, An_T, At_T, b_Tr, T, m-1, n-1, l-1, solver_tol, maxit, fault)
      end if
    end do
  elseif (solver .eq. 3) then
    call solver3d_bicg(Ab_T, As_T, Aw_T, Ap_Tr, Ae_T, An_T, At_T, b_Tr, T, m-1, n-1, l-1, solver_tol, maxit)
  else
    call solver3d_tdma(Ab_T, As_T, Aw_T, Ap_Tr, Ae_T, An_T, At_T, b_Tr, T, m-1, n-1, l-1, solver_tol, maxit)
  end if

  return
end subroutine temperature3d_solve
