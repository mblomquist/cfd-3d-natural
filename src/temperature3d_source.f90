! temperature3d_source Subroutine for 3D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-06-27 (YYYY-MM-DD)
!
! This subourtine calculates the coefficients used in the solution for
! temperature in the SIMPLER method.
!

subroutine temperature3d_source

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k
  real(8) :: Fw, Fe, Fs, Fn, Fb, Ft, Dw, De, Dn, Ds, Db, Dt

  ! Solve for source coefficients
  do i = 2,m-2
    do j = 2,n-2
      do k = 2,l-2

        ! Update convective terms
        Fw = dy*dz*u(i,j,k)
        Fe = dy*dz*u(i+1,j,k)
        Fs = dz*dx*v(i,j,k)
        Fn = dz*dx*v(i,j+1,k)
        Fb = dx*dy*w(i,j,k)
        Ft = dx*dy*w(i,j,k+1)

        ! Update diffusion terms
        Dw = dy*dz/dx/(Ra/Pr)**(0.5)
        De = dy*dz/dx/(Ra/Pr)**(0.5)
        Ds = dz*dx/dy/(Ra/Pr)**(0.5)
        Dn = dz*dx/dy/(Ra/Pr)**(0.5)
        Db = dx*dy/dz/(Ra/Pr)**(0.5)
        Dt = dx*dy/dz/(Ra/Pr)**(0.5)

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
  	    b_T(i,j,k) = Su_T(i,j,k)*dx*dy*dz

	     end do
    end do
  end do

  call temperature3d_boundary

  return

end subroutine temperature3d_source
