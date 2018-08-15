! Solver file for 3d cfd Problems

subroutine solver3d_bicg(Ab, As, Aw, Ap, Ae, An, At, b, phi, m, n, l, tol, maxit)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, n, l
  integer, intent(in) :: maxit
  real(8), intent(in) :: tol
  real(8), dimension(m,n,l), intent(in) :: As, Aw, Ap, Ae, An, Ab, At, b
  real(8), dimension(m,n,l), intent(inout) :: phi

  ! Define internal variables
  integer :: i, j, k, itr
  real(8), dimension(m*n*l,7) :: A_values
  integer, dimension(7) :: A_distance
  real(8), dimension(m*n*l) :: r, rt, u, ut, c, Axx, Atut, x, b_values
  real(8) :: rho, rho1, gamma, beta, alpha, r_norm

  A_distance = (/-m*n, -m, -1, 0, 1, m, m*n/)

  ! Convert values into CDS Format
  do k = 1,l
    do j = 1,n
      do i = 1,m

        ! Compress stiffness matrix values
		    A_values(i+(j-1)*m+(k-1)*m*n,1) = -Ab(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,2) = -As(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,3) = -Aw(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,4) = Ap(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,5) = -Ae(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,6) = -An(i,j,k)
		    A_values(i+(j-1)*m+(k-1)*m*n,7) = -At(i,j,k)

        ! Compress right-hand side values
        b_values(i+(j-1)*m+(k-1)*m*n) = b(i,j,k)+1.0e-8

        ! Compress preconditioning values
        x(i+(j-1)*m+(k-1)*m*n) = phi(i,j,k)

	    end do
    end do
  end do

  ! ================================================================= !
  ! ====================== Start BiCG Algoritm ====================== !
  ! ================================================================= !

  ! Set x
  !if (sum(abs(x)) .eq. 0.0) then
  !  x = 1.
  !end if

  ! Compute r0
  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Axx)
  r = b_values - Axx

  ! Set rt
  rt = r

  ! Set rho1
  rho1 = 1.

  ! Set u and ut
  u = 0.
  ut = 0.

  ! Start BiCG Loop
  do itr = 1,maxit

    rho = ddot(m*n*l, r, 1, rt, 1)

	  beta = -rho/rho1

	  u = r - beta * u

	  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, u, c)

	  ut = rt - beta * ut

	  gamma = ddot(m*n*l, c, 1, rt, 1)

	  alpha = rho / gamma

	  x = x + alpha * u

	  r = r - alpha * c

	  r_norm = dnrm2(m*n*l, r, 1)

    !print *, 'Iteration: ', itr
    !print *, 'Relative residual: ', r_norm

	  if (r_norm .le. tol) then

      !print *, 'BiCG Algorithm successfully converged!'
      !print *, 'Number of Iterations: ', itr
      !print *, 'Relative residual: ', r_norm

      ! Update phi with the solution
      do k = 1,l
        do j = 1,n
          do i = 1,m
            phi(i,j,k) = x(i+(j-1)*m+(k-1)*m*n)
        end do
        end do
      end do

	    return

    elseif (itr .eq. maxit) then
      !print *, 'BiCG Algorithm did not converge!'
      !print *, 'Number of Iterations: ', itr
      !print *, 'Relative residual: ', r_norm
	  end if

	  call mkl_ddiagemv('T', m*n*l, A_values, m*n*l, A_distance, 7, ut, Atut)

	  rt = rt - alpha * Atut

	  rho1 = rho

  end do

  ! ================================================================= !
  ! ======================= End BiCG Algoritm ======================= !
  ! ================================================================= !

  return

end subroutine solver3d_bicg

subroutine solver3d_bicgstab(Ab, As, Aw, Ap, Ae, An, At, b, phi, m, n, l, tol, maxit)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, n, l
  integer, intent(in) :: maxit
  real(8), intent(in) :: tol
  real(8), dimension(m,n,l), intent(in) :: As, Aw, Ap, Ae, An, Ab, At, b
  real(8), dimension(m,n,l), intent(inout) :: phi

  ! Define internal variables
  integer :: i, j, k, itr
  real(8), dimension(m*n*l,7) :: A_values
  integer, dimension(7) :: A_distance
  real(8), dimension(m*n*l) :: r, r0, r1, x, p, Axp, Axs, s, Axx, b_values
  real(8) :: alpha, omega, beta, r_norm, s_norm


  A_distance = (/-m*n, -m, -1, 0, 1, m, m*n/)

  ! Convert values into CDS Format
  do k = 1,l
    do j = 1,n
      do i = 1,m

        ! Compress stiffness matrix values
		    A_values(i+(j-1)*m+(k-1)*m*n,1) = -Ab(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,2) = -As(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,3) = -Aw(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,4) = Ap(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,5) = -Ae(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,6) = -An(i,j,k)
		    A_values(i+(j-1)*m+(k-1)*m*n,7) = -At(i,j,k)

        ! Compress right-hand side values
        b_values(i+(j-1)*m+(k-1)*m*n) = b(i,j,k)+1.0e-8

        ! Compress preconditioning values
        x(i+(j-1)*m+(k-1)*m*n) = phi(i,j,k)

	  end do
    end do
  end do


  ! ================================================================= !
  ! ==================== Start BiCGStab Algoritm ==================== !
  ! ================================================================= !

  ! Set preconditioning array to 1
  !if (sum(abs(x)) .eq. 0.0) then
  !  x = 1.
  !end if

  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Axx)

  r = b_values - Axx

  r0 = 1

  if (m .eq. n) then
    r0(1) = 2
  end if

  p = r

  do itr = 1,maxit

    call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, p, Axp)

    alpha = ddot(m*n*l,r,1,r0,1)/ddot(m*n*l,Axp,1,r0,1)

    s = r - alpha*Axp

    call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, s, Axs)

    omega = ddot(m*n*l,Axs,1,s,1) / ddot(m*n*l,Axs,1,Axs,1)

    x = x + alpha*p + omega*s

    r1 = s - omega*Axs

    beta = ddot(m*n*l,r1,1,r0,1) / ddot(m*n*l,r,1,r0,1) * alpha/omega

    p = r1 + beta*(p-omega*Axp)

    r = r1

    ! Check convergence of BiCG
    !call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Axx)
	  !r_norm = abs(dnrm2(m*n*l, b_values-Axx, 1))
    r_norm = abs(dnrm2(m*n*l, r, 1))

    if (r_norm < tol) then
        print *, 'BiCGStab Algorithm successfully converged!'
        print *, 'Number of Iterations: ', itr
        print *, 'Relative residual: ', r_norm

        ! Update phi with the solution
        do k = 1,l
          do j = 1,n
            do i = 1,m
              phi(i,j,k) = x(i+(j-1)*m+(k-1)*m*n)
	        end do
          end do
        end do

        exit
    end if

    if (itr .eq. maxit) then
      print *, '************************************'
      print *, '************************************'
      print *, 'BiCGStab Algorithm did not converge!'
      print *, 'Number of Iterations: ', itr
      print *, 'Relative residual: ', r_norm
      print *, '************************************'
      print *, '************************************'
    else
      !print *, '************************************'
      !print *, '************************************'
      !print *, 'Number of Iterations: ', itr
      !print *, 'Relative residual: ', r_norm
      !print *, '************************************'
      !print *, '************************************'
    end if

  end do

  ! ================================================================= !
  ! ====================== End BiCG Algoritm ======================== !
  ! ================================================================= !

  ! Update phi with the solution
  do k = 1,l
    do j = 1,n
      do i = 1,m
        phi(i,j,k) = x(i+(j-1)*m+(k-1)*m*n)
	    end do
    end do
  end do

  return

end subroutine solver3d_bicgstab

subroutine solver3d_bicgstab2(Ab, As, Aw, Ap, Ae, An, At, b, phi, m, n, l, tol, maxit)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, n, l, maxit
  real(8), intent(in) :: tol
  real(8), dimension(m,n,l), intent(in) :: As, Aw, Ap, Ae, An, Ab, At, b
  real(8), dimension(m,n,l), intent(inout) :: phi

  ! Define internal variables
  integer :: i, j, k, itr
  integer, dimension(7) :: A_distance
  real(8) :: alpha, beta, gamma, mu, nu, rho, rho_1, tau, omega_1, omega_2, r_norm
  real(8), dimension(m*n*l) :: Ax, p, r, r0, r0_hat, s, t, w, v, x, x0, b_values
  real(8), dimension(m*n*l, 7) :: A_values

  !  Set A_distance
  A_distance = (/-m*n, -m, -1, 0, 1, m, m*n/)

  ! Convert values into CDS Format
  do k = 1,l
    do j = 1,n
      do i = 1,m

        ! Compress stiffness matrix values
		    A_values(i+(j-1)*m+(k-1)*m*n,1) = -Ab(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,2) = -As(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,3) = -Aw(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,4) = Ap(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,5) = -Ae(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,6) = -An(i,j,k)
		    A_values(i+(j-1)*m+(k-1)*m*n,7) = -At(i,j,k)

        ! Compress right-hand side values
        b_values(i+(j-1)*m+(k-1)*m*n) = b(i,j,k)+1.0e-8

        ! Compress preconditioning values
        x0(i+(j-1)*m+(k-1)*m*n) = phi(i,j,k)

	    end do
    end do
  end do

  ! ======================================================================== !
  ! ========== Start Bi-conjugate Gradients Stabilized (2) Method ========== !
  ! ======================================================================== !

  ! Check x0
  !if (sum(abs(x0)) .eq. 0.0) then
  !  x0 = 1.
  !end if

    ! Check inital guess
    call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x0, Ax)
    r0 = b_values - Ax

    r_norm = abs(dnrm2(m*n*l, r0, 1))

    if (r_norm < tol) then
      print *, 'Initial guess is a sufficient solution'
  	  print *, 'relative residual: ', r_norm
      return
    end if

    ! Set r0_hat
    r0_hat = r0

    ! Check that dot(r0, r0_hat) .ne. 0
    rho = ddot(m*n*l, r0, 1, r0_hat, 1)
    if (rho .eq. 0) then
      r0_hat = r0 + 1
    end if

    ! Set scalars
    rho = 1
    alpha = 1
    omega_1 = 1
    omega_2 = 1

    ! Set vectors
    w = 0
    v = 0
    p = 0

    ! Update r and x
    x = x0
    r = r0

    ! Start BiCGSTAB(2) Loop
    do itr = 1, maxit+1, 2

      rho_1 = -omega_2*rho

  	  ! Even Bi-CG step:
  	  rho = ddot(m*n*l,r,1,r0_hat,1)
  	  beta = alpha * rho / rho_1
  	  rho_1 = rho
  	  p = r - beta * (p - omega_1 * v - omega_2 * w)
  	  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, p, v)
  	  gamma = ddot(m*n*l, v, 1, r0_hat, 1)
  	  alpha = rho / gamma
  	  r = r - alpha * v
  	  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, r, s)
  	  x = x + alpha * p

  	  ! Check solution
  	  !call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Ax)
  	  !r_norm = abs(dnrm2(m*n*l, b_values - Ax, 1))
      r_norm = abs(dnrm2(m*n*l, r, 1))

      !print *, 'r_norm(0):', r_norm

  	  if (r_norm < tol) then
        !print *, 'BiCGSTAB(2) Algorithm successfully converged!(mid)'
        !print *, 'Number of Iterations: ', itr
        !print *, 'Relative residual: ', r_norm

        do k = 1,l
          do j = 1,n
            do i = 1,m
              phi(i,j,k) = x(i+(j-1)*m+(k-1)*m*n)
      	    end do
          end do
        end do

        return
      end if

  	  ! Odd Bi-CG step:
  	  rho = ddot(m*n*l, s, 1, r0_hat, 1)

  	  beta = alpha * rho / rho_1

  	  rho_1 = rho
  	  v = s - beta * v

  	  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, v, w)

  	  gamma = ddot(m*n*l, w, 1, r0_hat, 1)

  	  alpha = rho/gamma


  	  p = r - beta * p

  	  r = r - alpha * v

  	  s = s - alpha * w

  	  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, s, t)

  	  ! GMRES(2)-part
  	  omega_1 = ddot(m*n*l, r, 1, s, 1)

  	  mu = ddot(m*n*l, s, 1, s, 1)

  	  nu = ddot(m*n*l, s, 1, t, 1)

  	  tau = ddot(m*n*l, t, 1, t, 1)

  	  omega_2 = ddot(m*n*l, r, 1, t, 1)

  	  tau = tau - nu**2 / mu

  	  omega_2 = (omega_2 - nu * omega_1 / mu) / tau

  	  omega_1 = (omega_1 - nu*omega_2) / mu

  	  x = x + alpha * p + omega_1 * r + omega_2 * s

  	  r = r - omega_1 * s - omega_2 * t


  	  ! Check solution
  	  !call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Ax)
  	  !r_norm = abs(dnrm2(m*n*l, b_values - Ax, 1))
      r_norm = abs(dnrm2(m*n*l, r, 1))

  	  if (r_norm < tol) then
        !print *, 'BiCGSTAB(2) Algorithm successfully converged! (end)'
        !print *, 'Number of Iterations: ', itr+1
        !print *, 'Relative residual: ', r_norm

        do k = 1,l
          do j = 1,n
            do i = 1,m
              phi(i,j,k) = x(i+(j-1)*m+(k-1)*m*n)
      	    end do
          end do
        end do

        return
      end if

  	  if (itr .gt. maxit) then
        !print *, '************************************'
        !print *, '************************************'
        !print *, 'BiCGStab Algorithm did not converge!'
        !print *, 'Number of Iterations: ', itr
        !print *, 'Relative residual: ', r_norm
        !print *, '************************************'
        !print *, '************************************'
      else
        !print *, 'Number of Iterations: ', itr
        !print *, 'Relative residual: ', r_norm
      end if

    end do


  ! ======================================================================== !
  ! ============ End Bi-conjugate Gradients Stabilized (2) Method ========== !
  ! ======================================================================== !

  ! Update phi with the solution
  do k = 1,l
    do j = 1,n
      do i = 1,m
        phi(i,j,k) = x(i+(j-1)*m+(k-1)*m*n)
	    end do
    end do
  end do

  return

end subroutine solver3d_bicgstab2

subroutine solver3d_gmres(Ab, As, Aw, Ap, Ae, An, At, b, phi, m, n, l, tol, maxit, fault)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, n, l
  integer, intent(in) :: maxit
  integer, intent(inout) :: fault
  real(8), intent(in) :: tol
  real(8), dimension(m,n,l), intent(in) :: As, Aw, Ap, Ae, An, Ab, At, b
  real(8), dimension(m,n,l), intent(inout) :: phi

  ! Define internal variables :: Matrix Conversion
  integer :: i, j, k, ii, jj, kk, itr_init
  real(8), dimension(m*n*l,7) :: A_values
  integer, dimension(7) :: A_distance
  real(8), dimension(m*n*l) :: b_values

  ! Define internal variables :: GMRES

  real(8) :: b_norm, err, r_norm, mult1, temp_tri_sol
  real(8), dimension(m*n*l) :: r, Axx, x, e1, Qy, beta
  real(8), dimension(maxit) :: sn, cs, y, cs2, sn2
  real(8), dimension(m*n*l, maxit+1) :: Q
  real(8), dimension(maxit+1, maxit) :: H



  ! ================================================================= !
  ! ==================== Start Matrix Conversion ==================== !
  ! ================================================================= !

  ! Set fault
  fault = 0

  ! Set input integer
  itr_init = 0

  do i = maxit-1, 3, -1
    itr_init = itr_init + i
  end do

  ! Define distance between diagonals
  A_distance = (/-m*n, -m, -1, 0, 1, m, m*n/)

  mult1 = 1.0

  ! Convert values into CDS Format
  do k = 1,l
    do j = 1,n
      do i = 1,m

        ! Compress stiffness matrix values
		    A_values(i+(j-1)*m+(k-1)*m*n,1) = -Ab(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,2) = -As(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,3) = -Aw(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,4) = Ap(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,5) = -Ae(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,6) = -An(i,j,k)
		    A_values(i+(j-1)*m+(k-1)*m*n,7) = -At(i,j,k)

        ! Compress right-hand side values
        b_values(i+(j-1)*m+(k-1)*m*n) = b(i,j,k)

        ! Compress preconditioning values
        x(i+(j-1)*m+(k-1)*m*n) = phi(i,j,k)

	  end do
    end do
  end do

  ! ================================================================= !
  ! ==================== End Matrix Conversion ====================== !
  ! ================================================================= !

  ! ================================================================= !
  ! ====================== Start GMRES Algoritm ===================== !
  ! ================================================================= !

  ! Check x
  !if (sum(abs(x)) .eq. 0.0) then
  !  x = 1.
  !end if

  ! Compute residual vector
  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Axx)
  r = b_values - Axx

  ! Compute residual norm
  r_norm = dnrm2(m*n*l, r, 1)

  ! Compute norm of values vector
  b_norm = dnrm2(m*n*l, b_values, 1)

  ! Compute inital error
  err = r_norm / b_norm

  ! Initialize vectors
  sn = 0.
  cs = 0.
  e1 = 0.
  e1(1) = 1.

  ! Set Q to inital residual
  Q(:,1) = r / r_norm

  ! Set inital beta
  beta = r_norm * e1

  ! Start GMRES Loop
  do k = 1,maxit


    cs2 = cs
    sn2 = sn

    ! Run Arnoldi Loop
	  call arnoldi(A_values, A_distance, Q, H, m*n*l, k, maxit)

	  ! Eliminate the last element in H ith row and update rotation matrix
	  call givens_rotation(H, cs2, sn2, m*n*l, k, maxit)

    cs = cs2
    sn = sn2

	  ! Update the residual vector
	  beta(k+1) = -sn(k) * beta(k)
	  beta(k) = cs(k) * beta(k)
	  err = abs(beta(k+1)) / b_norm

	  if (err .le. tol) then

      y(1:k) = beta(1:k)

      call dtrsv('U', 'N', 'N', k, H(1:k,1:k), k, y(1:k), 1)

      Qy = 0.

      do i = 1,m*n*l
        do j = 1,k

          Qy(i) = Qy(i) + Q(i,j)*y(j)

        end do
      end do

      x = x + Qy

      do kk = 1,l
        do jj = 1,n
          do ii = 1,m
            phi(ii,jj,kk) = x(ii+(jj-1)*m+(kk-1)*m*n)
          end do
        end do
      end do

      !print *, 'GMRES Algorithm successfully converged!'
      !print *, 'Number of Iterations: ', k
      !print *, 'Relative residual: ', err
      fault = 1

      return

	  elseif (k .eq. maxit) then

      y(1:k) = beta(1:k)

      call dtrsv('U', 'N', 'N', k, H(1:k,1:k), k, y(1:k), 1)

      Qy = 0.

      do i = 1,m*n*l
        do j = 1,k

          Qy(i) = Qy(i) + Q(i,j)*y(j)

        end do
      end do

      x = x + Qy

      do kk = 1,l
        do jj = 1,n
          do ii = 1,m
            phi(ii,jj,kk) = x(ii+(jj-1)*m+(kk-1)*m*n)
          end do
        end do
      end do

      !print *, 'GMRES Algorithm did not converge!'
      !print *, 'Number of Iterations: ', k
      !print *, 'Relative residual: ', err
      fault = 0

      return

	  end if

  end do

  ! ================================================================= !
  ! ======================== End GMRES Algoritm ===================== !
  ! ================================================================= !



  return

end subroutine solver3d_gmres


subroutine arnoldi(A_values, A_distance, Q, H, m, k, maxit)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, k, maxit
  integer, dimension(7) :: A_distance
  real(8), dimension(m,7) :: A_values
  real(8), dimension(m,maxit+1) :: Q
  real(8), dimension(maxit+1,maxit) :: H

  ! Define internal variables
  integer :: i
  real(8), dimension(m) :: q_vec

  call mkl_ddiagemv('N', m, A_values, m, A_distance, 7, Q(:,k), q_vec)

  do i = 1,k

    H(i,k) =  ddot(m, q_vec, 1, Q(:,i), 1)
	  q_vec = q_vec - H(i,k) * Q(:,i)

  end do

  H(k+1,k) = dnrm2(m, q_vec, 1)
  q_vec = q_vec / H(k+1,k)
  Q(:,k+1) = q_vec

  return

end subroutine arnoldi

subroutine givens_rotation(H, cs, sn, m, k, maxit)

  ! Define implicit
  implicit none

  ! Define input variables
  integer, intent(in) :: m, k, maxit
  real(8), dimension(maxit) :: sn, cs
  real(8), dimension(maxit+1, maxit) :: H

  ! Define internal variables
  integer :: i
  real(8) :: temp

  do i = 1,k-1

    temp = cs(i) * H(i,k) + sn(i)*H(i+1,k)

	  H(i+1,k) = -sn(i)*H(i,k) + cs(i)*H(i+1,k)
	  H(i,k) = temp

  end do

  call calc_rotation(H, cs, sn, m, k, maxit)

  H(k,k) = cs(k)*H(k,k) + sn(k)*H(k+1,k)
  H(k+1,k) = 0.

  return

end subroutine givens_rotation


subroutine calc_rotation(H, cs, sn, m, k, maxit)

  ! Define implicit
  implicit none

  ! Define Input Variables
  integer, intent(in) :: m, k, maxit
  real(8), dimension(maxit) :: cs, sn
  real(8), dimension(maxit+1, maxit) :: H

  ! Define internal variables
  real(8) :: temp

  if (H(k,k) .eq. 0.) then

    cs(k) = 0.
	  sn(k) = 1.

  else

    temp = (H(k,k)**2.0+H(k+1,k)**2.0)**(0.5)
	  cs(k) = abs(H(k,k)) / temp
	  sn(k) = cs(k) * H(k+1,k) / H(k,k)

  end if

  return
end subroutine calc_rotation

subroutine solver3d_tdma(Ab, As, Aw, Ap, Ae, An, At, b, phi, m, n, l, tol, maxit)

  integer, intent(in) :: m, n, l, maxit
  real(8), dimension(m,n,l), intent(in) :: Aw, Ae, As, An, At, Ab, Ap, b
  real(8), dimension(m,n,l), intent(inout) :: phi

  integer :: i, j, k, itr
  real(8), dimension(m) :: awe, bwe, cwe, dwe, phiwe
  real(8), dimension(n) :: asn, bsn, csn, dsn, phisn
  real(8), dimension(l) :: abt, bbt, cbt, dbt, phibt
  real(8), dimension(m,n,l) :: r
  real(8) :: r_sum, tol

  do itr = 1,maxit

    ! ==================== West - East ==================== !
	do k = 1,l
	  do j = 1,n
	    do i = 1,m

		  awe(i) = Ap(i,j,k)
		  bwe(i) = -Ae(i,j,k)
		  cwe(i) = -Aw(i,j,k)

		  if (j .eq. 1) then

		    if (k .eq. 1) then
			  dwe(i) = b(i,j,k)+An(i,j,k)*phi(i,j+1,k)+At(i,j,k)*phi(i,j,k+1)
			elseif (k .eq. l) then
			  dwe(i) = b(i,j,k)+An(i,j,k)*phi(i,j+1,k)+Ab(i,j,k)*phi(i,j,k-1)
			else
			  dwe(i) = b(i,j,k)+An(i,j,k)*phi(i,j+1,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
			end if

		  elseif (j .eq. n) then

		    if (k .eq. 1) then
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+At(i,j,k)*phi(i,j,k+1)
			elseif (k .eq. l) then
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+Ab(i,j,k)*phi(i,j,k-1)
			else
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
			end if

		  else

		    if (k .eq. 1) then
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)+At(i,j,k)*phi(i,j,k+1)
			elseif (k .eq. l) then
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)+Ab(i,j,k)*phi(i,j,k-1)
			else
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
			end if

		  end if

		end do

    call solver1d_tdma(awe, bwe, cwe, dwe, phiwe, m)

    phi(:,j,k) = phiwe(:)

	  end do
	end do

	! =================== South - North =================== !
	do i = 1,m
	  do k = 1,l
	    do j = 1,n

		  asn(j) = Ap(i,j,k)
		  bsn(j) = -An(i,j,k)
		  csn(j) = -As(i,j,k)

		  if (k .eq. 1) then

		    if (i .eq. 1) then
		      dsn(j) = b(i,j,k)+Ae(i,j,k)*phi(i+1,j,k)+At(i,j,k)*phi(i,j,k+1)
		    elseif (i .eq. m) then
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+At(i,j,k)*phi(i,j,k+1)
		    else
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i+1,j,k)+At(i,j,k)*phi(i,j,k+1)
		    end if

		  elseif (k .eq. l) then

		    if (i .eq. 1) then
		      dsn(j) = b(i,j,k)+Ae(i,j,k)*phi(i+1,j,k)+Ab(i,j,k)*phi(i,j,k-1)
		    elseif (i .eq. m) then
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ab(i,j,k)*phi(i,j,k-1)
		    else
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i+1,j,k)+Ab(i,j,k)*phi(i,j,k-1)
		    end if

		  else

		    if (i .eq. 1) then
		      dsn(j) = b(i,j,k)+Ae(i,j,k)*phi(i+1,j,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
		    elseif (i .eq. m) then
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
		    else
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i+1,j,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
		    end if

		  end if

		end do

    call solver1d_tdma(asn, bsn, csn, dsn, phisn, n)

    phi(i,:,k) = phisn(:)

	  end do
	end do


	! ===================== Bottom - Top ================== !
	do j = 1,n
	  do i = 1,m
	    do k = 1,l

		  abt(k) = Ap(i,j,k)
		  bbt(k) = -At(i,j,k)
		  cbt(k) = -Ab(i,j,k)

		  if (j .eq. 1) then
  			if (i .eq. 1) then
  			  dbt(k) = b(i,j,k)+Ae(i,j,k)*phi(i+1,j,k)+An(i,j,k)*phi(i,j+1,k)
  			elseif (i .eq. m) then
  			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+An(i,j,k)*phi(i,j+1,k)
  			else
  			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i+1,j,k)+An(i,j,k)*phi(i,j+1,k)
  			end if
		  elseif (j .eq. n) then
  			if (i .eq. 1) then
  			  dbt(k) = b(i,j,k)+Ae(i,j,k)*phi(i+1,j,k)+As(i,j,k)*phi(i,j-1,k)
  			elseif (i .eq. m) then
  			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+As(i,j,k)*phi(i,j-1,k)
  			else
  			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i+1,j,k)+As(i,j,k)*phi(i,j-1,k)
  			end if
		  else
  			if (i .eq. 1) then
  			  dbt(k) = b(i,j,k)+Ae(i,j,k)*phi(i+1,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)
  			elseif (i .eq. m) then
  			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)
  			else
  			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i+1,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)
  			end if
		  end if

		end do

    call solver1d_tdma(abt, bbt, cbt, dbt, phibt, l)

    phi(i,j,:) = phibt(:)

	  end do
	end do

	! =================== Check Solution ================== !
  r = 0.

	do i = 2,m-1
	  do j = 2,n-1
	    do k = 2,l-1

			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))

        end do
      end do
    end do

	r_sum = 0.

	do i = 2,m-1
	  do j = 2,n-1
	    do k = 2,l-1
		    r_sum = r_sum + (r(i,j,k))**2.0
		  end do
	  end do
	end do

  r_sum = r_sum**(0.5)

	if (r_sum .le. tol) then
	  !print *, "TDMA Compelete."
    !print *, "r_sum:", r_sum
    !print *, "itrs:", itr
      return
  else
    !print *, "r_sum:", r_sum
    !print *, "itrs:", itr
	end if

  end do

  return

end subroutine solver3d_tdma

subroutine solver1d_tdma(a, b, c, d, phi, n)

  ! Define input / output variables
  integer, intent(in) :: n
  real(8), intent(in) :: a(n), b(n), c(n), d(n)
  real(8), intent(out) :: phi(n)

  ! Define internal variables
  real(8), dimension(n) :: P, Q

  integer :: i

  ! Start forward-substitution
  P(1) = b(1)/a(1)
  Q(1) = d(1)/a(1)

  do i = 2, n, 1
    P(i) = b(i)/(a(i)-c(i)*P(i-1))
    Q(i) = (d(i)-c(i)*Q(i-1))/(a(i)-c(i)*P(i-1))

  end do

  phi(n) = Q(n)

  ! Start backward-substitution
  do i = n-1, 1, -1

    phi(i) = Q(i)-P(i)*phi(i+1)

  end do

  return

end subroutine solver1d_tdma
