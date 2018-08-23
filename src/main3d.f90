! 3D CFD Solver for Natural Convection
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This program solves a three-dimensional, steady state CFD problem
! using the SIMPLER method.

program main3d

  implicit none

  ! Include standard variable header
  !   Standard variable header establishes parameters and global values
  !   subroutines will call these modify these values.
  include "var3d.dec"

  ! Initialize Problem
  !   The initialize subroutine modifies global vaules to incorporate
  !   boundary conditions, calculate global values, etc.
  call initialize3d

  ! Output Problem Parameters to Terminal
  print *, 'Problem Initialization Complete.'
  print *, ''
  print *, 'Grid size: ', m, n, l
  print *, 'delta_T:', delta_T
  print *, 'length:', length
  print *, 'u0:', u0
  print *, 'Rayleigh Number: ', Ra
  print *, 'Prandtl Number: ', Pr
  print *, 'dx, dy, dz:', dx, dy, dz
  print *, ''

  call cpu_time(start_time)

  ! Start SIMPLER Algorithm
  !   The SIMPLER Algorithm solves the problem using the iterative method.
  print *, 'Starting SIMPLER Algorithm...'
  call simpler3d
  print *, 'SIMPLER Algorithm Complete.'

  call cpu_time(end_time)

  print *, "SIMPLER Algorithm Duration:", end_time-start_time

  ! Output results
  !   The output subroutine writes the final values of P, T, U, V, W,
  !   and convergance outputs.
  print *, 'Writing results to file.'
  call output3d_results
  print *, 'Program complete.'

! End program
end program main3d
