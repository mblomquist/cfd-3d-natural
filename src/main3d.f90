! 3D CFD Solver for Natural Convection
!
! Written by Matt Blomquist
! Last Update: 2018-09-07 (YYYY-MM-DD)
!
! This program solves a three-dimensional, steady state
! Rayleigh-Benard, Natural Convection problem
!
program main3d

  ! Define implicit
  implicit none

  ! Pull in standard variable header
  include "var3d.dec"

  ! Initalize problem
  call initialize3d

  ! Output Problem Parameters to Terminal
  print *, "3D Natural Convection"
  print *, ""
  print *, "Grid Size (m, n, l):"
  print *, m, n, l
  print *, ""
  print *, "Dimensionless Quantities"
  print *, "Ra:", Ra
  print *, "Pr:", Pr
  print *, ""
  print *, "Characteristic Properties"
  print *, "Velocity scale:", u0
  print *, "Length scale", depth
  print *, "Pressure scale", depth/rho/u0**(2.0)
  print *, ""

  ! Start SIMPLER Algorithm
  print *, "Starting SIMPLER Algorithm..."

  ! Start clock
  call cpu_time(t0)

  call simpler3d

  ! Stop clock
  call cpu_time(tf)

  print *, "SIMPLER Algorithm duration: ", (tf-t0)/4.0

  ! Output results
  print *, "Writing results to file..."
  call output3d
  print *, "Program complete."

 ! End program
 end program main3d
