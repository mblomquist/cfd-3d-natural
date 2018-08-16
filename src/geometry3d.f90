! geometry2d Subroutine for 2D CFD Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-18 (YYYY-MM-DD)
!
! This subroutine calculates the geometry properties and associated
! non-dimensional quantities for use in the SIMPLER method.

subroutine geometry3d

  ! Pull in standard variable header
  include "var3d.dec"

  ! Calculate non-dimensional quantities for SIMPLER
  dx = length/m/depth
  dy = width/n/depth
  dz = depth/l/depth

  return

end subroutine geometry3d
