!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.

module mod_gravity

use mod_precision,only:rk => realkind

implicit none

private

public :: gravityClairaut
public :: gravitational_constant

real(kind=rk),parameter :: gravitational_constant = 6.67408e-11_rk

!===============================================================================
contains



!-------------------------------------------------------------------------------
pure elemental real(kind=rk)&
  function gravityClairaut(latitude) result(grav)
  !! Returns the gravitational acceleration at the Earth's surface as function
  !! of latitude, based on Clairaut's formula.
  !!
  !! TODO Reference
  real(kind=rk),intent(in) :: latitude !! Latitude [rad]
  grav = 9.780327_rk*(1.0026454_rk                &
                     -0.0026512_rk*cos(2*latitude)&
                     +0.0000058_rk*cos(2*latitude)**2)
endfunction gravityClairaut
!-------------------------------------------------------------------------------
endmodule mod_gravity
