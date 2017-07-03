!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.
!
!===============================================================================
module mod_linear_wave_theory

use mod_precision,only:intkind,realkind

implicit none

private

public :: elevation
public :: pressure
public :: horizontalAcceleration
public :: horizontalVelocity
public :: verticalAcceleration
public :: verticalVelocity

!===============================================================================
contains



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function elevation(x,t,a,k,omega)

  !! Returns the elevation [m] of a sinusoid wave given its amplitude [m],
  !! wavenumber [rad/m], and frequency [Hz].

  real(kind=realkind),intent(in) :: x
    !! Horizontal space [m]
  real(kind=realkind),intent(in) :: t
    !! Time [s]
  real(kind=realkind),intent(in) :: a
    !! Wave amplitude [m]
  real(kind=realkind),intent(in) :: k
    !! Wavenumber [rad/m]
  real(kind=realkind),intent(in) :: omega
    !! Angular frequency [rad]

  elevation = a*sin(k*x-omega*t)

endfunction elevation
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function pressure(x,z,t,a,k,omega,rho,grav)

  !! Returns the pressure [Pa] at depth z (negative downward) for a sinusoid
  !! wave given its amplitude [m], wavenumber [rad/m], and frequency [Hz].

  real(kind=realkind),intent(in) :: x
    !! Horizontal space [m]
  real(kind=realkind),intent(in) :: z
    !! Vertical displacement [m] from the surface, negative downward
  real(kind=realkind),intent(in) :: t
    !! Time [s]
  real(kind=realkind),intent(in) :: a
    !! Wave amplitude [m]
  real(kind=realkind),intent(in) :: k
    !! Wavenumber [rad/m]
  real(kind=realkind),intent(in) :: omega
    !! Angular frequency [rad]
  real(kind=realkind),intent(in) :: rho
    !! Water density [kg/m^3]
  real(kind=realkind),intent(in) :: grav
    !! Gravitational acceleration [m/s^2]

  pressure = -rho*grav*(elevation(x,t,a,k,omega)-z)

endfunction pressure
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function horizontalAcceleration(x,z,t,a,k,omega)

  !! Returns the horizontal acceleration of a water particle under a sinusoid wave,
  !! given its amplitude, wavenumber, and frequency.

  real(kind=realkind),intent(in) :: x
    !! Horizontal space [m]
  real(kind=realkind),intent(in) :: z
    !! Vertical space, negative downward [m]
  real(kind=realkind),intent(in) :: t
    !! Time [s]
  real(kind=realkind),intent(in) :: a
    !! Wave amplitude [m]
  real(kind=realkind),intent(in) :: k
    !! Wavenumber [rad/m]
  real(kind=realkind),intent(in) :: omega
    !! Angular frequency [rad]

  horizontalAcceleration = -a*omega**2*cos(k*x-omega*t)*exp(k*z)

endfunction horizontalAcceleration
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function horizontalVelocity(x,z,t,a,k,omega)

  !! Returns the horizontal velocity of a water particle under a sinusoid wave,
  !! given its amplitude, wavenumber, and frequency.

  real(kind=realkind),intent(in) :: x
    !! Horizontal space [m]
  real(kind=realkind),intent(in) :: z
    !! Vertical space, negative downward [m]
  real(kind=realkind),intent(in) :: t
    !! Time [s]
  real(kind=realkind),intent(in) :: a
    !! Wave amplitude [m]
  real(kind=realkind),intent(in) :: k
    !! Wavenumber [rad/m]
  real(kind=realkind),intent(in) :: omega
    !! Angular frequency [rad]

  horizontalVelocity = a*omega*sin(k*x-omega*t)*exp(k*z)

endfunction horizontalVelocity
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind)&
  function verticalAcceleration(x,z,t,a,k,omega)

  !! Returns the vertical acceleration of a water particle under a sinusoid wave,
  !! given its amplitude, wavenumber, and frequency.

  real(kind=realkind),intent(in) :: x
    !! Horizontal space [m]
  real(kind=realkind),intent(in) :: z
    !! Vertical space, negative downward [m]
  real(kind=realkind),intent(in) :: t
    !! Time [s]
  real(kind=realkind),intent(in) :: a
    !! Wave amplitude [m]
  real(kind=realkind),intent(in) :: k
    !! Wavenumber [rad/m]
  real(kind=realkind),intent(in) :: omega
    !! Angular frequency [rad]

  verticalAcceleration = -a*omega**2*sin(k*x-omega*t)*exp(k*z)

endfunction verticalAcceleration
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function verticalVelocity(x,z,t,a,k,omega)

  !! Returns the vertical velocity of a water particle under a sinusoid wave,
  !! given its amplitude, wavenumber, and frequency.

  real(kind=realkind),intent(in) :: x
    !! Horizontal space [m]
  real(kind=realkind),intent(in) :: z
    !! Vertical space, negative downward [m]
  real(kind=realkind),intent(in) :: t
    !! Time [s]
  real(kind=realkind),intent(in) :: a
    !! Wave amplitude [m]
  real(kind=realkind),intent(in) :: k
    !! Wavenumber [rad/m]
  real(kind=realkind),intent(in) :: omega
    !! Angular frequency [rad]

  verticalVelocity = -a*omega*cos(k*x-omega*t)*exp(k*z)

endfunction verticalVelocity
!-------------------------------------------------------------------------------
endmodule mod_linear_wave_theory
