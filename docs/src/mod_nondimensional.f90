!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.

module mod_nondimensional

use mod_precision,only:realkind

implicit none

private

public :: nondimensionalDepth
public :: nondimensionalEnergy
public :: nondimensionalFetch
public :: nondimensionalFrequency
public :: nondimensionalRoughness_S1974
public :: nondimensionalRoughness_H1986
public :: nondimensionalTime
public :: waveAge

contains



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function waveAge(cp,u)
  !! Returns wave age, the ratio of phase speed and friction velocity or wind
  !! speed, depending on the caller's definition of wave age.
  real(kind=realkind),intent(in) :: cp
    !! Phase speed [m/s]
  real(kind=realkind),intent(in) :: u
    !! Friction velocity or wind speed [m/s]
  waveAge = cp / u
endfunction waveAge
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function nondimensionalDepth(wspd,depth,grav)
  !! Returns nondimensional depth based on input wind speed [m/s], mean water
  !! depth [m], and gravitational acceleration [m/s^2].
  real(kind=realkind),intent(in) :: wspd
    !! Wind speed at reference height [m/s]
  real(kind=realkind),intent(in) :: depth
    !! Mean water depth [m]
  real(kind=realkind),intent(in) :: grav
    !! Gravitational acceleration [m/s^2]
  nondimensionalDepth = grav*depth/wspd**2
endfunction nondimensionalDepth
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function nondimensionalEnergy(wspd,sigma,&
  grav)
  !! Returns nondimensional energy based on input wind speed, RMS of wave
  !! variance, and gravitational acceleration.
  real(kind=realkind),intent(in) :: wspd
    !! Wind speed at reference height [m/s]
  real(kind=realkind),intent(in) :: sigma
    !! Root mean square of wave variance
  real(kind=realkind),intent(in) :: grav
    !! Gravitational acceleration [m/s^2]
  nondimensionalEnergy = sigma**2*grav**2/wspd**4
endfunction nondimensionalEnergy
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function nondimensionalFetch(wspd,fetch,grav)
  !! Returns nondimensional energy based on input wind speed, RMS of wave
  !! variance, and gravitational acceleration.
  real(kind=realkind),intent(in) :: wspd
    !! Wind speed at reference height [m/s]
  real(kind=realkind),intent(in) :: fetch
    !! Fetch [m]
  real(kind=realkind),intent(in) :: grav
    !! Gravitational acceleration [m/s^2]
  nondimensionalFetch = grav*fetch/wspd**2
endfunction nondimensionalFetch
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function nondimensionalFrequency(wspd,fpeak,&
  grav)
  !! Returns nondimensional frequency based on input wind speed, peak frequency,
  !! and gravitational acceleration.
  real(kind=realkind),intent(in) :: wspd
    !! Wind speed at reference height [m/s]
  real(kind=realkind),intent(in) :: fpeak
    !! Peak frequency [Hz]
  real(kind=realkind),intent(in) :: grav
    !! Gravitational acceleration [m/s^2]
  nondimensionalFrequency = fpeak*wspd/grav
endfunction nondimensionalFrequency
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function nondimensionalRoughness_S1974(z0,&
  ustar,grav)
  !! Returns the aerodynamic roughness length scaled by friction velocity
  !! squared and gravitational acceleration, after Stewart (1974).
  !!
  !! References:
  !!
  !! TODO reference
  real(kind=realkind),intent(in) :: z0
    !! Roughness length [m]
  real(kind=realkind),intent(in) :: ustar
    !! Friction velocity [m/s]
  real(kind=realkind),intent(in) :: grav
    !! Gravitational acceleration [m/s^2]
  nondimensionalRoughness_S1974 = grav*z0/ustar**2
endfunction nondimensionalRoughness_S1974
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function nondimensionalRoughness_H1986(z0,hs)
  !! Returns the aerodynamic roughness length scaled by significant wave height,
  !! after Huang (1986).
  !!
  !! TODO reference
  real(kind=realkind),intent(in) :: z0
    !! Roughness length [m]
  real(kind=realkind),intent(in) :: hs
    !! Significant wave height [m]
  nondimensionalRoughness_H1986 = z0 / hs
endfunction nondimensionalRoughness_H1986
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=realkind) function nondimensionalTime(wspd,time,grav)
  !! Returns nondimensional time (duration) based on input wind speed, duration,
  !! and gravitational acceleration.
  real(kind=realkind),intent(in) :: wspd
    !! Wind speed at reference height [m/s]
  real(kind=realkind),intent(in) :: time
    !! Time [s]
  real(kind=realkind),intent(in) :: grav
    !! Gravitational acceleration [m/s^2]
  nondimensionalTime = grav*time/wspd
endfunction nondimensionalTime
!-------------------------------------------------------------------------------
endmodule mod_nondimensional
