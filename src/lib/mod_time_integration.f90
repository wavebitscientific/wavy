!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD 3-clause license. See LICENSE for details.

module mod_time_integration

use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_domain,only:domain_type

implicit none

private

public :: backward_euler
public :: exact_exponential
public :: forward_euler
public :: integrate

interface integrate
  module procedure :: integrate_spectrum
  module procedure :: integrate_domain
endinterface integrate

interface exact_exponential
  module procedure :: exact_exponential_spectrum
  module procedure :: exact_exponential_domain
endinterface exact_exponential

interface backward_euler
  module procedure :: backward_euler_spectrum
  module procedure :: backward_euler_domain
endinterface backward_euler

interface forward_euler
  module procedure :: forward_euler_spectrum
  module procedure :: forward_euler_domain
endinterface forward_euler

contains

!-------------------------------------------------------------------------------
pure type(spectrum_type) function integrate_spectrum(func,initial,tendency,dt)
  !! Integrates spectrum forward in time using a time integration method
  !! provided as the argument `func`.
  interface  
    pure type(spectrum_type) function func(initial,tendency,dt)
      import :: spectrum_type,rk
      type(spectrum_type),intent(in) :: initial
      type(spectrum_type),intent(in) :: tendency
      real(kind=rk),intent(in) :: dt
    endfunction func
  endinterface  
  type(spectrum_type),intent(in) :: initial !! Initial spectrum instance
  type(spectrum_type),intent(in) :: tendency !! Spectrum tendency instance
  real(kind=rk),intent(in) :: dt !! Time step [s]
  integrate_spectrum = func(initial,tendency,dt)
endfunction integrate_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure type(domain_type) function integrate_domain(func,initial,tendency,dt)
  !! Integrates domain forward in time using a time integration method
  !! provided as the argument `func`.
  interface
    pure type(domain_type) function func(initial,tendency,dt)
      import :: domain_type,rk
      type(domain_type),intent(in) :: initial
      type(domain_type),intent(in) :: tendency
      real(kind=rk),intent(in) :: dt
    endfunction func
  endinterface
  type(domain_type),intent(in) :: initial !! Initial domain instance
  type(domain_type),intent(in) :: tendency !! Spectrum domain instance
  real(kind=rk),intent(in) :: dt !! Time step [s]
  integrate_domain = func(initial,tendency,dt)
endfunction integrate_domain
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function backward_euler_spectrum(initial,&
  tendency,dt) result(spec)
  !! Integrates a spectrum forward in time using a 1st order implicit backward 
  !! Euler integration scheme.
  type(spectrum_type),intent(in) :: initial !! Initial spectrum instance
  type(spectrum_type),intent(in) :: tendency !! Spectrum tendency instance
  real(kind=rk),intent(in) :: dt !! Time step [s]
  spec = initial / (- tendency * dt + 1._rk)
endfunction backward_euler_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(domain_type) function backward_euler_domain(initial,&
  tendency,dt) result(domain)
  !! Integrates a domain instance forward in time using a 1st order implicit 
  !! backward Euler integration scheme.
  type(domain_type),intent(in) :: initial !! Initial spectrum instance
  type(domain_type),intent(in) :: tendency !! Spectrum tendency instance
  real(kind=rk),intent(in) :: dt !! Time step [s]
  domain = initial
  call domain % setSpectrum(backward_euler(initial % getSpectrum(),&
    tendency % getSpectrum(),dt))
endfunction backward_euler_domain
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function exact_exponential_spectrum(initial,&
  tendency,dt) result(spec)
  !! Integrates a spectrum instance forward in time using the exact exponential.
  type(spectrum_type),intent(in) :: initial !! Initial spectrum instance
  type(spectrum_type),intent(in) :: tendency !! Spectrum tendency instance
  real(kind=rk),intent(in) :: dt !! Time step [s]
  spec = initial
  call spec % setSpectrum(exp(tendency % getSpectrum() * dt))
  !spec = (initial * exp(tendency % getSpectrum() * dt))
  spec = initial * spec
endfunction exact_exponential_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(domain_type) function exact_exponential_domain(initial,&
  tendency,dt) result(domain)
  !! Integrates a domain instance forward in time using the exact exponential.
  type(domain_type),intent(in) :: initial !! Initial spectrum instance
  type(domain_type),intent(in) :: tendency !! Spectrum tendency instance
  real(kind=rk),intent(in) :: dt !! Time step [s]
  domain = initial
  domain = exact_exponential(initial % getSpectrum(),&
    tendency % getSpectrum(),dt)
endfunction exact_exponential_domain
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function forward_euler_spectrum(initial,&
  tendency,dt) result(spec)
  !! Integrates a spectrum forward in time using a 1st order Euler integration
  !! scheme.
  type(spectrum_type),intent(in) :: initial !! Initial spectrum instance
  type(spectrum_type),intent(in) :: tendency !! Spectrum tendency instance
  real(kind=rk),intent(in) :: dt !! Time step [s]
  spec = initial + tendency * dt
endfunction forward_euler_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(domain_type) function forward_euler_domain(initial,&
  tendency,dt) result(domain)
  !! Integrates a domain instance forward in time using a 1st order Euler
  !! integration scheme.
  type(domain_type),intent(in) :: initial !! Initial spectrum instance
  type(domain_type),intent(in) :: tendency !! Spectrum tendency instance
  real(kind=rk),intent(in) :: dt !! Time step [s]
  domain = initial
  domain = forward_euler(initial % getSpectrum(),tendency % getSpectrum(),dt)
endfunction forward_euler_domain
!-------------------------------------------------------------------------------
endmodule mod_time_integration
