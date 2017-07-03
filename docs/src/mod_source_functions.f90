!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.

module mod_source_functions

use mod_precision,only:ik => intkind,rk => realkind

implicit none

private

public :: sin_DCCM2012
public :: sds_DCCM2012
public :: sdt_DCCM2012
public :: sbf_DCCM2012
public :: sbf_JONSWAP
public :: snl_DCCM2012

contains


!-------------------------------------------------------------------------------
pure elemental function sin_DCCM2012(spectrum,wspd,wdir,input_height,ustar,&
  vonkarman) result(tendency)

  !! TODO implement currents averaged over the effective depth layer for
  !! modulation of phase speed.

  !! Returns a spectrum instance with the wave growth ($S_{in}$) tendency
  !! formulated by Donelan et al. (2012) and based on the sheltering hypothesis
  !! by Jeffreys (1924, 1925).
  !!
  !! The result instance has the units of 1/s. This source function must be
  !! re-evaluated if any of the input parameters change.
  !!
  !! References:
  !!
  !! Donelan, M. A., M. Curcic, S. S. Chen, and A. K. Magnusson, 2012: Modeling
  !! waves and wind stress, *J. Geophys. Res. Oceans*, **117**, C00J23,
  !! doi:10.1029/2011JC007787.
  !!
  !! Jeffreys, H., 1924: On the formation of waves by wind, *Proc. R. Soc. A*,
  !! **107**, 189–206.
  !!
  !! Jeffreys, H., 1925: On the formation of waves by wind, II, *Proc. R. Soc.
  !! A*, **110**, 341–347.

  use mod_spectrum,only:spectrum_type
  use mod_aerodynamic_drag,only:windAtReferenceHeight
  use mod_const,only:twopi

  type(spectrum_type),intent(in) :: spectrum
    !! Input spectrum instance
  real(kind=rk),intent(in) :: wspd
    !! Input wind speed [m/s]
  real(kind=rk),intent(in) :: wdir
    !! Input wind direction [rad], mathematical convention
  real(kind=rk),intent(in) :: input_height
    !! Height of input wind speed [m/s]
  real(kind=rk),intent(in) :: ustar
    !! Air-side friction velocity [m/s]
  real(kind=rk),intent(in) :: vonkarman
    !! Von Karman constant

  type(spectrum_type) :: tendency

  real(kind=rk),dimension(:,:),allocatable :: s_in

  real(kind=rk),dimension(:),allocatable :: f
  real(kind=rk),dimension(:),allocatable :: th
  real(kind=rk),dimension(:),allocatable :: k
  real(kind=rk),dimension(:),allocatable :: cp
  real(kind=rk),dimension(:),allocatable :: omega
  real(kind=rk),dimension(:),allocatable :: half_wavelength
  real(kind=rk),dimension(:),allocatable :: wspd_input

  real(kind=rk) :: grav
  real(kind=rk) :: rho_air
  real(kind=rk) :: rho_water

  real(kind=rk),parameter :: a1_windsea = 0.11_rk
  real(kind=rk),parameter :: a1_swell = 0.01_rk
  real(kind=rk),parameter :: a1_opposed = 0.10_rk

  real(kind=rk),parameter :: field_scale_negative = a1_opposed / a1_windsea
  real(kind=rk),parameter :: field_scale_swell = a1_swell / a1_opposed

  real(kind=rk),dimension(:,:),allocatable :: sheltering_coefficient

  integer :: nfreq,nfreqs
  integer :: ndir,ndirs

  tendency = spectrum

  grav = spectrum % getGravity()
  rho_air = spectrum % getAirDensity()
  rho_water = spectrum % getWaterDensity()
  f = spectrum % getFrequency()
  th = spectrum % getDirections()
  k = spectrum % getWavenumber()
  cp = spectrum % getPhaseSpeed()
  half_wavelength = 0.5_rk*spectrum % getWavelength()
  omega = twopi*f

  nfreqs = size(f)
  ndirs = size(th)

  ! Evaluate wind speed at height of half-wavelength of each wave component
  wspd_input = windAtReferenceHeight(wspd,input_height,half_wavelength,ustar,&
    vonkarman)

  allocate(s_in(nfreqs,ndirs))
  allocate(sheltering_coefficient(nfreqs,ndirs))

  ! Set the initial sheltering coefficient to a1_windsea everywhere
  sheltering_coefficient = a1_windsea
  do concurrent(nfreq = 1:nfreqs,ndir = 1:ndirs)
    ! If input is negative, adjust the sheltering coefficient to a1_opposed
    if(wspd_input(nfreq)*cos(wdir-th(ndir))-cp(nfreq) < 0)then
      sheltering_coefficient(nfreq,ndir) = sheltering_coefficient(nfreq,ndir)&
        * field_scale_negative
      ! If input is negative but has positive misalignment, adjust the
      ! sheltering coefficient to a1_swell
      if(cos(wdir-th(ndir)) > 0)then
        sheltering_coefficient(nfreq,ndir) = sheltering_coefficient(nfreq,ndir)&
        * field_scale_swell
      endif
    endif
  enddo

  do concurrent(ndir=1:ndirs)
    s_in(:,ndir) = sheltering_coefficient(:,ndir)*rho_air/rho_water            &
      *(wspd_input*cos(wdir-th(ndir))-cp)*abs(wspd_input*cos(wdir-th(ndir))-cp)&
      *omega*k/grav
  enddo

  tendency = s_in
  deallocate(s_in,sheltering_coefficient)

endfunction sin_DCCM2012
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function sds_DCCM2012(spectrum,sds_coefficient,sds_power,&
  mss_coefficient) result(tendency)

  !! Returns a spectrum instance with the wave dissipation ($S_{ds}$) tendency
  !! formulated by Donelan et al. (2012).
  !!
  !! The result instance has the units of 1/s. This source function must be
  !! re-evaluated every time the spectrum is updated.
  !!
  !! References:
  !!
  !! Donelan, M. A., B. K. Haus, W. J. Plant, and O. Troianowski, 2010:
  !! Modulation of short wind waves by long waves, *J. Geophys. Res. Oceans*,
  !! **115**, C10003, doi:10.1029/2009JC005794.
  !!
  !! Donelan, M. A., M. Curcic, S. S. Chen, and A. K. Magnusson, 2012: Modeling
  !! waves and wind stress, *J. Geophys. Res. Oceans*, **117**, C00J23,
  !! doi:10.1029/2011JC007787.

  use mod_spectrum,only:spectrum_type
  use mod_const,only:twopi

  type(spectrum_type),intent(in) :: spectrum
    !! Linear coefficient of the dissipation function
  real(kind=rk),intent(in) :: sds_coefficient
    !! Linear coefficient of the dissipation function
  real(kind=rk),intent(in) :: sds_power
    !! The exponent of the saturation spectrum
  real(kind=rk),intent(in) :: mss_coefficient
    !! Linear coefficient of the mean square slope adjustment to Sds
  type(spectrum_type) :: tendency
    !! Result tendency instance

  tendency = sds_coefficient * twopi * spectrum % getFrequency2d()&
    * spectrum % saturationSpectrum()**sds_power                  &
    * (1 + mss_coefficient * spectrum % meanSquareSlopeDirectional())**2

endfunction sds_DCCM2012
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function sdt_DCCM2012(spectrum,sdt_coefficient,ustar)&
  result(tendency)

  !! Returns a spectrum instance with the wave dissipation due to turbulence
  !! ($S_{dt}$) tendency formulated by Donelan et al. (2012).
  !!
  !! The result instance has the units of 1/s. This source function can be
  !! evaluated once and stored if wavenumber array, air and water densities, and
  !! friction velocity are held constant.
  !!
  !! References:
  !!
  !! Donelan, M. A., M. Curcic, S. S. Chen, and A. K. Magnusson, 2012: Modeling
  !! waves and wind stress, *J. Geophys. Res. Oceans*, **117**, C00J23,
  !! doi:10.1029/2011JC007787.

  use mod_spectrum,only:spectrum_type

  type(spectrum_type),intent(in) :: spectrum
    !! Spectrum instance
  real(kind=rk),intent(in) :: sdt_coefficient
    !! Linear coefficient of the turbulent dissipation function
  real(kind=rk),intent(in) :: ustar
    !! Air-side friction velocity [m/s]
  type(spectrum_type) :: tendency
    !! Result tendency instance

  tendency = sdt_coefficient * sqrt(spectrum % getAirDensity()&
    / spectrum % getWaterDensity()) * ustar * spectrum % getWavenumber2d()

endfunction sdt_DCCM2012
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function sbf_DCCM2012(spectrum,friction_coefficient,&
  percolation_coefficient) result(tendency)

  !! Returns a spectrum instance with the wave dissipation tendency due to
  !! bottom friction and percolation, formulated by Donelan et al. (2012).
  !!
  !! The result instance has the units of 1/s. If the mean water depth and
  !! wavenumber arrays are constant, this source function can be evaluated once
  !! and stored. Otherwise, the function must be re-evaluated if the mean water
  !! depth or wavenumber array change.
  !!
  !! References:
  !!
  !! Donelan, M. A., M. Curcic, S. S. Chen, and A. K. Magnusson, 2012: Modeling
  !! waves and wind stress, *J. Geophys. Res. Oceans*, **117**, C00J23,
  !! doi:10.1029/2011JC007787.
  !!
  !! Shemdin, P., K. Hasselmann, S. V. Hsiao, and K. Herterich, 1978: Non-linear
  !! and linear bottom interaction effects in shallow water, p347-372 in:
  !! Turbulent fluxes through the sea surface, wave dynamics and prediction, A.
  !! Favre and K. Hasselmann (eds), Plenum, New York, 677p.

  use mod_spectrum,only:spectrum_type

  type(spectrum_type),intent(in) :: spectrum
    !! Spectrum instance
  real(kind=rk),intent(in) :: friction_coefficient
    !! Bottom friction coefficient [m/s]
  real(kind=rk),intent(in) :: percolation_coefficient
    !! Bottom permeability coefficient [m/s]
  type(spectrum_type) :: tendency
    !! Result tendency instance

  real(kind=rk) :: d
  real(kind=rk),dimension(:,:),allocatable :: k

  d = spectrum % getDepth()
  k = spectrum % getWavenumber2d()

  tendency = friction_coefficient * k / sinh(2*k*d)&
           + percolation_coefficient / cosh(k*d)**2

endfunction sbf_DCCM2012
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function sbf_JONSWAP(spectrum,friction_coefficient)&
  result(tendency)

  !! Returns a spectrum instance with the bottom friction ($S_{bot}$) tendency
  !! based on JONSWAP field data (Hasselmann et al., 1973). It is also the
  !! default parameterization scheme used in the WAM model (WAMDIG, 1988).
  !!
  !! References:
  !!
  !! Hasselmann, K. et al., 1973. Measurements of wind-wave growth and swell
  !! decay during the Joint North Sea Wave Project (JONSWAP). *Dtsch. Hydrogh.
  !! Z.*, Suppl. A, **8**, 12, 95pp.
  !!
  !! WAMDI Group, 1988. The WAM model – a third generation ocean wave prediction
  !!  model. *J. Phys. Oceanogr.*, **18**, 1775–1810.

  use mod_spectrum,only:spectrum_type

  type(spectrum_type),intent(in) :: spectrum
    !! Spectrum instance
  real(kind=rk),intent(in) :: friction_coefficient
    !! Bottom friction coefficient [m/s]
  type(spectrum_type) :: tendency
    !! Result tendency instance

  tendency = 2 * friction_coefficient                                       &
    * (spectrum % getPhaseSpeed2d() / spectrum % getGroupSpeed2d() - 0.5_rk)&
    / (spectrum % getGravity()*spectrum % getDepth())

endfunction sbf_JONSWAP
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function snl_DCCM2012(spectrum,sds_tendency,snl_coefficient)&
  result(tendency)

  !! Returns a spectrum instance with the non-linear wave-wave energy transfer
  !! ($S_{nl}$) tendency formulated by Donelan et al. (2012).
  !!
  !! References:
  !!
  !! Donelan, M. A., M. Curcic, S. S. Chen, and A. K. Magnusson, 2012: Modeling
  !! waves and wind stress, *J. Geophys. Res. Oceans*, **117**, C00J23,
  !! doi:10.1029/2011JC007787.

  use mod_spectrum,only:spectrum_type
  use mod_const,only:twopi

  type(spectrum_type),intent(in) :: spectrum
    !! Spectrum instance
  type(spectrum_type),intent(in) :: sds_tendency
    !! Spectral dissipation tendency instance
  real(kind=rk),intent(in) :: snl_coefficient
    !! Linear coefficient of the dissipation function
  type(spectrum_type) :: tendency
    !! Result tendency instance

  real(kind=rk),dimension(:,:),allocatable :: sds_spectrum
  real(kind=rk),dimension(:,:),allocatable :: s_nl

  real(kind=rk),dimension(:),allocatable :: f
  real(kind=rk),dimension(:),allocatable :: th
  real(kind=rk),dimension(:),allocatable :: k
  real(kind=rk),dimension(:),allocatable :: dk

  real(kind=rk),dimension(:),allocatable :: w1,w2
  real(kind=rk) :: bf1,bf1a,bf2,dlnf

  integer :: nfreq,nfreqs
  integer :: ndir,ndirs

  f = spectrum % getFrequency()
  k = spectrum % getWavenumber()
  dk = spectrum % getWavenumberSpacing()
  th = spectrum % getDirections()

  nfreqs = size(f)
  ndirs = size(th)

  tendency = spectrum

  sds_spectrum = sds_tendency % getSpectrum()*spectrum % getSpectrum()

  allocate(s_nl(nfreqs,ndirs))
  s_nl = 0

  allocate(w1(nfreqs),w2(nfreqs))
  w1 = 0
  w2 = 0

  dlnf = (log(f(size(f)))-log(f(1)))/float(nfreqs-1)
  bf1 = exp(-16*dlnf**2)
  bf2 = exp(-64*dlnf**2)
  bf1a = bf1/(bf1+bf2)
  bf2  = bf2/(bf1+bf2)
  bf1  = bf1a

  do nfreq = 1,nfreqs-2
    w1(nfreq) = snl_coefficient*bf1*k(nfreq+1)*dk(nfreq+1)/(k(nfreq)*dk(nfreq))
    w2(nfreq) = snl_coefficient*bf2*k(nfreq+2)*dk(nfreq+2)/(k(nfreq)*dk(nfreq))
  enddo

  do concurrent(ndir=1:ndirs)
    do nfreq = 1,nfreqs-2
      s_nl(nfreq,ndir) = w1(nfreq)*sds_spectrum(nfreq+1,ndir)&
                       + w2(nfreq)*sds_spectrum(nfreq+2,ndir)&
                       - snl_coefficient*sds_spectrum(nfreq,ndir)
    enddo
  enddo

  tendency = s_nl
  deallocate(w1,w2,s_nl,sds_spectrum)

endfunction snl_DCCM2012
!-------------------------------------------------------------------------------
endmodule mod_source_functions
