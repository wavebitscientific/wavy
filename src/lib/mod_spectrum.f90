!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.
!
module mod_spectrum

use,intrinsic :: iso_c_binding,only:c_int,c_float
use mod_precision,only:ik => intkind,rk => realkind
use mod_utility,only:diff,diff_periodic
use mod_const,only:eps,pi,twopi,stderr,stdout
use datetime_module,only:datetime,timedelta
use json_module,only:json_core,json_file,json_value 

implicit none

private

public :: spectrum_type

type :: spectrum_type

  !! Spectrum class.

  private

  type(datetime) :: start_time !! Simulation start time
  type(datetime) :: end_time   !! Simulation end time
  type(timedelta) :: time_step !! Time step [s]

  real(kind=rk),dimension(:,:),allocatable :: spec !! 2-d spectrum

  real(kind=rk),dimension(:),allocatable :: f   !! Frequency [Hz]
  real(kind=rk),dimension(:),allocatable :: df  !! Frequency spacing [Hz]
  real(kind=rk),dimension(:),allocatable :: k   !! Wavenumber [rad/m]
  real(kind=rk),dimension(:),allocatable :: dk  !! Wavenumber spacing [rad/m]
  real(kind=rk),dimension(:),allocatable :: th  !! Direction [rad]
  real(kind=rk),dimension(:),allocatable :: dth !! Directional spacing [rad]
  real(kind=rk),dimension(:),allocatable :: cp  !! Phase speed [m/s]
  real(kind=rk),dimension(:),allocatable :: cg  !! Group speed [m/s]
  real(kind=rk),dimension(:),allocatable :: u   !! Mean current velocity in x-direction [m/s]
  real(kind=rk),dimension(:),allocatable :: v   !! Mean current velocity in y-direction [m/s]
  real(kind=rk),dimension(:),allocatable :: z   !! Depth levels for current array [m]

  real(kind=rk) :: air_density     !! Air density [kg/m^3]
  real(kind=rk) :: depth           !! Mean water depth [m]
  real(kind=rk) :: elevation       !! Mean surface elevation [m]
  real(kind=rk) :: grav            !! Gravitational acceleration [m/s^2]
  real(kind=rk) :: surface_tension !! Surface tension [N/m]
  real(kind=rk) :: water_density   !! Water density [kg/m^3]

  contains

  ! Public type-bound methods
  procedure,public,pass(self) :: frequencyMoment
  procedure,public,pass(self) :: getAirDensity
  procedure,public,pass(self) :: getAmplitude
  procedure,public,pass(self) :: getCurrent_u
  procedure,public,pass(self) :: getCurrent_v
  procedure,public,pass(self) :: getDepth
  procedure,public,pass(self) :: getDepthLevels
  procedure,public,pass(self) :: getDirections
  procedure,public,pass(self) :: getDirections2d
  procedure,public,pass(self) :: getElevation
  procedure,public,pass(self) :: getFrequency
  procedure,public,pass(self) :: getFrequency2d
  procedure,public,pass(self) :: getGravity
  procedure,public,pass(self) :: getGroupSpeed
  procedure,public,pass(self) :: getGroupSpeed2d
  procedure,public,pass(self) :: getWaveAction
  procedure,public,pass(self) :: getWavelength
  procedure,public,pass(self) :: getWavenumber
  procedure,public,pass(self) :: getWavenumberSpacing
  procedure,public,pass(self) :: getWavenumber2d
  procedure,public,pass(self) :: getPhaseSpeed
  procedure,public,pass(self) :: getPhaseSpeed2d
  procedure,public,pass(self) :: getSpectrum
  procedure,public,pass(self) :: getSurfaceTension
  procedure,public,pass(self) :: getWaterDensity
  procedure,public,pass(self) :: isAllocated
  procedure,public,pass(self) :: isMonochromatic
  procedure,public,pass(self) :: isOmnidirectional
  procedure,public,pass(self) :: meanPeriod
  procedure,public,pass(self) :: meanPeriodZeroCrossing
  procedure,public,pass(self) :: meanSquareSlope
  procedure,public,pass(self) :: meanSquareSlopeDirectional
  procedure,public,pass(self) :: momentum_x
  procedure,public,pass(self) :: momentum_y
  procedure,public,pass(self) :: momentumFlux_xx
  procedure,public,pass(self) :: momentumFlux_xy
  procedure,public,pass(self) :: momentumFlux_yy
  procedure,public,pass(self) :: omnidirectionalSpectrum
  procedure,public,pass(self) :: peakedness
  procedure,public,pass(self) :: peakFrequency
  procedure,public,pass(self) :: peakFrequencyDiscrete
  procedure,public,pass(self) :: saturationSpectrum
  procedure,public,pass(self) :: setAirDensity
  procedure,public,pass(self) :: setCurrent1d
  procedure,public,pass(self) :: setCurrent2d
  procedure,public,pass(self) :: setDepth
  procedure,public,pass(self) :: setElevation
  procedure,public,pass(self) :: setGravity
  procedure,public,pass(self) :: setSurfaceTension
  procedure,public,pass(self) :: setWaterDensity
  procedure,public,pass(self) :: setSpectrum1d
  procedure,public,pass(self) :: setSpectrum2d
  procedure,public,pass(self) :: significantWaveHeight
  procedure,public,pass(self) :: significantSurfaceOrbitalVelocity
  procedure,public,pass(self) :: stokesDrift
  procedure,public,pass(self) :: stokesDrift2d
  procedure,public,pass(self) :: ursellNumber
  procedure,public,pass(self) :: wavenumberMoment
  procedure,public,pass(self) :: wavenumberSpectrum
  procedure,public,pass(self) :: readJSON
  procedure,public,pass(self) :: writeJSON

  ! Private methods used to overload arithmetic operators
  procedure,private,pass(self) :: assign_array_1d
  procedure,private,pass(self) :: assign_array_2d
  procedure,private,pass(self) :: real_add_spectrum
  procedure,private,pass(self) :: real_sub_spectrum
  procedure,private,pass(self) :: real_mult_spectrum
  procedure,private,pass(self) :: real_div_spectrum
  procedure,private,pass(self) :: real2d_mult_spectrum
  procedure,private,pass(self) :: spectrum_add_spectrum
  procedure,private,pass(self) :: spectrum_add_real
  procedure,private,pass(self) :: spectrum_sub_spectrum
  procedure,private,pass(self) :: spectrum_sub_real
  procedure,private,pass(self) :: spectrum_mult_spectrum
  procedure,private,pass(self) :: spectrum_mult_real
  procedure,private,pass(self) :: spectrum_mult_real2d
  procedure,private,pass(self) :: spectrum_div_spectrum
  procedure,private,pass(self) :: spectrum_div_real
  procedure,private,pass(self) :: spectrum_unary_minus
  procedure,private,pass(self) :: eq
  procedure,private,pass(self) :: neq
  procedure,private,pass(self) :: gt
  procedure,private,pass(self) :: ge
  procedure,private,pass(self) :: lt
  procedure,private,pass(self) :: le

  ! Generic procedures
  generic,public :: setCurrent => setCurrent1d,setCurrent2d
  generic,public :: setSpectrum => setSpectrum1d,setSpectrum2d

  ! Generic operators
  generic :: assignment(=) => assign_array_1d,&
                              assign_array_2d
  generic :: operator(+)   => spectrum_add_spectrum,&
                              spectrum_add_real,&
                              real_add_spectrum
  generic :: operator(-)   => spectrum_sub_spectrum,&
                              spectrum_sub_real,&
                              real_sub_spectrum,&
                              spectrum_unary_minus
  generic :: operator(*)   => spectrum_mult_spectrum,&
                              spectrum_mult_real,&
                              spectrum_mult_real2d,&
                              real_mult_spectrum,&
                              real2d_mult_spectrum
  generic :: operator(/)   => spectrum_div_spectrum,&
                              spectrum_div_real,&
                              real_div_spectrum
  generic :: operator(==)  => eq
  generic :: operator(/=)  => neq
  generic :: operator(>)   => gt
  generic :: operator(>=)  => ge
  generic :: operator(<)   => lt
  generic :: operator(<=)  => le

endtype spectrum_type

interface spectrum_type
  module procedure :: constructor
endinterface spectrum_type

contains

!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function constructor(fmin,fmax,df,ndirs,&
  depth,grav,air_density,water_density,surface_tension) result(spectrum)

  !! Constructor function for the spectrum object.

  real(kind=rk),intent(in) :: fmin
    !! Minimum frequency bin [Hz]
  real(kind=rk),intent(in) :: fmax
    !! Maximum frequency bin [Hz]
  real(kind=rk),intent(in) :: df
    !! Frequency increment, df = f(n+1)/f(n)
  integer,intent(in) :: ndirs
    !! Number of directional bins
  real(kind=rk),intent(in) :: depth
    !! Mean water depth [m]
  real(kind=rk),intent(in),optional :: grav
    !! Gravitational acceleration [m/s^2]
  real(kind=rk),intent(in),optional :: air_density
    !! Air density [kg/m^3]
  real(kind=rk),intent(in),optional :: water_density
    !! Water density [kg/m^3]
  real(kind=rk),intent(in),optional :: surface_tension
    !! Surface tension [N/m]

  integer :: n
  integer :: nfreqs

  if(present(grav))then
    spectrum % grav = grav
  else
    spectrum % grav = 9.8_rk
  endif

  if(present(air_density))then
    spectrum % air_density = air_density
  else
    spectrum % air_density = 1.2_rk
  endif

  if(present(water_density))then
    spectrum % water_density = water_density
  else
    spectrum % water_density = 1e3_rk
  endif

  if(present(surface_tension))then
    spectrum % surface_tension = surface_tension
  else
    spectrum % surface_tension = 0.07_rk
  endif

  spectrum % depth = depth

  if(fmin == fmax)then
    !! monochromatic
    nfreqs = 1
  else
    nfreqs = int((log(fmax)-log(fmin))/log(df))
  endif

  allocate(spectrum % spec(nfreqs,ndirs))
  spectrum % spec = 0

  allocate(spectrum % f(nfreqs))
  allocate(spectrum % df(nfreqs))
  allocate(spectrum % k(nfreqs))
  allocate(spectrum % dk(nfreqs))
  allocate(spectrum % cp(nfreqs))
  allocate(spectrum % cg(nfreqs))
  allocate(spectrum % th(ndirs))

  if(nfreqs == 1)then
    spectrum % f(n) = fmin
  else
    do concurrent(n = 1:nfreqs)
      spectrum % f(n) = exp(log(fmin)+(n-1)*log(df))
    enddo
  endif

  do concurrent(n = 1:nfreqs)
    spectrum % k(n) = wavenumber(spectrum % f(n),         &
                                 spectrum % depth,        &
                                 spectrum % water_density,&
                                 spectrum % grav,         &
                                 spectrum % surface_tension)
  enddo

  spectrum % cp = twopi * spectrum % f / spectrum % k
  spectrum % cg = twopi * diff(spectrum % f) / diff(spectrum % k)

  do concurrent(n = 1:ndirs)
    spectrum % th(n) = (n-0.5*(ndirs+1))*twopi/ndirs
  enddo

  spectrum % df = diff(spectrum % f)
  spectrum % dk = diff(spectrum % k)

  if(ndirs > 1)then
    spectrum % dth = diff_periodic(spectrum % th)
    spectrum % dth(1) = spectrum % dth(2)
    spectrum % dth(ndirs) = spectrum % dth(ndirs-1)
  else
    spectrum % dth = [1]
  endif

  call spectrum % setCurrent2d([0._rk],[0._rk],[0._rk])

endfunction constructor
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function getAirDensity(self) result(air_density)
  !! Returns the air_density [kg/m^3] of the `spectrum` instance.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk) :: air_density !! Air density [kg/m^3]
  air_density = self % air_density
endfunction getAirDensity
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental logical function isAllocated(self)
  !! Returns the allocation status of the spectrum array.
  class(spectrum_type),intent(in) :: self !! `domain` instance
  isAllocated = allocated(self % spec)
endfunction isAllocated
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function isMonochromatic(self)
  !! Returns `.true.` if only one frequency bin is allocated,
  !! and `.false.` otherwise.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  logical :: isMonochromatic !! return value (boolean)
  if(size(self % f) == 1)then
    isMonochromatic = .true.
  else
    isMonochromatic = .false.
  endif
endfunction isMonochromatic
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function isOmnidirectional(self)
  !! Returns `.true.` if only one direction bin is allocated,
  !! and `.false.` otherwise.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  logical :: isOmnidirectional !! return value (boolean)
  if(size(self % th) == 1)then
    isOmnidirectional = .true.
  else
    isOmnidirectional = .false.
  endif
endfunction isOmnidirectional
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getFrequency(self) result(f)
  !! Returns the frequency [Hz] array of the spectrum instance.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:),allocatable :: f !! Frequency [Hz]
  f = self % f
endfunction getFrequency
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getFrequency2d(self) result(f)
  !! Returns the frequency [Hz] array of the spectrum instance, reshaped to
  !! match the spectrum array shape. This method is most useful for conforming
  !! shape array in 2-d spectrum computations.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:,:),allocatable :: f !! Frequency [Hz]
  integer :: ndirs,nfreqs
  integer :: ndir
  nfreqs = size(self % f)
  ndirs = size(self % th)
  allocate(f(nfreqs,ndirs))
  do concurrent(ndir = 1:ndirs)
    f(:,ndir) = self % f
  enddo
endfunction getFrequency2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getWavenumber(self) result(k)
  !! Returns the wavenumber [rad/m] array of the spectrum instance.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:),allocatable :: k !! Wavenumber [rad/m]
  k = self % k
endfunction getWavenumber
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getWavenumberSpacing(self) result(dk)
  !! Returns the wavenumber spacing [rad/m] array of the spectrum instance.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:),allocatable :: dk !! Wavenumber spacing [rad/m]
  dk = self % dk
endfunction getWavenumberSpacing
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getWavenumber2d(self) result(k)
  !! Returns the wavenumber [rad/m] array of the spectrum instance, reshaped to
  !! match the spectrum array shape. This method is most useful for conforming
  !! shape array in 2-d spectrum computations.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:,:),allocatable :: k !! Wavenumber [rad/m]
  integer :: ndirs,nfreqs
  integer :: ndir
  nfreqs = size(self % f)
  ndirs = size(self % th)
  allocate(k(nfreqs,ndirs))
  do concurrent(ndir = 1:ndirs)
    k(:,ndir) = self % k
  enddo
endfunction getWavenumber2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getWavelength(self) result(lambda)
  !! Returns the wavelength [m] array of the spectrum instance.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:),allocatable :: lambda !! Wavelength [m]
  lambda = twopi / self % k
endfunction getWavelength
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getDirections(self) result(th)
  !! Returns the directions [rad] array of the spectrum instance.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:),allocatable :: th !! Directions [rad]
  th = self % th
endfunction getDirections
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getDirections2d(self) result(th)
  !! Returns the directions [rad] array of the spectrum instance, reshaped to
  !! match the spectrum array shape. This method is most useful for conforming
  !! shape array in 2-d spectrum computations.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:,:),allocatable :: th !! Directions [rad]
  integer :: ndirs,nfreqs
  integer :: nfreq
  nfreqs = size(self % f)
  ndirs = size(self % th)
  allocate(th(nfreqs,ndirs))
  do concurrent(nfreq = 1:nfreqs)
    th(nfreq,:) = self % th
  enddo
endfunction getDirections2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getPhaseSpeed(self) result(cp)
  !! Returns the phase speed [m/s] array of the spectrum instance.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:),allocatable :: cp !! Phase speed [m/s]
  cp = self % cp
endfunction getPhaseSpeed
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getPhaseSpeed2d(self) result(cp)
  !! Returns the phase speed [m/s] array of the spectrum instance.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:,:),allocatable :: cp !! Phase speed [m/s]
  integer :: ndirs,nfreqs
  integer :: ndir
  nfreqs = size(self % f)
  ndirs = size(self % th)
  allocate(cp(nfreqs,ndirs))
  do concurrent(ndir = 1:ndirs)
    cp(:,ndir) = self % cp
  enddo
endfunction getPhaseSpeed2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getGroupSpeed(self) result(cg)
  !! Returns the phase speed [m/s] array of the spectrum instance.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:),allocatable :: cg !! Group speed [m/s]
  cg = self % cg
endfunction getGroupSpeed
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getGroupSpeed2d(self) result(cg)
  !! Returns the group speed [m/s] array of the spectrum instance.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:,:),allocatable :: cg !! Group speed [m/s]
  integer :: ndirs,nfreqs
  integer :: ndir
  nfreqs = size(self % f)
  ndirs = size(self % th)
  allocate(cg(nfreqs,ndirs))
  do concurrent(ndir = 1:ndirs)
    cg(:,ndir) = self % cg
  enddo
endfunction getGroupSpeed2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getSpectrum(self) result(spec)
  !! Returns the spectrum array.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:,:),allocatable :: spec !! Spectrum array
  spec = self % spec
endfunction getSpectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getWaveAction(self) result(wave_action)
  !! Returns the wave action spectrum, which corresponds to the the wave
  !! variance spectrum normalized by the intrinsic frequency.
  class(spectrum_type),intent(in) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:,:),allocatable :: wave_action
    !! Wave action array
  wave_action = self % spec / self % getFrequency2d()
endfunction getWaveAction
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getAmplitude(self) result(a)
  !! Returns the amplitude array.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:,:),allocatable :: a !! Amplitude [m]
  integer :: ndir,ndirs
  if(self % isMonochromatic())then
    a = sqrt(2*self % spec)
  else
    ndirs = size(self % getDirections())
    do concurrent(ndir = 1:ndirs)
      a(:,ndir) = sqrt(2*self % spec(:,ndir)*self % df)
    enddo
  endif
endfunction getAmplitude
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getCurrent_u(self) result(u)
  !! Returns the current velocity in x-direction.
  class(spectrum_type),intent(in) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:),allocatable :: u
    !! Mean current velocity in x-direction [m/s]
  u = self % u
endfunction getCurrent_u
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getCurrent_v(self) result(v)
  !! Returns the current velocity in y-direction.
  class(spectrum_type),intent(in) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:),allocatable :: v
    !! Mean current velocity in y-direction [m/s]
  v = self % v
endfunction getCurrent_v
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function getDepthLevels(self) result(z)
  !! Returns the depth levels at which the current arrays are defined.
  class(spectrum_type),intent(in) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:),allocatable :: z
    !! Depth levels of current fields [m]
  z = self % z
endfunction getDepthLevels
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function getDepth(self) result(depth)
  !! Returns the mean water depth [m].
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk) :: depth !! Mean water depth [m]
  depth = self % depth
endfunction getDepth
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function getElevation(self) result(elevation)
  !! Returns the mean surface elevation anomaly [m].
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk) :: elevation !! Mean surface elevation anomaly [m]
  elevation = self % elevation
endfunction getElevation
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function getGravity(self) result(grav)
  !! Returns the gravitational acceleration [m/s^2].
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk) :: grav !! Gravitational acceleration [m/s^2]
  grav = self % grav
endfunction getGravity
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function getSurfaceTension(self) result(surface_tension)
  !! Returns the surface tension [N/m].
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk) :: surface_tension !! Surface tension [N/m]
  surface_tension = self % surface_tension
endfunction getSurfaceTension
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function getWaterDensity(self) result(water_density)
  !! Returns the water density [kg/m^3].
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk) :: water_density !! Water density [kg/m^3]
  water_density = self % water_density
endfunction getWaterDensity
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function omnidirectionalSpectrum(self) result(spec)
  !! Returns the omnidirectional spectrum that corresponds to the input
  !! directional spectrum, integrated over all directions.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:),allocatable :: spec !! Spectrum array
  integer(kind=ik) :: ndir,ndirs,nfreqs
  nfreqs = size(self % spec,dim=1)
  ndirs = size(self % spec,dim=2)
  allocate(spec(nfreqs))
  spec = 0
  do ndir = 1,ndirs
    spec(:) = spec(:)+self % spec(:,ndir)*self % dth(ndir)
  enddo
endfunction omnidirectionalSpectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure subroutine setSpectrum1d(self,spec)
  !! Sets the 2-d spectrum array. This procedure is overloaded by the
  !! generic procedure setSpectrum.
  class(spectrum_type),intent(inout) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:),intent(in) :: spec
    !! Input 1-d spectrum array
  self % spec(:,1) = spec
endsubroutine setSpectrum1d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure subroutine setSpectrum2d(self,spec)
  !! Sets the 2-d spectrum array. This procedure is overloaded by the
  !! generic procedure setSpectrum.
  class(spectrum_type),intent(inout) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:,:),intent(in) :: spec
    !! Input 2-d spectrum array
  self % spec = spec
endsubroutine setSpectrum2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure subroutine setCurrent1d(self,u,z)
  !! Sets the 1-d current velocity field. This procedure is overloaded by the
  !! generic procedure setCurrent.
  class(spectrum_type),intent(inout) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:),intent(in) :: u
    !! Current velocity in x-direction [m/s]
  real(kind=rk),dimension(:),intent(in) :: z
    !! Depth levels for the velocity array [m]
  self % u = u
  self % z = z
  allocate(self % v(size(u)))
  self % v = 0
endsubroutine setCurrent1d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure subroutine setCurrent2d(self,u,v,z)
  !! Sets the 2-d current velocity field. This procedure is overloaded by the
  !! generic procedure setCurrent.
  class(spectrum_type),intent(inout) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:),intent(in) :: u
    !! Current velocity in x-direction [m/s]
  real(kind=rk),dimension(:),intent(in) :: v
    !! Current velocity in y-direction [m/s]
  real(kind=rk),dimension(:),intent(in) :: z
    !! Depth levels for the velocity array [m]
  self % u = u
  self % v = v
  self % z = z
endsubroutine setCurrent2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental subroutine setDepth(self,depth)
  !! Sets the mean surface elevation value.
  class(spectrum_type),intent(inout) :: self !! Spectrum instance
  real(kind=rk),intent(in) :: depth !! Mean water depth [m]
  self % depth = depth
endsubroutine setDepth
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental subroutine setElevation(self,elevation)
  !! Sets the mean surface elevation value.
  class(spectrum_type),intent(inout) :: self
    !! Spectrum instance
  real(kind=rk),intent(in) :: elevation
    !! Mean surface elevation anomaly [m]
  self % elevation = elevation
endsubroutine setElevation
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental subroutine setGravity(self,grav)
  !! Sets the gravitational acceleration [m/s^2].
  class(spectrum_type),intent(inout) :: self !! Spectrum instance
  real(kind=rk),intent(in) :: grav !! Gravitational acceleration [m/s^2]
  self % grav = grav
endsubroutine setGravity
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental subroutine setSurfaceTension(self,surface_tension)
  !! Sets the surface tension [N/m].
  class(spectrum_type),intent(inout) :: self !! Spectrum instance
  real(kind=rk),intent(in) :: surface_tension !! Surface tension [N/m]
  self % surface_tension = surface_tension
endsubroutine setSurfaceTension
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental subroutine setAirDensity(self,air_density)
  !! Sets the air density [kg/m^3].
  class(spectrum_type),intent(inout) :: self !! Spectrum instance
  real(kind=rk),intent(in) :: air_density !! Air density [kg/m^3]
  self % air_density = air_density
endsubroutine setAirDensity
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental subroutine setWaterDensity(self,water_density)
  !! Sets the water density [kg/m^3].
  class(spectrum_type),intent(inout) :: self !! Spectrum instance
  real(kind=rk),intent(in) :: water_density !! Water density [kg/m^3]
  self % water_density = water_density
endsubroutine setWaterDensity
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function meanSquareSlope(self)
  !! Returns the mean square slope of the spectrum, which is the second
  !! moment of the wavenumber spectrum.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  meanSquareSlope = self % wavenumberMoment(2)
endfunction meanSquareSlope
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function meanSquareSlopeDirectional(self) result(mss)

  !! For each directional frequency bin, computes the mean square slope of all
  !! all waves longer than that bin, projected to the direction of that bin.

  class(spectrum_type),intent(in) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:,:),allocatable :: mss
    !! Directional mean square slope

  real(kind=rk),dimension(:,:),allocatable :: wavenumber_spectrum
  real(kind=rk),dimension(:,:),allocatable :: dir_projection

  integer :: nfreq,nfreqs
  integer :: ndir,ndirs

  nfreqs = size(self % f)
  ndirs = size(self % th)

  wavenumber_spectrum = self % wavenumberSpectrum()

  ! Compute projection of each wave direction onto every other direction
  allocate(dir_projection(ndirs,ndirs))
  do concurrent(ndir = 1:ndirs)
    dir_projection(:,ndir) = abs(cos(self % th(ndir)-self % th(:)))
  enddo

  allocate(mss(nfreqs,ndirs))
  mss = 0

  do ndir = 1,ndirs
    do nfreq = 2,nfreqs
      mss(nfreq,ndir) = mss(nfreq-1,ndir)                           &
        + sum(wavenumber_spectrum(nfreq-1,:)*dir_projection(:,ndir))&
        * self % k(nfreq-1)**2 * self % dk(nfreq-1)
    enddo
  enddo

  deallocate(dir_projection,wavenumber_spectrum)

endfunction meanSquareSlopeDirectional
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function momentum_x(self)
  !! Returns total wave momentum [kg/m/s] in x-direction.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  integer(kind=ik) :: n
  integer(kind=ik) :: nfreqs
  nfreqs = size(self % f)
  momentum_x = 0
  do n = 1,nfreqs
    momentum_x = momentum_x + sum(self % spec(n,:)*self % dth*cos(self % th))&
      * self % df(n) / self % cp(n)
  enddo
  momentum_x = momentum_x * self % water_density * self % grav
endfunction momentum_x
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function momentum_y(self)
  !! Returns total wave momentum [kg/m/s] in y-direction.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  integer(kind=ik) :: n
  integer(kind=ik) :: nfreqs
  nfreqs = size(self % f)
  momentum_y = 0
  do n = 1,nfreqs
    momentum_y = momentum_y + sum(self % spec(n,:)*self % dth*sin(self % th))&
      * self % df(n) / self % cp(n)
  enddo
  momentum_y = momentum_y * self % water_density * self % grav
endfunction momentum_y
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function momentumFlux_xx(self)
  !! Returns total advective flux [kg/m^2/s^2] in y-direction of momentum in
  !! y-direction.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  integer(kind=ik) :: n
  integer(kind=ik) :: nfreqs
  nfreqs = size(self % f)
  momentumFlux_xx = 0
  do n = 1,nfreqs
    momentumFlux_xx = momentumFlux_xx                     &
      + sum(self % spec(n,:)*self % dth*cos(self % th)**2)&
      * self % df(n) * self % cg(n) / self % cp(n)
  enddo
  momentumFlux_xx = momentumFlux_xx * self % water_density * self % grav
endfunction momentumFlux_xx
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function momentumFlux_xy(self)
  !! Returns total advective flux [kg/m^2/s^2] in x-direction of momentum in
  !! y-direction and vice versa (flux in y-direction of momentum in
  !! y-direction), because \int{Cgx*My} == \int{Cgy*Mx}.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  integer(kind=ik) :: n
  integer(kind=ik) :: nfreqs
  nfreqs = size(self % f)
  momentumFlux_xy = 0
  do n = 1,nfreqs
    momentumFlux_xy = momentumFlux_xy                                 &
      + sum(self % spec(n,:)*self % dth*cos(self % th)*sin(self % th))&
      * self % df(n) * self % cg(n) / self % cp(n)
  enddo
  momentumFlux_xy = momentumFlux_xy * self % water_density * self % grav
endfunction momentumFlux_xy
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function momentumFlux_yy(self)
  !! Returns total advective flux [kg/m^2/s^2] in y-direction of momentum in
  !! y-direction.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  integer(kind=ik) :: n
  integer(kind=ik) :: nfreqs
  nfreqs = size(self % f)
  momentumFlux_yy = 0
  do n = 1,nfreqs
    momentumFlux_yy = momentumFlux_yy                     &
      + sum(self % spec(n,:)*self % dth*sin(self % th)**2)&
      * self % df(n) * self % cg(n) / self % cp(n)
  enddo
  momentumFlux_yy = momentumFlux_yy * self % water_density * self % grav
endfunction momentumFlux_yy
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function frequencyMoment(self,n)
  !! Returns the spectral frequency moment of order n.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  integer,intent(in) :: n !! Moment order
  !frequencyMoment = sum(self % f**n*sum(self % spec,dim=2)*self % df)
  frequencyMoment = sum(self % f**n*self % omnidirectionalSpectrum()*self % df)
endfunction frequencyMoment
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function peakedness(self)
  !! Returns the peakedness parameter that quantifies the sharpness of the
  !! spectral peak, following Goda (1970).
  !!
  !! References:
  !!
  !! Goda, Y., 1970. Numerical experiments on waves statistics with spectral
  !! simulation. *Report. Port and Harbour Research Institute*, Japan, **9**,
  !! 3-57.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  peakedness = 2*sum(self % f*self % omnidirectionalSpectrum()**2*self % df)&
             / (self % frequencyMoment(0)**2+eps)
endfunction peakedness
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function peakFrequency(self)
  !! Returns the peak frequency based on Young (1995).
  !!
  !! References:
  !!
  !! Young, I, 1995. The determination of confidence limits associated with
  !! estimates of the spectral peak frequency. *Ocean Engng.*, **22**, 669-686.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  peakFrequency = sum(self % f*self % omnidirectionalSpectrum()**4*self % df)&
                / (sum(self % omnidirectionalSpectrum()**4*self % df)+eps)
endfunction peakFrequency
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function peakFrequencyDiscrete(self)
  !! Returns the peak frequency based on simple discrete maximum location of
  !! the spectrum array.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  peakFrequencyDiscrete = &
    self % f(maxloc(self % omnidirectionalSpectrum(),dim=1))
endfunction peakFrequencyDiscrete
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function wavenumberMoment(self,n)
  !! Returns the spectral wavenumber moment of order n.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  integer,intent(in) :: n !! Moment order
  wavenumberMoment = sum(self % k**n*sum(self % wavenumberSpectrum(),dim=2)&
                                    *self % dk)
endfunction wavenumberMoment
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function wavenumberSpectrum(self) result(spec)
  !! Returns the wavenumber spectrum array of the spectrum instance.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  real(kind=rk),dimension(:,:),allocatable :: spec !! Spectrum array
  integer :: nfreq,nfreqs
  integer :: ndirs
  nfreqs = size(self % f)
  ndirs = size(self % th)
  allocate(spec(nfreqs,ndirs))
  do concurrent(nfreq = 1:nfreqs)
    spec(nfreq,:) = self % spec(nfreq,:) * self % cg(nfreq) / twopi
  enddo
endfunction wavenumberSpectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function saturationSpectrum(self)
  !! Returns the saturation spectrum B(k) = F(k)k^4.
  class(spectrum_type),intent(in) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:,:),allocatable :: saturationSpectrum
    !! Saturation spectrum result
  real(kind=rk),dimension(:,:),allocatable :: wavenumber_spectrum
  integer :: nfreqs,ndirs
  integer :: n
  nfreqs = size(self % f)
  ndirs = size(self % th)
  allocate(saturationSpectrum(nfreqs,ndirs))
  wavenumber_spectrum = self % wavenumberSpectrum()
  do concurrent(n=1:ndirs)
    saturationSpectrum(:,n) = wavenumber_spectrum(:,n) * self % k**4
  enddo
endfunction saturationSpectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function significantWaveHeight(self)
  !! Returns the significant wave height [m].
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  significantWaveHeight = 4*sqrt(self % frequencyMoment(0))
endfunction significantWaveHeight
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function &
  significantSurfaceOrbitalVelocity(self) result(uorb)
  !! Returns the significant surface orbital velocity [m/s].
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  uorb = 2*sqrt(sum((twopi*self % f)**2*sum(self % spec,dim=2)*self % df))
endfunction significantSurfaceOrbitalVelocity
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function stokesDrift(self,z)
  !! Exact solution of Stokes drift based on linear wave theory, given input 
  !! omnidirectional spectrum and distance from surface `z` [m], negative 
  !! downward.
  class(spectrum_type),intent(in) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:),intent(in) :: z
    !! Distance from surface [m], negative downward
  real(kind=rk),dimension(:),allocatable :: stokesDrift
    !! Stokes drift array [m/s]
  integer(kind=ik) :: n
  allocate(stokesDrift(size(z)))
  do concurrent(n = 1:size(z))
    stokesDrift(n) = sum(self % omnidirectionalSpectrum()*self % k&
                        *exp(2*self % k*z(n))*self % df)
  enddo
endfunction stokesDrift
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function stokesDrift2d(self,z)
  !! Exact solution of Stokes drift based on linear wave theory, given input 
  !! directional spectrum and distance from surface `z` [m], negative downward.
  class(spectrum_type),intent(in) :: self
    !! Spectrum instance
  real(kind=rk),dimension(:),intent(in) :: z
    !! Distance from surface [m], negative downward
  real(kind=rk),dimension(:,:),allocatable :: stokesDrift2d
    !! Stokes drift array [m/s]
  integer(kind=ik) :: n,ndir,ndirs
  ndirs = size(self % getDirections())
  allocate(stokesDrift2d(size(z),2))
  stokesDrift2d = 0
  do n = 1,size(z)
    stokesDrift2d(n,:) = 0
    do ndir = 1,ndirs
      ! x-component of Stokes drift
      stokesDrift2d(n,1) = stokesDrift2d(n,1)&
                         + sum(self % spec(:,ndir)*cos(self % th(ndir))&
                              *self % k*exp(2*self % k*z(n))&
                              *self % df*self % dth(ndir))
      ! y-component of Stokes drift
      stokesDrift2d(n,2) = stokesDrift2d(n,2)&
                         + sum(self % spec(:,ndir)*sin(self % th(ndir))&
                              *self % k*exp(2*self % k*z(n))&
                              *self % df*self % dth(ndir))
    enddo
  enddo
endfunction stokesDrift2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function meanPeriod(self)
  !! Returns the mean wave period [s].
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  meanPeriod = self % frequencyMoment(0) / (self % frequencyMoment(1) + eps)
endfunction meanPeriod
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function meanPeriodZeroCrossing(self)
  !! Returns the zero-crossing mean wave period [s]:
  !!
  !! Tm02 = \sqrt(m_0 / m_2)
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  meanPeriodZeroCrossing = sqrt(self % frequencyMoment(0)&
                              / (self % frequencyMoment(2) + eps))
endfunction meanPeriodZeroCrossing
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function ursellNumber(self)
  !! Returns the Ursell number.
  class(spectrum_type),intent(in) :: self !! Spectrum instance
  ursellNumber = self % grav / (8*sqrt(2._rk)*pi**2)              &
               * self % significantWaveHeight() * self % meanPeriod()**2&
               / self % depth**2
endfunction ursellNumber
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure subroutine assign_array_1d(self,array)
  !! Assigns a 1-d array of reals to a `spectrum` instance. This procedure
  !! overloads the assignment ('=') operator.
  class(spectrum_type),intent(inout) :: self
    !! l.h.s. `spectrum` instance
  real(kind=rk),dimension(:),intent(in) :: array
    !! r.h.s. array of reals
  call self % setSpectrum(array)
endsubroutine assign_array_1d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure subroutine assign_array_2d(self,array)
  !! Assigns a 2-d array of reals to a `spectrum` instance. This procedure
  !! overloads the assignment ('=') operator.
  class(spectrum_type),intent(inout) :: self
    !! l.h.s. `spectrum` instance
  real(kind=rk),dimension(:,:),intent(in) :: array
    !! r.h.s. array of reals
  call self % setSpectrum(array)
endsubroutine assign_array_2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental logical function eq(self,s2)
  !! Logical equality comparison function. Overloads the `==` operator.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  class(spectrum_type),intent(in) :: s2 !! r.h.s. spectrum instance
  eq = all(self % getSpectrum() == s2 % getSpectrum())
endfunction eq
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental logical function neq(self,s2)
  !! Logical inequality comparison function. Overloads the `/=` operator.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  class(spectrum_type),intent(in) :: s2 !! r.h.s. spectrum instance
  neq = .not. self == s2
endfunction neq
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental logical function gt(self,s2)
  !! Logical greater than comparison function. Overloads the `>` operator.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  class(spectrum_type),intent(in) :: s2 !! r.h.s. spectrum instance
  gt = all(self % getSpectrum() > s2 % getSpectrum())
endfunction gt
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental logical function lt(self,s2)
  !! Logical less than comparison function. Overloads the `<` operator.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  class(spectrum_type),intent(in) :: s2 !! r.h.s. spectrum instance
  lt = all(self % getSpectrum() < s2 % getSpectrum())
endfunction lt
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental logical function ge(self,s2)
  !! Logical greater than or equal comparison function. Overloads the `>=` 
  !! operator.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  class(spectrum_type),intent(in) :: s2 !! r.h.s. spectrum instance
  ge = all(self % getSpectrum() >= s2 % getSpectrum())
endfunction ge
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental logical function le(self,s2)
  !! Logical less than or equal comparison function. Overloads the `<=` 
  !! operator.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  class(spectrum_type),intent(in) :: s2 !! r.h.s. spectrum instance
  le = all(self % getSpectrum() <= s2 % getSpectrum())
endfunction le
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function spectrum_add_spectrum(self,s2)&
  result(spec)
  !! Returns a spectrum instance with the spectrum array values being the sum of
  !! the two input spectrum instances.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  class(spectrum_type),intent(in) :: s2 !! r.h.s. spectrum instance
  spec = self
  spec = self % getSpectrum() + s2 % getSpectrum()
endfunction spectrum_add_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function spectrum_sub_spectrum(self,s2)&
  result(spec)
  !! Subtracts one spectrum instance from another and returns the resulting
  !! spectrum instance.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  class(spectrum_type),intent(in) :: s2 !! r.h.s. spectrum instance
  spec = self
  spec = self % getSpectrum() - s2 % getSpectrum()
endfunction spectrum_sub_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function spectrum_mult_spectrum(self,s2)&
  result(spec)
  !! Returns a product of two spectrum instances.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  class(spectrum_type),intent(in) :: s2 !! r.h.s. spectrum instance
  spec = self
  spec = self % getSpectrum() * s2 % getSpectrum()
endfunction spectrum_mult_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function spectrum_div_spectrum(self,s2)&
  result(spec)
  !! Returns a division of two spectrum instances.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  class(spectrum_type),intent(in) :: s2 !! r.h.s. spectrum instance
  spec = self
  spec = self % getSpectrum() / s2 % getSpectrum()
endfunction spectrum_div_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function spectrum_add_real(self,a)&
  result(spec)
  !! Returns a sum of a spectrum instance and a real number.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  real(kind=rk),intent(in) :: a !! r.h.s. real number
  spec = self
  spec = self % getSpectrum() + a
endfunction spectrum_add_real
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function spectrum_sub_real(self,a)&
  result(spec)
  !! Returns a difference of a spectrum instance and a real number.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  real(kind=rk),intent(in) :: a !! r.h.s. real number
  spec = self
  spec = self % getSpectrum() - a
endfunction spectrum_sub_real
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function spectrum_mult_real(self,a)&
  result(spec)
  !! Returns a product of a spectrum instance and a real number.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  real(kind=rk),intent(in) :: a !! r.h.s. real number
  spec = self
  spec = self % getSpectrum() * a
endfunction spectrum_mult_real
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure type(spectrum_type) function spectrum_mult_real2d(self,a)&
  result(spec)
  !! Returns a product of a spectrum instance and a real number.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  real(kind=rk),dimension(:,:),intent(in) :: a !! r.h.s. real 2-d array
  spec = self
  spec = self % getSpectrum() * a
endfunction spectrum_mult_real2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function spectrum_div_real(self,a) result(spec)
  !! Returns a division of a spectrum instance and a real number.
  class(spectrum_type),intent(in) :: self !! l.h.s. spectrum instance
  real(kind=rk),intent(in) :: a !! r.h.s. real number
  spec = self
  spec = self % getSpectrum() / a
endfunction spectrum_div_real
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function real_add_spectrum(a,self)&
  result(spec)
  !! Returns a sum of a real number and a spectrum instance.
  real(kind=rk),intent(in) :: a !! l.h.s. real number
  class(spectrum_type),intent(in) :: self !! r.h.s. spectrum instance
  spec = self
  spec = self % getSpectrum() + a
endfunction real_add_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function real_sub_spectrum(a,self)&
  result(spec)
  !! Returns a difference between a real number and a spectrum instance.
  real(kind=rk),intent(in) :: a !! l.h.s. real number
  class(spectrum_type),intent(in) :: self !! r.h.s. spectrum instance
  spec = self
  spec = a - self % getSpectrum()
endfunction real_sub_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function real_mult_spectrum(a,self)&
  result(spec)
  !! Returns a product of a real number and a spectrum instance.
  real(kind=rk),intent(in) :: a !! l.h.s. real number
  class(spectrum_type),intent(in) :: self !! r.h.s. spectrum instance
  spec = self
  spec = self % getSpectrum() * a
endfunction real_mult_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure type(spectrum_type) function real2d_mult_spectrum(a,self)&
  result(spec)
  !! Returns a product of a real number and a spectrum instance.
  real(kind=rk),dimension(:,:),intent(in) :: a !! l.h.s. real 2-d array
  class(spectrum_type),intent(in) :: self !! r.h.s. spectrum instance
  spec = self
  spec = self % getSpectrum() * a
endfunction real2d_mult_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function real_div_spectrum(a,self) result(spec)
  !! Returns a division of a real number and a spectrum instance.
  real(kind=rk),intent(in) :: a !! l.h.s. real number
  class(spectrum_type),intent(in) :: self !! r.h.s. spectrum instance
  spec = self
  spec = a / self % getSpectrum()
endfunction real_div_spectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental type(spectrum_type) function spectrum_unary_minus(self)&
  result(spec)
  !! Returns a negative value of the spectrum instance.
  class(spectrum_type),intent(in) :: self !! r.h.s. spectrum instance
  spec = self
  spec = - self % getSpectrum()
endfunction spectrum_unary_minus
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function wavenumber(f,depth,water_density,grav,surface_tension)

  !! Solves the linear water wave dispersion relationship using a
  !! Newton-Raphson iteration loop.

  real(kind=rk),intent(in) :: f
    !! Intrinsic frequency [Hz]
  real(kind=rk),optional,intent(in) :: depth
    !! Mean water depth [m]
  real(kind=rk),optional,intent(in) :: water_density
    !! Water density [kg/m^3]
  real(kind=rk),optional,intent(in) :: grav
    !! Gravitational acceleration [m/s^2]
  real(kind=rk),optional,intent(in) :: surface_tension
    !! Surface tension [N/m]

  real(kind=rk) :: wavenumber,dk,b,fnd,t
  integer :: counter

  associate(k => wavenumber)

  fnd = twopi*f*sqrt(depth/grav)
  k = fnd**2

  b = surface_tension/(water_density*grav*depth**2)

  counter = 1
  dk = 2e-3_rk
  newton_raphson: do
    t = tanh(k)
    dk = -(fnd**2-k*t*(1+b*k**2))&
         /(3*b*k**2*t+t+k*(1+b*k**2)*(1-t**2))
    k = k-dk
    if(abs(dk) < eps .or. counter > 100)then
      exit newton_raphson
    endif
    counter = counter+1
  enddo newton_raphson

  k = k / depth

  endassociate

endfunction wavenumber
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine readJSON(self,filename)
  !! Read a spectrum instance from a JSON file.
  class(spectrum_type),intent(inout) :: self !! `spectrum` instance
  character(len=*),intent(in) :: filename !! JSON file name
  type(json_file) :: json
  logical :: found 
  integer(kind=ik) :: nfreqs
  integer(kind=ik) :: ndirs
  real(kind=rk),dimension(:),allocatable :: arr
  call json % initialize()
  call json % load_file(trim(filename))
  call json % get('frequency',self % f,found)
  call json % get('wavenumber',self % k,found)
  call json % get('directions',self % th,found)
  call json % get('spectrum',arr,found)
  nfreqs = size(self % f)
  ndirs = size(self % th)
  self % spec = reshape(arr,[nfreqs,ndirs])
  call json % get('depth',self % depth,found)
  call json % get('elevation',self % elevation,found)
  call json % get('gravity',self % grav,found)
  call json % get('air_density',self % air_density,found)
  call json % get('water_density',self % water_density,found)
  call json % get('surface_tension',self % surface_tension,found)
  call json % get('u-velocity',self % u,found)
  call json % get('v-velocity',self % v,found)
  call json % get('z',self % z,found)
  call json % destroy()
  self % df = diff(self % f)
  self % dk = diff(self % k)
  if(ndirs > 1)then
    self % dth = diff_periodic(self % th)
  else
    self % dth = [1]
  endif
  self % cp = twopi * self % f / self % k
  self % cg = twopi * self % df / self % dk
endsubroutine readJSON
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine writeJSON(self,filename,minify)
  !! Writes a spectrum instance to a JSON file.
  class(spectrum_type),intent(in) :: self !! `spectrum` instance
  character(len=*),intent(in) :: filename !! JSON file name
  logical,intent(in) :: minify !! Logical switch to minify the JSON file
  type(json_core) :: json
  type(json_value),pointer :: ptr
  call json % initialize(no_whitespace=minify,real_format='ES')
  call json % create_object(ptr,'')
  call json % add(ptr,'frequency',self % getFrequency())
  call json % add(ptr,'wavenumber',self % getWavenumber())
  call json % add(ptr,'directions',self % getDirections())
  call json % add(ptr,'spectrum',pack(self % getSpectrum(),.true.))
  call json % add(ptr,'depth',self % getDepth())
  call json % add(ptr,'elevation',self % getElevation())
  call json % add(ptr,'gravity',self % getGravity())
  call json % add(ptr,'air_density',self % getAirDensity())
  call json % add(ptr,'water_density',self % getWaterDensity())
  call json % add(ptr,'surface_tension',self % getSurfaceTension())
  call json % add(ptr,'u-velocity',self % getCurrent_u())
  call json % add(ptr,'v-velocity',self % getCurrent_v())
  call json % add(ptr,'z',self % getDepthLevels())
  call json % print(ptr,trim(filename))
  call json % destroy(ptr)
endsubroutine writeJSON
!-------------------------------------------------------------------------------
endmodule mod_spectrum
