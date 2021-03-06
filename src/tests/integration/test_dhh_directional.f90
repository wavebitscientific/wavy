!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3-clause license. See LICENSE for details.
!
program test_dhh_directional

use mod_const
use mod_utility
use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_spectral_shapes,only:donelanHamiltonHui,&
                             donelanHamiltonHuiDirectionalSpectrum,&
                             donelanHamiltonHuiDirectionalSpreading,&
                             jonswapPeakFrequency
use mod_nondimensional,only:nondimensionalFrequency
use mod_functional,only:arange

implicit none

type(spectrum_type) :: spec
real(kind=rk),dimension(:),allocatable :: wspd
real(kind=rk),dimension(:),allocatable :: fetch
real(kind=rk),dimension(:),allocatable :: dirs
real(kind=rk),dimension(:,:),allocatable :: spreading
real(kind=rk),parameter :: grav = 9.8_rk
integer :: m,n,ndir,ndirs,nfreqs

write(*,*)
write(*,*)'Test omnidirectional spectrum with Donelan et al. (1985) shape;'
write(*,*)'Vary wind speed from 5 to 60 m/s with 5 m/s increments.'
write(*,*)'Vary fetch between 1, 10, and 100 km.'
write(*,*)'Peak frequency based on JONSWAP.'
write(*,*)

! Initialize spectrum
spec = spectrum_type(fmin=0.0313_rk,fmax=2._rk,df=1.1_rk,ndirs=36,depth=1e3_rk)

dirs = spec % getDirections()

ndirs  = size(spec % getDirections())
nfreqs = size(spec % getFrequency())

allocate(spreading(nfreqs,ndirs))
spreading = 0

! Set array of wind speed and fetch values:
wspd = arange(5,60,5)
fetch = [1e3_rk,1e4_rk,1e5_rk]

do m = 1,size(fetch)

  write(*,*)
  write(*,fmt='(a,f5.1,a)')'fetch = ',fetch(m)/1e3_rk,' km'
  write(*,fmt='(a)')'   wspd      Hs      Tp       Tm1      Tm2      mss      m0(f)'&
                  //'    m1(f)    m2(f)'
  write(*,fmt='(a)')'----------------------------------------------------------'&
                  //'----------------------'
  do n = 1,size(wspd)

    ! DHH directional spectral shape:
    spec = donelanHamiltonHuiDirectionalSpectrum(f = spec % getFrequency(),&
                              theta = spec % getDirections(),&
                              wspd  = wspd(n),&
                              fpeak = jonswapPeakFrequency(wspd(n),fetch(m),grav),&
                              theta_mean = 0._rk,&
                              grav  = grav)

    write(*,fmt='(9(f8.4,1x))')wspd(n),spec % significantWaveHeight(),&
      1./spec % peakFrequency(),spec % meanPeriod(),&
      spec % meanPeriodZeroCrossing(),spec % meanSquareSlope(),&
      spec % frequencyMoment(0),spec % frequencyMoment(1),&
      spec % frequencyMoment(2)

  enddo

enddo

endprogram test_dhh_directional
