!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.
!
program test_phillips

use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_spectral_shapes,only:phillips,jonswapPeakFrequency
use mod_utility,only:range

implicit none

type(spectrum_type) :: spec
real(kind=rk),dimension(:),allocatable :: wspd
real(kind=rk),dimension(:),allocatable :: fetch
real(kind=rk),parameter :: grav = 9.8_rk
integer :: m,n

write(*,*)
write(*,*)'Test omnidirectional spectrum with Phillips (1958) shape;'
write(*,*)'Vary fetch between 1, 10, and 100 km to compute the peak frequency'&
          //' based on the JONSWAP spectrum.'
write(*,*)

! Initialize spectrum
spec = spectrum_type(fmin=0.0313_rk,fmax=2._rk,df=1.1_rk,ndirs=1,depth=1e3_rk)

! Set array of wind speed and fetch values
wspd = range(5,60,5)
fetch = [1e3_rk,1e4_rk,1e5_rk]

do m = 1,size(fetch)

  write(*,*)
  write(*,fmt='(a,f5.1,a)')'fetch = ',fetch(m)/1e3_rk,' km'
  write(*,fmt='(a)')'   wspd      Hs      Tp       Tm1      Tm2      mss      m0(f)'&
                  //'    m1(f)    m2(f)'
  write(*,fmt='(a)')'----------------------------------------------------------'&
                  //'----------------------'
  do n = 1,size(wspd)
    ! Set the spectrum field to Phillips parametric shape
    call spec % setSpectrum(phillips(spec % getFrequency(),&
      jonswapPeakFrequency(wspd(n),fetch(m),grav),grav))
    write(*,fmt='(9(f8.4,1x))')wspd(n),spec % significantWaveHeight(),&
      1./spec % peakFrequency(),spec % meanPeriod(),&
      spec % meanPeriodZeroCrossing(),spec % meanSquareSlope(),&
      spec % frequencyMoment(0),spec % frequencyMoment(1),&
      spec % frequencyMoment(2)
  enddo

enddo

endprogram test_phillips
