!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.
!
program test_time_integration_spectrum

use mod_const
use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_spectral_shapes,only:jonswap
use mod_source_functions,only:sin_DCCM2012,sds_DCCM2012,sdt_DCCM2012,snl_DCCM2012
use mod_time_integration,only:forward_euler,exact_exponential,backward_euler

implicit none

type(spectrum_type) :: spec,spec2,sin_tend,sds_tend,sdt_tend,snl_tend
real(kind=rk) :: wspd
real(kind=rk) :: fetch
real(kind=rk),parameter :: grav = 9.8_rk
integer :: m,n,stat

write(*,*)
write(*,*)'Initialize omnidirectional spectrum with JONSWAP shape;'
write(*,*)

wspd = 1e1_rk
fetch = 1e5_rk

! Initialize spectrum
spec = spectrum_type(fmin=0.0313_rk,fmax=2._rk,df=1.1_rk,ndirs=1,depth=1e3_rk)
spec = jonswap(spec % getFrequency(),wspd=wspd,fetch=fetch,grav=grav)

spec2 = spectrum_type(fmin=0.0313_rk,fmax=2._rk,df=1.1_rk,ndirs=1,depth=1e3_rk)

call spec % writeJSON('spectrum.json',minify=.false.)
call spec2 % readJSON('spectrum.json')
call spec % writeJSON('spectrum2.json',minify=.false.)

write(*,*)spec == spec2
write(*,*)all(spec % getSpectrum() == spec2 % getSpectrum())

write(*,fmt='(a)')'   wspd      Hs      Tp       Tm1      Tm2      mss'&
                //'      m0(f)    m1(f)    m2(f)'
write(*,fmt='(a)')'----------------------------------------------------------'&
                //'----------------------'
do n = 1,61

  write(*,fmt='(9(f8.4,1x))')wspd,spec % significantWaveHeight(),&
    1./spec % peakFrequency(),spec % meanPeriod(),&
    spec % meanPeriodZeroCrossing(),spec % meanSquareSlope(),&
    spec % frequencyMoment(0),spec % frequencyMoment(1),&
    spec % frequencyMoment(2)

  sin_tend = sin_DCCM2012(spectrum=spec,wspd=wspd,wdir=0._rk,input_height=10._rk,&
                          ustar=1e-2_rk,vonkarman=0.4_rk)
  sds_tend = sds_DCCM2012(spec,42._rk,2.5_rk,120._rk)
  sdt_tend = sdt_DCCM2012(spec,1e-2_rk,1e-2_rk)
  snl_tend = snl_DCCM2012(spec,sds_tend,0.005_rk)

  !spec = forward_euler(spec,(sin_tend-sds_tend-sdt_tend)*spec,dt=60._rk)
  !spec = backward_euler(spec,(sin_tend-sds_tend-sdt_tend)*spec,dt=60._rk)
  spec = exact_exponential(spec,sin_tend-sds_tend-sdt_tend+snl_tend,dt=60._rk)
  !spec = exact_exponential(spec,sin_tend-sds_tend-sdt_tend,dt=60._rk)

enddo

endprogram test_time_integration_spectrum
