!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.
!
program test_time_integration_domain

use mod_const
use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_domain,only:domain_type
use mod_grid,only:grid_type
use mod_spectral_shapes,only:jonswap
use mod_source_functions,only:sin_DCCM2012,sds_DCCM2012,sdt_DCCM2012
use mod_time_integration,only:forward_euler,exact_exponential
use mod_utility,only:range

implicit none

type(spectrum_type),dimension(:,:),allocatable :: spec
type(domain_type) :: domain,sin_dom,sds_dom,sdt_dom

real(kind=rk) :: wspd
real(kind=rk) :: fetch
real(kind=rk),parameter :: grav = 9.8_rk
integer :: i,m,n,idm

idm = 100

write(*,*)
write(*,*)'Initialize omnidirectional spectrum with JONSWAP shape;'
write(*,*)

wspd = 1e1_rk
fetch = 1e5_rk

! Initialize domain
domain = domain_type(grid = grid_type(1,idm,range(0._rk,1e4_rk,1e3_rk)),&
                     spectrum = spectrum_type(0.0313_rk,2._rk,1.1_rk,1,1e3_rk))
sin_dom = domain
sds_dom = domain
sdt_dom = domain

spec = domain % getSpectrum()

! Initialize spectrum arrays
do concurrent(i = 1:idm)
  spec(i,1) = jonswap(spec(i,1) % getFrequency(),wspd,fetch,grav)
enddo
domain = spec

write(*,fmt='(a)')'   wspd      Hs      Tp       Tm1      Tm2      mss'&
                //'      m0(f)    m1(f)    m2(f)'
write(*,fmt='(a)')'----------------------------------------------------------'&
                //'----------------------'
do n = 1,61

  spec = domain % getSpectrum()
  write(*,fmt='(9(f8.4,1x))')wspd,spec(idm,1) % significantWaveHeight(),&
    1./spec(idm,1) % peakFrequency(),spec(idm,1) % meanPeriod(),&
    spec(idm,1) % meanPeriodZeroCrossing(),spec(idm,1) % meanSquareSlope(),&
    spec(idm,1) % frequencyMoment(0),spec(idm,1) % frequencyMoment(1),&
    spec(idm,1) % frequencyMoment(2)

  sin_dom = sin_DCCM2012(domain % getSpectrum(),wspd,wdir=0._rk,&
                         input_height=10._rk,ustar=1e-2_rk,vonkarman=0.4_rk)

  sds_dom = sds_DCCM2012(domain % getSpectrum(),sds_coefficient=42._rk,&
                         sds_power=2.5_rk,mss_coefficient=120._rk)

  sdt_dom = sdt_DCCM2012(domain % getSpectrum(),1e-2_rk,1e-2_rk)

  !domain = forward_euler(domain,(sin_dom-sds_dom)*domain,dt=1._rk)
  domain = exact_exponential(domain,sin_dom-sds_dom-sdt_dom,dt=60._rk)
  !domain = exact_exponential(domain,tend,dt=60._rk)

enddo

endprogram test_time_integration_domain
