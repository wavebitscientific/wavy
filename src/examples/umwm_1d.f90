!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2016, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.
!
program example_umwm_1d

use mod_const
use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_domain,only:domain_type
use mod_grid,only:grid_type
use mod_spectral_shapes,only:jonswap
use mod_source_functions,only:sin_DCCM2012,sds_DCCM2012,sdt_DCCM2012,snl_DCCM2012
use mod_time_integration,only:forward_euler,exact_exponential
use mod_functional,only:arange
use mod_advection

implicit none

type(grid_type) :: grid
type(spectrum_type) :: spectrum
type(spectrum_type),dimension(:,:),allocatable :: spec
type(domain_type) :: domain
type(domain_type) :: sin_dom,sds_dom,sdt_dom,sbf_dom,snl_dom

real(kind=rk) :: wspd,fetch,grav
integer :: i,m,n,idm

write(*,*)
write(*,*)'Initialize omnidirectional spectrum with JONSWAP shape;'
write(*,*)

idm = size(arange(0._rk,1e4_rk,1e3_rk))

! initialize 1-d grid, from 0 to 100 km with 1 km grid size
grid = grid_type(1,idm,arange(0._rk,1e4_rk,1e3_rk))

! initialize omnidirectional spectrum
spectrum = spectrum_type(0.0313_rk,2._rk,1.1_rk,1,1e3_rk)

! initialize domain
domain = domain_type(grid,spectrum)

! initialize tendency fields
sin_dom = domain
sds_dom = domain
snl_dom = domain
sdt_dom = domain
sbf_dom = domain

! initialize spectra to JONSWAP shape everywhere
wspd = 1e1_rk
fetch = 1e5_rk
grav = 9.8_rk
spec = domain % getSpectrum()
!do concurrent(i = 1:idm)
do concurrent(i = 1:1)
  spec(i,1) = jonswap(spec(i,1) % getFrequency(),wspd,fetch,grav)
enddo
!do concurrent(i = 2:idm)
!  spec(i,1) = jonswap(spec(i,1) % getFrequency(),1._rk,fetch,grav)
!enddo
domain = spec

write(*,fmt='(a)')'   wspd      Hs      Tp       Tm1      Tm2      mss'&
                //'      m0(f)    m1(f)    m2(f)'
write(*,fmt='(a)')'----------------------------------------------------------'&
                //'----------------------'
do n = 1,1441

  spec = domain % getSpectrum()
  !write(*,fmt='(9(f8.4,1x))')(n-1)*60./86400,spec(idm,1) % significantWaveHeight(),&
  !  1./spec(idm,1) % peakFrequency(),spec(idm,1) % meanPeriod(),&
  !  spec(idm,1) % meanPeriodZeroCrossing(),spec(idm,1) % meanSquareSlope(),&
  !  spec(idm,1) % frequencyMoment(0),spec(idm,1) % frequencyMoment(1),&
  !  spec(idm,1) % frequencyMoment(2)

  ! wind input
  sin_dom = sin_DCCM2012(domain % getSpectrum(),wspd,wdir=0._rk,&
                         input_height=10._rk,ustar=1e-2_rk,vonkarman=0.4_rk)

  ! wave dissipation
  sds_dom = sds_DCCM2012(domain % getSpectrum(),sds_coefficient=42._rk,&
                         sds_power=2.5_rk,mss_coefficient=120._rk)

  ! dissipation from turbulence
  sdt_dom = sdt_DCCM2012(domain % getSpectrum(),1e-2_rk,1e-2_rk)

  ! non-linear wave-wave transfer
  snl_dom = snl_DCCM2012(domain % getSpectrum(),sds_dom % getSpectrum(),5._rk)

  !domain = forward_euler(domain,(sin_dom-sds_dom-sdt_dom)*domain+snl_dom,dt=1._rk)
  domain = exact_exponential(domain,sin_dom-sds_dom-sdt_dom,dt=10._rk)
  !domain = forward_euler(domain,snl_dom,dt=10._rk)

  domain = forward_euler(domain,&
      domain % advect(advectUpwind1stOrder1dRank1,1,WAVY_OMNIDIRECTIONAL),&
      dt=10._rk)

  write(*,fmt='(11(f7.4,1x))')domain % significantWaveHeight()
  !write(*,fmt='(40(f7.4,1x))')domain % meanPeriod()

enddo

endprogram example_umwm_1d
