!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.
!
program test_advection

use mod_const
use mod_nondimensional
use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_domain,only:domain_type
use mod_grid,only:grid_type
use mod_spectral_shapes,only:jonswap,jonswapPeakFrequency
use mod_time_integration,only:backward_euler,forward_euler,exact_exponential
use mod_advection
use mod_utility,only:range,ones

implicit none

type(spectrum_type),dimension(:,:),allocatable :: spec
!type(spectrum_type) :: spectrum
type(domain_type) :: domain,tend
!type(grid_type) :: grid

real(kind=rk) :: wspd
real(kind=rk) :: fetch
real(kind=rk),parameter :: grav = 9.8_rk
integer :: i,m,n,idm,nfreqs,ip

real(kind=rk),dimension(:),allocatable :: freq
real(kind=rk),dimension(:,:),allocatable :: f
real(kind=rk),dimension(:,:),allocatable :: cg
real(kind=rk),dimension(:),allocatable :: dx
real(kind=rk),dimension(:,:),allocatable :: adv

idm = 20 

write(*,*)
write(*,*)'Initialize omnidirectional spectrum with JONSWAP shape;'
write(*,*)

wspd = 1e1_rk
fetch = 1e5_rk

! Initialize domain
domain = domain_type(grid = grid_type(1,idm,range(0._rk,2e4_rk,1e3_rk)),&
                     spectrum = spectrum_type(0.04_rk,2._rk,1.1_rk,1,1e3_rk))
spec = domain % getSpectrum()
tend = domain

! Initialize spectrum arrays
!do concurrent(i = 1:idm)
  spec(10,1) = jonswap(spec(1,1) % getFrequency(),wspd,fetch,grav)
!enddo
domain = spec

write(*,fmt='(a)')'   wspd      Hs      Tp       Tm1      Tm2      mss'&
                //'      m0(f)    m1(f)    m2(f)'
write(*,fmt='(a)')'----------------------------------------------------------'&
                //'----------------------'
do n = 1,60

  ip = 10

  !spec = domain % getSpectrum()
  !write(*,fmt='(9(f8.4,1x))')wspd,spec(ip,1) % significantWaveHeight(),&
  !  1./spec(ip,1) % peakFrequency(),spec(ip,1) % meanPeriod(),&
  !  spec(ip,1) % meanPeriodZeroCrossing(),spec(ip,1) % meanSquareSlope(),&
  !  spec(ip,1) % frequencyMoment(0),spec(ip,1) % frequencyMoment(1),&
  !  spec(ip,1) % frequencyMoment(2)

  !domain = forward_euler(domain,&
  !  domain % advect(advectUpwind1stOrder1dRank1,1,&
  !  WAVY_OMNIDIRECTIONAL),dt=30._rk)

  !write(*,fmt='(20(f6.4,1x))')domain % significantWaveHeight()
  write(*,fmt='(20(f7.4,1x))')domain % meanPeriod()

  domain = forward_euler(domain,&
    domain % advect(advectUpwind1stOrder1dRank2,1,&
    WAVY_DIRECTIONAL),dt=60._rk)

enddo
call domain % writeJSON('domain.json',minify=.false.)


endprogram test_advection
