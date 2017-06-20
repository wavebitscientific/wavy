!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.
!
program test_advection_2d

use mod_const
use mod_nondimensional
use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_domain,only:domain_type
use mod_grid,only:grid_type
use mod_spectral_shapes
use mod_time_integration,only:backward_euler,forward_euler,exact_exponential
use mod_advection
use mod_utility,only:tile,range,ones

implicit none

type(spectrum_type),dimension(:,:),allocatable :: spec
!type(spectrum_type) :: spectrum
type(domain_type) :: domain,tend
!type(grid_type) :: grid

real :: t0,t1

real(kind=rk) :: wspd
real(kind=rk) :: fetch
real(kind=rk),parameter :: grav = 9.8_rk
integer :: i,m,n,idm,nfreqs,jp,ip

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
domain = domain_type(grid = grid_type([1,1],[50,20],&
  dx = 1e3_rk*tile(ones(50,WAVY_INT),20),&
  dy = 1e3_rk*tile(ones(50,WAVY_INT),20)),&
  spectrum = spectrum_type(0.0313_rk,2._rk,1.1_rk,36,1e3_rk))

spec = domain % getSpectrum()
tend = domain

! Initialize spectrum arrays
!do concurrent(i = 1:idm)
  spec(5,5) = donelanHamiltonHuiDirectionalSpectrum(f = domain % getFrequency(),&
                              theta = domain % getDirections(),&
                              wspd  = wspd,&
                              fpeak = jonswapPeakFrequency(wspd,fetch,grav),&
                              theta_mean = 0._rk,&
                              grav  = grav)

!enddo
domain = spec

call cpu_time(t0)
call domain % writeJSON('domain.json',minify=.true.)
call cpu_time(t1)
write(*,*)'domain.json elapsed',t1-t0,'seconds'

write(*,fmt='(a)')'   wspd      Hs      Tp       Tm1      Tm2      mss'&
                //'      m0(f)    m1(f)    m2(f)'
write(*,fmt='(a)')'----------------------------------------------------------'&
                //'----------------------'

do n = 1,60

  ip = 1 
  jp = 5 

  spec = domain % getSpectrum()
  write(*,fmt='(9(f8.4,1x))')wspd,spec(ip,jp) % significantWaveHeight(),&
    1./spec(ip,jp) % peakFrequency(),spec(ip,jp) % meanPeriod(),&
    spec(ip,jp) % meanPeriodZeroCrossing(),spec(ip,jp) % meanSquareSlope(),&
    spec(ip,jp) % frequencyMoment(0),spec(ip,jp) % frequencyMoment(1),&
    spec(ip,jp) % frequencyMoment(2)

  !domain = forward_euler(domain,&
  !  domain % advect(advectUpwind1stOrder1dRank1,1,&
  !  WAVY_OMNIDIRECTIONAL),dt=30._rk)

  !write(*,fmt='(20(f6.4,1x))')domain % significantWaveHeight()
  !write(*,fmt='(20(f7.4,1x))')domain % meanPeriod()

  domain = forward_euler(domain,&
    domain % advect(advectUpwind1stOrder2dRank2,[1,1]),dt=60._rk)

enddo
!call domain % writeJSON('domain.json',minify=.false.)


endprogram test_advection_2d
