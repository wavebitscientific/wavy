!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.

program test_domain

use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_spectral_shapes,only:jonswap
use mod_domain,only:domain_type
use mod_grid,only:grid_type
use mod_testing,only:assert,initialize_tests,report_tests
use mod_utility,only:range,ones
use mod_const,only:WAVY_INT,WAVY_REAL

implicit none

type(grid_type) :: grid
type(domain_type) :: domain,d1,d2
type(spectrum_type),dimension(:,:),allocatable :: spec
type(spectrum_type) :: spectrum

integer,parameter :: stdout = 6

logical,dimension(:),allocatable :: tests
logical :: test_failed
integer :: i,idm
integer :: ntests
integer :: n,norder

n = 1

ntests = 12

call initialize_tests(tests,ntests)

idm = 101

grid = grid_type(1,idm,range(0._rk,1e4_rk,1e3_rk))
spectrum = spectrum_type(0.03_rk,2._rk,1.1_rk,36,1000._rk)
domain = domain_type(grid,spectrum)

call domain % setDepth(reshape(2e3_rk*ones(idm,WAVY_REAL),[idm,1]))
tests(n) = assert(all(domain % getDepth() == 2e3_rk),&
                  'setter/getter: depth')
n = n + 1

call domain % setElevation(reshape(ones(idm,WAVY_REAL),[idm,1]))
tests(n) = assert(all(domain % getElevation() == 1),&
                  'setter/getter: elevation')
n = n + 1

call domain % setGravity(reshape(9.9_rk*ones(idm,WAVY_REAL),[idm,1]))
tests(n) = assert(all(domain % getGravity() == 9.9_rk),&
                  'setter/getter: gravitational acceleration')
n = n + 1

call domain % setAirDensity(reshape(1.25_rk*ones(idm,WAVY_REAL),[idm,1]))
tests(n) = assert(all(domain % getAirDensity() == 1.25_rk),&
                  'setter/getter: air density')
n = n + 1

call domain % setWaterDensity(reshape(1.03e3_rk*ones(idm,WAVY_REAL),[idm,1]))
tests(n) = assert(all(domain % getWaterDensity() == 1.03e3_rk),&
                  'setter/getter: water density')
n = n + 1

call domain % setSurfaceTension(reshape(0.074_rk*ones(idm,WAVY_REAL),[idm,1]))
tests(n) = assert(all(domain % getSurfaceTension() == 0.074_rk),&
                  'setter/getter: surface tension')
n = n + 1

! Test overloaded operators for domain

d1 = domain_type(grid,spectrum)
d2 = domain_type(grid,spectrum)

spec = d1 % getSpectrum()
! Initialize spectrum arrays
do concurrent(i = 1:idm)
  spec(i,1) = jonswap(spec(i,1) % getFrequency(),1e1_rk,1e5_rk,9.8_rk)
enddo
call d1 % setSpectrum(spec)
call d2 % setSpectrum(spec)

d2 = d1 + 10._rk
tests(n) = assert(d2 == d1 + 10._rk,&
                  'overloaded operators: addition with reals')
n = n + 1

d2 = d1 - 10._rk
tests(n) = assert(d2 == d1 - 10._rk,&
                  'overloaded operators: subtraction with reals')
n = n + 1

d2 = 2._rk*d1
tests(n) = assert(d2 == 2._rk * d1,&
                  'overloaded operators: multiplication with reals')
n = n + 1

d2 = d1 / 2._rk
tests(n) = assert(d2 == d1 / 2._rk,&
                  'overloaded operators: division with reals')
n = n + 1

! Test getter and setter for spectrum array:

domain = domain_type(grid,spectrum)
spec = domain % getSpectrum()
do concurrent(i = 1:idm)
  spec(i,1) = jonswap(spec(i,1) % getFrequency(),1e1_rk,1e5_rk,9.8_rk)
enddo
domain = spec

tests(n) = assert(all(spec == domain % getSpectrum()),&
                  'getter and setter for spectrum attribute')
n = n + 1

d1 = domain
call d2 % setSpectrumArray(d1 % getSpectrumArray([0,0],.false.))
tests(n) = assert(d1 == d2,'getter and setter with spectrum arrays')
n = n + 1

write(unit=stdout,fmt='(80("-"))')

test_failed = .false.
call report_tests(tests,test_failed)
if(test_failed)stop 1

endprogram test_domain
