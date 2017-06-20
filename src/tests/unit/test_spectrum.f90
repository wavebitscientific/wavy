!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.

program test_spectrum

use mod_precision,only:rk => realkind
use mod_spectrum
use mod_spectral_shapes
use mod_testing,only:assert,initialize_tests,report_tests

implicit none

type(spectrum_type) :: s
type(spectrum_type) :: s1,s2,s3

integer,parameter :: stdout = 6

logical,dimension(:),allocatable :: tests
logical :: test_failed
integer :: ntests
integer :: n,norder

n = 1

ntests = 29

call initialize_tests(tests,ntests)

write(unit=stdout,fmt='(80("-"))')

s = spectrum_type(fmin = 0.0313_rk,fmax=2._rk,df=1.1_rk,ndirs=1,depth=1000._rk)

tests(n) = assert(size(s % getFrequency()) == 43,&
                  'expected number of frequencies')
n = n + 1

tests(n) = assert(.not. s % isMonochromatic(),&
                  'spectrum is not monochromatic')
n = n + 1

tests(n) = assert(s % isOmnidirectional(),&
                  'spectrum is omnidirectional')
n = n + 1

tests(n) = assert(all(s % getSpectrum() == 0),&
                  'initial spectrum instance has zero values')
n = n + 1

tests(n) = assert(s % significantWaveHeight() == 0,&
                  'sig. wave height == 0')
n = n + 1

tests(n) = assert(s % meanPeriod() == 0,&
                  'mean wave period == 0')
n = n + 1

tests(n) = assert(s % frequencyMoment(0) == 0,&
                  'frequency moment 0 == 0')
n = n + 1

tests(n) = assert(s % frequencyMoment(1) == 0,&
                  'frequency moment 1 == 0')
n = n + 1

tests(n) = assert(s % frequencyMoment(2) == 0,&
                  'frequency moment 2 == 0')
n = n + 1

tests(n) = assert(s % frequencyMoment(3) == 0,&
                  'frequency moment 3 == 0')
n = n + 1

tests(n) = assert(s % wavenumberMoment(0) == 0,&
                  'wavenumber moment 0 == 0')
n = n + 1

tests(n) = assert(s % wavenumberMoment(1) == 0,&
                  'wavenumber moment 1 == 0')
n = n + 1

tests(n) = assert(s % wavenumberMoment(2) == 0,&
                  'wavenumber moment 2 == 0')
n = n + 1

tests(n) = assert(s % wavenumberMoment(3) == 0,&
                  'wavenumber moment 3 == 0')
n = n + 1

call s % setDepth(2e3_rk)
tests(n) = assert(s % getDepth() == 2e3_rk,&
                  'setter/getter: depth')
n = n + 1

call s % setAirDensity(1.25_rk)
tests(n) = assert(s % getAirDensity() == 1.25_rk,&
                  'setter/getter: air density')
n = n + 1

call s % setWaterDensity(1.03e3_rk)
tests(n) = assert(s % getWaterDensity() == 1.03e3_rk,&
                  'setter/getter: water density')
n = n + 1

call s % setSurfaceTension(0.074_rk)
tests(n) = assert(s % getSurfaceTension() == 0.074_rk,&
                  'setter/getter: surface tension')
n = n + 1

call s % setGravity(9.806_rk)
tests(n) = assert(s % getGravity() == 9.806_rk,&
                  'setter/getter: gravitational acceleration')
n = n + 1

s1 = spectrum_type(fmin = 0.0313_rk,fmax=2._rk,df=1.1_rk,ndirs=1,depth=1000._rk)
s1 = jonswap(s1 % getFrequency(),10._rk,1e5_rk,9.8_rk)
s2 = s1 + 10._rk
tests(n) = assert(s2 == s1 + 10._rk,&
                  'overloaded operators: addition with reals')
n = n + 1

s2 = s1 - 10._rk
tests(n) = assert(s2 == s1 - 10._rk,&
                  'overloaded operators: subtraction with reals')
n = n + 1

s2 = 2._rk*s1
tests(n) = assert(s2 == 2._rk * s1,&
                  'overloaded operators: multiplication with reals')
n = n + 1

s2 = s1 / 2._rk
tests(n) = assert(s2 == s1 / 2._rk,&
                  'overloaded operators: division with reals')
n = n + 1

s2 = s1
tests(n) = assert(s2 == s1,'overloaded operators: equals')
n = n + 1

s2 = s1 + 10._rk
tests(n) = assert(s2 /= s1,'overloaded operators: not equals')
n = n + 1

tests(n) = assert(s2 > s1,'overloaded operators: greater than')
n = n + 1

tests(n) = assert(s2 >= s1,'overloaded operators: greater than or equal')
n = n + 1

tests(n) = assert(s1 < s2,'overloaded operators: less than')
n = n + 1

tests(n) = assert(s1 <= s2,'overloaded operators: less than or equal')
n = n + 1

write(unit=stdout,fmt='(80("-"))')

test_failed = .false.
call report_tests(tests,test_failed)
if(test_failed)stop 1

endprogram test_spectrum
