!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.

program test_utility

use mod_testing,only:assert,initialize_tests,report_tests
use mod_precision,only:rk => realkind,ik => intkind
use mod_utility,only:diff,diff_periodic,zeros,ones,range,tile
use mod_const,only:WAVY_REAL,WAVY_INT

implicit none

integer,parameter :: ndm = 1000

integer,parameter :: stdout = 6

logical,dimension(:),allocatable :: tests
logical :: test_failed
integer :: ntests
integer :: n,norder

real(kind=rk),dimension(10) :: ones1d
real(kind=rk),dimension(10,10) :: ones2d
real(kind=rk),dimension(10,10,10) :: ones3d
real(kind=rk),dimension(10,10,10,10) :: ones4d

integer(kind=ik),dimension(10) :: ones1d_int
integer(kind=ik),dimension(10,10) :: ones2d_int
integer(kind=ik),dimension(10,10,10) :: ones3d_int
integer(kind=ik),dimension(10,10,10,10) :: ones4d_int

n = 1

! Edit this number when adding new tests
ntests = 41

call initialize_tests(tests,ntests)

tests(n) = assert(all(zeros(1,WAVY_INT) == [0]),&
                  'all(zeros(1,WAVY_INT) %  == [0])')
n = n + 1

tests(n) = assert(all(zeros(5,WAVY_INT) == [0,0,0,0,0]),&
                  'all(zeros(5,WAVY_INT) == [0,0,0,0,0]')
n = n + 1

tests(n) = assert(size(zeros(0,WAVY_INT)) == 0,&
                  'size(zeros(0,WAVY_INT)) == 0')
n = n + 1

tests(n) = assert(size(zeros(-5,WAVY_INT)) == 0,&
                  'size(zeros(-5,WAVY_INT)) == 0')
n = n + 1

tests(n) = assert(all(zeros(1,WAVY_REAL) == [0]),&
                  'all(zeros(1,WAVY_REAL) %  == [0])')
n = n + 1

tests(n) = assert(all(zeros(5,WAVY_REAL) == [0,0,0,0,0]),&
                  'all(zeros(5,WAVY_REAL) == [0,0,0,0,0]')
n = n + 1

tests(n) = assert(size(zeros(0,WAVY_REAL)) == 0,&
                  'size(zeros(0,WAVY_REAL)) == 0')
n = n + 1

tests(n) = assert(size(zeros(-5,WAVY_REAL)) == 0,&
                  'size(zeros(-5,WAVY_REAL)) == 0')
n = n + 1

tests(n) = assert(all(zeros(5,WAVY_REAL) == zeros(5,WAVY_INT)),&
                  'all(zeros(5,WAVY_REAL) == zeros(5,WAVY_INT))')
n = n + 1

tests(n) = assert(all(ones(5,WAVY_REAL) == ones(5,WAVY_INT)),&
                  'all(ones(5,WAVY_REAL) == ones(5,WAVY_INT))')
n = n + 1

tests(n) = assert(all(range(1,5,1) == [1,2,3,4,5]),&
                  'all(range(1,5,1) == [1,2,3,4,5])')
n = n + 1

tests(n) = assert(all(range(1,5) == [1,2,3,4,5]),&
                  'all(range(1,5) == [1,2,3,4,5])')
n = n + 1

tests(n) = assert(all(range(1,5) == range(1,5,1)),&
                  'all(range(1,5) == range(1,5,1))')
n = n + 1

tests(n) = assert(all(range(1,5,2) == [1,3,5]),&
                  'all(range(1,5,2) == [1,3,5])')
n = n + 1

tests(n) = assert(all(range(5,1,-1) == [5,4,3,2,1]),&
                  'all(range(5,1,-1) == [5,4,3,2,1])')
n = n + 1

tests(n) = assert(all(range(1._rk,5._rk,1._rk) == [1,2,3,4,5]),&
                  'all(range(1._rk,5._rk,1._rk) == [1,2,3,4,5])')
n = n + 1

tests(n) = assert(all(range(1._rk,5._rk) == [1,2,3,4,5]),&
                  'all(range(1._rk,5._rk) == [1,2,3,4,5])')
n = n + 1

tests(n) = assert(all(range(1._rk,5._rk) == range(1._rk,5._rk,1._rk)),&
                  'all(range(1._rk,5._rk) == range(1._rk,5._rk,1._rk))')
n = n + 1

tests(n) = assert(all(range(1._rk,5._rk,2._rk) == [1._rk,3._rk,5._rk]),&
                  'all(range(1._rk,5._rk,2._rk) == [1._rk,3._rk,5._rk])')
n = n + 1

tests(n) = assert(all(range(5._rk,1._rk,-1._rk)             &
                        == [5._rk,4._rk,3._rk,2._rk,1._rk]),&
                  'all(range(5._rk,1._rk,-1._rk)'           &
                   //'== [5._rk,4._rk,3._rk,2._rk,1._rk])')
n = n + 1

tests(n) = assert(all(range(1._rk,5._rk) == range(1,5)),&
                  'all(range(1._rk,5._rk) == range(1,5))')
n = n + 1

tests(n) = assert(all(diff(zeros(10,WAVY_REAL)) == zeros(10,WAVY_REAL)),&
                  'all(diff(zeros(10,WAVY_REAL)) == zeros(10,WAVY_REAL))')
n = n + 1

tests(n) = assert(all(diff(ones(10,WAVY_REAL)) == zeros(10,WAVY_REAL)),&
                  'all(diff(ones(10,WAVY_REAL)) == zeros(10,WAVY_REAL))')
n = n + 1

tests(n) = assert(all(diff(zeros(1,WAVY_REAL)) == zeros(1,WAVY_REAL)),&
                  'all(diff(zeros(1,WAVY_REAL)) == zeros(1,WAVY_REAL))')
n = n + 1

tests(n) = assert(all(diff(ones(1,WAVY_REAL)) == zeros(1,WAVY_REAL)),&
                  'all(diff(ones(1,WAVY_REAL)) == zeros(1,WAVY_REAL))')
n = n + 1

tests(n) = assert(all(diff(ones(0,WAVY_REAL)) == zeros(0,WAVY_REAL)),&
                  'all(diff(ones(0,WAVY_REAL)) == zeros(0,WAVY_REAL))')
n = n + 1

tests(n) = assert(all(diff_periodic(zeros(10,WAVY_REAL)) == zeros(10,WAVY_REAL)),&
                  'all(diff_periodic(zeros(10,WAVY_REAL)) == zeros(10,WAVY_REAL))')
n = n + 1

tests(n) = assert(all(diff_periodic(ones(10,WAVY_REAL)) == zeros(10,WAVY_REAL)),&
                  'all(diff_periodic(ones(10,WAVY_REAL)) == zeros(10,WAVY_REAL))')
n = n + 1

tests(n) = assert(all(diff_periodic(zeros(1,WAVY_REAL)) == zeros(1,WAVY_REAL)),&
                  'all(diff_periodic(zeros(1,WAVY_REAL)) == zeros(1,WAVY_REAL))')
n = n + 1

tests(n) = assert(all(diff_periodic(ones(1,WAVY_REAL)) == zeros(1,WAVY_REAL)),&
                  'all(diff_periodic(ones(1,WAVY_REAL)) == zeros(1,WAVY_REAL))')
n = n + 1

tests(n) = assert(all(diff_periodic(ones(0,WAVY_REAL)) == zeros(0,WAVY_REAL)),&
                  'all(diff_periodic(ones(0,WAVY_REAL)) == zeros(0,WAVY_REAL))')
n = n + 1

ones1d = 1
ones2d = 1
ones3d = 1
ones4d = 1

tests(n) = assert(all(tile(ones1d,10) == ones2d),'tile_1d_realkind')
n = n + 1

tests(n) = assert(all(tile(ones2d,10) == ones3d),'tile_2d_realkind')
n = n + 1

tests(n) = assert(all(tile(ones3d,10) == ones4d),'tile_3d_realkind')
n = n + 1

tests(n) = assert(size(tile(ones1d,0)) == 0,&
                  'tile_1d_realkind result of zero copies size')
n = n + 1

tests(n) = assert(all(shape(tile(ones1d,0)) == [10,0]),&
                  'tile_1d_realkind result of zero copies shape')
n = n + 1

ones1d_int = 1
ones2d_int = 1
ones3d_int = 1
ones4d_int = 1

tests(n) = assert(all(tile(ones1d_int,10) == ones2d_int),'tile_1d_intkind')
n = n + 1

tests(n) = assert(all(tile(ones2d_int,10) == ones3d_int),'tile_2d_intkind')
n = n + 1

tests(n) = assert(all(tile(ones3d_int,10) == ones4d_int),'tile_3d_intkind')
n = n + 1

tests(n) = assert(size(tile(ones1d_int,0)) == 0,&
                  'tile_1d_intkind result of zero copies size')
n = n + 1

tests(n) = assert(all(shape(tile(ones1d_int,0)) == [10,0]),&
                  'tile_1d_intkind result of zero copies shape')
n = n + 1

write(unit=stdout,fmt='(80("-"))')

test_failed = .false.
call report_tests(tests,test_failed)
if(test_failed)stop 1

endprogram test_utility
