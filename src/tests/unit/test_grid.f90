!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.

program test_grid

use mod_precision,only:rk => realkind
use mod_grid,only:grid_type
use mod_testing,only:assert,initialize_tests,report_tests
use mod_utility,only:range,ones,tile
use mod_const,only:WAVY_INT

implicit none

type(grid_type) :: grid
real(kind=rk),dimension(:,:),allocatable :: dx,dy
real(kind=rk),dimension(:,:),allocatable :: x,y

integer,parameter :: stdout = 6

logical,dimension(:),allocatable :: tests
logical :: test_failed
integer :: ntests
integer :: n,norder

n = 1

ntests = 8

call initialize_tests(tests,ntests)

write(unit=stdout,fmt='(80("-"))')

! 1-d grid tests

grid = grid_type(1,101,range(0._rk,1e5_rk,1e3_rk))
x = grid % getAxisX()
dx = grid % getGridSpacingX()

tests(n) = assert(all(x(:,1) == range(0._rk,1e5_rk,1e3_rk)),&
                  '1-d grid instance x-axis getter')
n = n + 1

tests(n) = assert(all(dx(:,1) == 1e3_rk*ones(101,WAVY_INT)),&
                  '1-d grid instance x-grid spacing getter')
n = n + 1

tests(n) = assert(all(grid % getLowerBounds() == [1,1]),&
                  '1-d grid instance lower bounds getter')
n = n + 1

tests(n) = assert(all(grid % getUpperBounds() == [101,1]),&
                  '1-d grid instance lower bounds getter')
n = n + 1

tests(n) = assert(size(grid % getAxisY()) == 0,&
                  '1-d grid instance y-axis is size zero')
n = n + 1

tests(n) = assert(size(grid % getGridSpacingY()) == 0,&
                  '1-d grid instance y-spacing is size zero')
n = n + 1

! 2-d grid tests

grid = grid_type([1,1],[101,101],dx = 2e3_rk*tile(ones(101,WAVY_INT),101),&
  dy = 1e3_rk*tile(ones(101,WAVY_INT),101))
x = grid % getAxisX()
y = grid % getAxisY()
dx = grid % getGridSpacingX()
dy = grid % getGridSpacingY()

tests(n) = assert(all(dx == 2e3_rk*tile(ones(101,WAVY_INT),101)),&
                  '2-d grid instance x-grid spacing getter')
n = n + 1

tests(n) = assert(all(dy == 1e3_rk*tile(ones(101,WAVY_INT),101)),&
                  '2-d grid instance y-grid spacing getter')
n = n + 1

write(unit=stdout,fmt='(80("-"))')

test_failed = .false.
call report_tests(tests,test_failed)
if(test_failed)stop 1

endprogram test_grid
