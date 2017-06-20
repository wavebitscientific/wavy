!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.

program test_json

use mod_precision,only:rk => realkind
use mod_testing,only:assert,initialize_tests,report_tests
use mod_utility,only:range
use mod_const,only:WAVY_REAL
use json_module,only:json_core,json_value,json_file

implicit none

integer,parameter :: stdout = 6

logical,dimension(:),allocatable :: tests
logical :: test_failed
integer :: ntests
integer :: n,norder

type(json_core) :: json
type(json_file) :: jfile
type(json_value),pointer :: ptr

logical :: found
real(kind=rk),dimension(:),allocatable :: arr

n = 1

ntests = 1

call initialize_tests(tests,ntests)

write(unit=stdout,fmt='(80("-"))')

call json % initialize()
call json % create_object(ptr,'')
call json % add(ptr,'value',range(1e2_rk,1e3_rk,1e2_rk))
call json % print(ptr,'test.json')
call json % destroy(ptr)

call jfile % initialize()
call jfile % load_file('test.json')
call jfile % get('value', arr, found)
call jfile % destroy()

tests(n) = assert(all(arr == range(1e2_rk,1e3_rk,1e2_rk)),&
                  'write/read array to a json file')
n = n + 1

write(unit=stdout,fmt='(80("-"))')

test_failed = .false.
call report_tests(tests,test_failed)
if(test_failed)stop 1

endprogram test_json
