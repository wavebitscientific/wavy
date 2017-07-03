!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.

module mod_const

use mod_precision,only:ik => intkind,rk => realkind

implicit none

private

integer(kind=ik),dimension(1),parameter,public :: WAVY_OMNIDIRECTIONAL = [1]
integer(kind=ik),dimension(1,1),parameter,public :: WAVY_DIRECTIONAL = reshape([1],[1,1])
integer(kind=ik),dimension(1,1,1),parameter,public :: WAVY_DIRECTIONAL_2D = reshape([1],[1,1,1])

real(kind=rk),parameter,public :: WAVY_REAL = 1._rk
integer(kind=ik),parameter,public :: WAVY_INT = 1

integer(kind=ik),parameter,public :: huge_int = huge(1_ik)
real(kind=rk),parameter,public :: tiny_real = tiny(1e0_rk)
real(kind=rk),parameter,public :: huge_real = huge(1e0_rk)

real(kind=rk),parameter,public :: eps = tiny(1e0_rk)
real(kind=rk),parameter,public :: pi = 3.14159265358979323846264338327950_rk
real(kind=rk),parameter,public :: twopi = 2*pi

integer(kind=ik),parameter,public :: stdout = 6
integer(kind=ik),parameter,public :: stderr = 0

endmodule mod_const
