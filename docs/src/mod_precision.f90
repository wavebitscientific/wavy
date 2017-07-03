!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.

module mod_precision

use iso_fortran_env,only:int16,int32,int64,real32,real64,real128

implicit none

private

public :: intkind,realkind

#ifdef REAL64
integer,parameter :: realkind = real64
#elif REAL128
integer,parameter :: realkind = real128
#else
integer,parameter :: realkind = real32
#endif

integer,parameter :: intkind = int32

endmodule mod_precision
