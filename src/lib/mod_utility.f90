!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.
!
module mod_utility

use mod_precision,only:ik => intkind,rk => realkind

implicit none

private

public :: diff
public :: diff_periodic
public :: ones
public :: range
public :: tile
public :: zeros

interface diff
  module procedure :: diff_1d
  module procedure :: diff_2d
endinterface diff

interface diff_periodic
  module procedure :: diff_periodic_1d
  module procedure :: diff_periodic_2d
endinterface diff_periodic

interface ones
  module procedure :: ones_int
  module procedure :: ones_real
endinterface ones

interface range
  module procedure :: range_int
  module procedure :: range_real
endinterface range

interface tile
  module procedure :: tile_1d_int
  module procedure :: tile_1d_real
  module procedure :: tile_2d_int
  module procedure :: tile_2d_real
  module procedure :: tile_3d_int
  module procedure :: tile_3d_real
endinterface tile

interface zeros
  module procedure :: zeros_int
  module procedure :: zeros_real
endinterface zeros

!===============================================================================
contains



!-------------------------------------------------------------------------------
pure function diff_1d(x) result(dx)
  !! Returns a centered-difference of a 1-d array, with first order
  !! differencing applied for the boundary points. This procedure is overloaded
  !! by the generic procedure `diff`.
  real(kind=rk),dimension(:),intent(in)  :: x !! Input array
  real(kind=rk),dimension(:),allocatable :: dx
  integer(kind=ik) :: idm
  idm = size(x)
  allocate(dx(idm))
  if(idm == 0)then
    return
  elseif(idm == 1)then
    dx = 0
    return
  endif
  dx(2:idm-1) = 0.5_rk*(x(3:idm)-x(1:idm-2))
  dx(1) = x(2)-x(1)
  dx(idm) = x(idm)-x(idm-1)
endfunction diff_1d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function diff_2d(x,dim) result(dx)
  !! Returns a centered-difference of a 2-d array along dimension dim, with
  !! first order differencing applied for the boundary points. This procedure is
  !! overloaded by the generic procedure `diff`.
  real(kind=rk),dimension(:,:),intent(in) :: x
    !! Input array
  integer(kind=ik),intent(in) :: dim
    !! Dimension along which to differentiate
  real(kind=rk),dimension(:,:),allocatable :: dx
  integer(kind=ik) :: idm,jdm
  idm = size(x,dim=1)
  jdm = size(x,dim=2)
  allocate(dx(idm,jdm))
  if(idm == 0)then
    return
  elseif(idm == 1)then
    dx = 0
    return
  endif
  if(dim == 1)then
    dx(2:idm-1,:) = 0.5_rk*(x(3:idm,:)-x(1:idm-2,:))
    dx(1,:) = x(2,:)-x(1,:)
    dx(idm,:) = x(idm,:)-x(idm-1,:)
  elseif(dim == 2)then
    dx(:,2:idm-1) = 0.5_rk*(x(:,3:idm)-x(:,1:idm-2))
    dx(:,1) = x(2,:)-x(:,1)
    dx(:,idm) = x(:,idm)-x(:,idm-1)
  else
    dx = 0
  endif
endfunction diff_2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function diff_periodic_1d(x) result(dx)
  !! Returns a centered-difference of a 1-d array with periodic boundary
  !! conditions. This procedure is overloaded by the generic procedure `diff`.
  real(kind=rk),dimension(:),intent(in)  :: x !! Input array
  real(kind=rk),dimension(:),allocatable :: dx
  integer(kind=ik) :: idm
  idm = size(x)
  allocate(dx(idm))
  if(idm == 0)then
    return
  elseif(idm == 1)then
    dx = 0
    return
  endif
  dx(2:idm-1) = 0.5_rk*(x(3:idm)-x(1:idm-2))
  dx(1) = 0.5_rk*(x(2)-x(idm))
  dx(idm) = 0.5_rk*(x(1)-x(idm-1))
endfunction diff_periodic_1d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function diff_periodic_2d(x,dim) result(dx)
  !! Returns a centered-difference of a 2-d array along dimension dim, with
  !! periodic boundary conditions. This procedure is overloaded by the generic
  !! procedure `diff`.
  real(kind=rk),dimension(:,:),intent(in) :: x
    !! Input array
  integer(kind=ik),intent(in) :: dim
    !! Dimension along which to differentiate
  real(kind=rk),dimension(:,:),allocatable :: dx
  integer(kind=ik) :: idm,jdm
  idm = size(x,dim=1)
  jdm = size(x,dim=2)
  allocate(dx(idm,jdm))
  if(idm == 0)then
    return
  elseif(idm == 1)then
    dx = 0
    return
  endif
  if(dim == 1)then
    dx(2:idm-1,:) = 0.5_rk*(x(3:idm,:)-x(1:idm-2,:))
    dx(1,:) = 0.5_rk*(x(2,:)-x(idm,:))
    dx(idm,:) = 0.5_rk*(x(1,:)-x(idm-1,:))
  elseif(dim == 2)then
    dx(:,2:idm-1) = 0.5_rk*(x(:,3:idm)-x(:,1:idm-2))
    dx(:,1) = 0.5_rk*(x(:,2)-x(:,idm))
    dx(:,idm) = 0.5_rk*(x(:,1)-x(:,idm-1))
  else
    dx = 0
  endif
endfunction diff_periodic_2d
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function ones_int(length,kindflag) result(ones)
  !! Returns a 1-d array of integer ones. This procedure is overloaded by the
  !! generic procedure `ones`.
  integer(kind=ik),intent(in) :: length   !! Array length
  integer(kind=ik),intent(in) :: kindflag !! Array type
  integer(kind=ik),dimension(:),allocatable :: ones
  allocate(ones(length))
  ones = 1
endfunction ones_int
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function ones_real(length,kindflag) result(ones)
  !! Returns a 1-d array of floating-point ones. This procedure is overloaded by
  !! the generic procedure `ones`.
  integer(kind=ik),intent(in) :: length !! Array length
  real(kind=rk),intent(in) :: kindflag !! Array type
  real(kind=rk),dimension(:),allocatable :: ones
  allocate(ones(length))
  ones = 1._rk
endfunction ones_real
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function range_int(start,end,increment) result(range)
  !! Returns an array of integers given start, end, and increment values. If the
  !! increment argument is not passed, default increment is 1. This procedure is
  !! overloaded by the generic procedure `range`.
  integer(kind=ik),intent(in) :: start
    !! Start value of the array
  integer(kind=ik),intent(in) :: end
    !! End value of the array
  integer(kind=ik),intent(in),optional :: increment
    !! Array increment
  integer(kind=ik),dimension(:),allocatable :: range
  integer(kind=ik) :: i
  integer(kind=ik) :: increment_
  integer(kind=ik) :: length
  if(present(increment))then
    increment_ = increment
  else
    increment_ = 1
  endif
  length = (end-start)/increment_+1
  allocate(range(length))
  do concurrent(i = 1:length)
    range(i) = start+(i-1)*increment_
  enddo
endfunction range_int
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function range_real(start,end,increment) result(range)
  !! Returns an array of reals given start, end, and increment values. If the
  !! increment argument is not passed, default increment is 1. This procedure is
  !! overloaded by the generic procedure `range`.
  real(kind=rk),intent(in) :: start
    !! Start value of the array
  real(kind=rk),intent(in) :: end
    !! End value of the array
  real(kind=rk),intent(in),optional :: increment
    !! Array increment
  real(kind=rk),dimension(:),allocatable :: range
  real(kind=rk) :: increment_
  integer(kind=ik) :: i
  integer(kind=ik) :: length
  if(present(increment))then
    increment_ = increment
  else
    increment_ = 1
  endif
  length = int((end-start)/increment_)+1
  allocate(range(length))
  do concurrent(i = 1:length)
    range(i) = start+(i-1)*increment_
  enddo
endfunction range_real
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function tile_1d_int(array,n) result(tiled_array)
  !! Tiles the input array `n` times. Returns a tiled array that has rank equal
  !! to `size(shape(array))+1` and that has values equal to values of `array`,
  !! repeated `n` times. This version is for 1-d input array of integers. This
  !! procedure is overloaded by the generic procedure `tile`.
  integer(kind=ik),dimension(:),intent(in) :: array !! Input array
  integer(kind=ik),intent(in) :: n !! Number of times to copy input array
  integer(kind=ik),dimension(:,:),allocatable :: tiled_array
  integer(kind=ik) :: i
  allocate(tiled_array(size(array),n))
  do concurrent(i=1:n)
    tiled_array(:,i) = array(:)
  enddo
endfunction tile_1d_int
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function tile_1d_real(array,n) result(tiled_array)
  !! Tiles the input array `n` times. Returns a tiled array that has rank equal
  !! to `size(shape(array))+1` and that has values equal to values of `array`,
  !! repeated `n` times. This version is for 1-d input array of reals. This
  !! procedure is overloaded by the generic procedure `tile`.
  real(kind=rk),dimension(:),intent(in) :: array !! Input array
  integer(kind=ik),intent(in) :: n !! Number of times to copy input array
  real(kind=rk),dimension(:,:),allocatable :: tiled_array
  integer(kind=ik) :: i
  allocate(tiled_array(size(array),n))
  do concurrent(i=1:n)
    tiled_array(:,i) = array(:)
  enddo
endfunction tile_1d_real
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function tile_2d_int(array,n) result(tiled_array)
  !! Tiles the input array `n` times. Returns a tiled array that has rank equal
  !! to `size(shape(array))+1` and that has values equal to values of `array`,
  !! repeated `n` times. This version is for 2-d input array of integers. This
  !! procedure is overloaded by the generic procedure `tile`.
  integer(kind=ik),dimension(:,:),intent(in) :: array !! Input array
  integer(kind=ik),intent(in) :: n !! Number of times to copy input array
  integer(kind=ik),dimension(:,:,:),allocatable :: tiled_array
  integer(kind=ik) :: i
  allocate(tiled_array(size(array,dim=1),size(array,dim=2),n))
  do concurrent(i=1:n)
    tiled_array(:,:,i) = array(:,:)
  enddo
endfunction tile_2d_int
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function tile_2d_real(array,n) result(tiled_array)
  !! Tiles the input array `n` times. Returns a tiled array that has rank equal
  !! to `size(shape(array))+1` and that has values equal to values of `array`,
  !! repeated `n` times. This version is for 2-d input array of reals. This
  !! procedure is overloaded by the generic procedure `tile`.
  real(kind=rk),dimension(:,:),intent(in) :: array !! Input array
  integer(kind=ik),intent(in) :: n !! Number of times to copy input array
  real(kind=rk),dimension(:,:,:),allocatable :: tiled_array
  integer(kind=ik) :: i
  allocate(tiled_array(size(array,dim=1),size(array,dim=2),n))
  do concurrent(i=1:n)
    tiled_array(:,:,i) = array(:,:)
  enddo
endfunction tile_2d_real
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function tile_3d_int(array,n) result(tiled_array)
  !! Tiles the input array `n` times. Returns a tiled array that has rank equal
  !! to `size(shape(array))+1` and that has values equal to values of `array`,
  !! repeated `n` times. This version is for 3-d input array of integers. This
  !! procedure is overloaded by the generic procedure `tile`.
  integer(kind=ik),dimension(:,:,:),intent(in) :: array !! Input array
  integer(kind=ik),intent(in) :: n !! Number of times to copy input array
  integer(kind=ik),dimension(:,:,:,:),allocatable :: tiled_array
  integer(kind=ik) :: i
  allocate(tiled_array(size(array,dim=1),size(array,dim=2),size(array,dim=3),n))
  do concurrent(i=1:n)
    tiled_array(:,:,:,i) = array(:,:,:)
  enddo
endfunction tile_3d_int
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function tile_3d_real(array,n) result(tiled_array)
  !! Tiles the input array `n` times. Returns a tiled array that has rank equal
  !! to `size(shape(array))+1` and that has values equal to values of `array`,
  !! repeated `n` times. This version is for 3-d input array of reals. This
  !! procedure is overloaded by the generic procedure `tile`.
  real(kind=rk),dimension(:,:,:),intent(in) :: array !! Input array
  integer(kind=ik),intent(in) :: n !! Number of times to copy input array
  real(kind=rk),dimension(:,:,:,:),allocatable :: tiled_array
  integer(kind=ik) :: i
  allocate(tiled_array(size(array,dim=1),size(array,dim=2),size(array,dim=3),n))
  do concurrent(i=1:n)
    tiled_array(:,:,:,i) = array(:,:,:)
  enddo
endfunction tile_3d_real
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function zeros_int(length,kindflag) result(zeros)
  !! Returns a 1-d array of integer zeros. This procedure is overloaded by the
  !! generic procedure `zeros`.
  integer(kind=ik),intent(in) :: length !! Array length
  integer(kind=ik),intent(in) :: kindflag !! Array type
  integer(kind=ik),dimension(:),allocatable :: zeros
  allocate(zeros(length))
  zeros = 0
endfunction zeros_int
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function zeros_real(length,kindflag) result(zeros)
  !! Returns a 1-d array of floating-point zeros. This procedure is overloaded by
  !! the generic procedure `zeros`.
  integer(kind=ik),intent(in) :: length !! Array length
  real(kind=rk),intent(in) :: kindflag !! Array type
  real(kind=rk),dimension(:),allocatable :: zeros
  allocate(zeros(length))
  zeros = 0._rk
endfunction zeros_real
!-------------------------------------------------------------------------------
endmodule mod_utility
