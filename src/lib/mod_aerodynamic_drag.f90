!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3-clause license. See LICENSE for details.

module mod_aerodynamic_drag

use mod_precision,only:ik => intkind,rk => realkind

implicit none

private

! Public constants
public :: charnockParameter_C1955
public :: charnockParameter_G1977
public :: charnockParameter_W1980

! Public functions
public :: charnockParameter
public :: dragCoefficient_AMV2012
public :: dragCoefficient_DW1962
public :: dragCoefficient_LP1981
public :: dragCoefficient_LY2004
public :: dragCoefficient_TY2001
public :: dragCoefficient_Y1988
public :: dragCoefficient_YT1996
public :: dragCoefficientFromRoughnessLength
public :: dragCoefficientFromFrictionVelocity
public :: skinFrictionVelocity
public :: windAtReferenceHeight

! Named constants
real(kind=rk),parameter :: charnockParameter_C1955 = 0.0120_rk
real(kind=rk),parameter :: charnockParameter_G1977 = 0.0144_rk
real(kind=rk),parameter :: charnockParameter_W1980 = 0.0185_rk

contains


!-------------------------------------------------------------------------------
pure elemental function charnockParameter(z0,ustar,grav) result(alpha)
  !! Returns the Charnock parameter $\alpha$ based on input roughness
  !! length [m], friction velocity [m/s], and gravitational
  !! acceleration [m/s^2].
  !!
  !! TODO add reference
  real(kind=rk),intent(in) :: z0
    !! Roughness length [m]
  real(kind=rk),intent(in) :: ustar
    !! Friction velocity [m/s]
  real(kind=rk),intent(in) :: grav
    !! Gravitational acceleration [m/s^2]
  real(kind=rk) :: alpha
    !! Charnock parameter
  alpha = grav * z0 / ustar**2
endfunction charnockParameter
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function dragCoefficient_AMV2012(wspd) result(cd)
  !! Returns the drag coefficient given the input wind speed at 10 m,
  !! based on Andreas, Maahrt, and Vickers (2012).
  !!
  !! TODO add reference
  real(kind=rk),intent(in) :: wspd
    !! Wind speed at 10 m height [m/s]
  real(kind=rk) :: cd
    !! Drag coefficient
  real(kind=rk) :: ustar
  if(wspd <= 9)then
    ustar = 5.13e-3_rk + 2.83e-2_rk*wspd
  else
    ustar = -2.43e-1_rk + 5.831e-2_rk*wspd
  endif
  cd = dragCoefficientFromFrictionVelocity(ustar,wspd)
endfunction dragCoefficient_AMV2012
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function dragCoefficient_DW1962(wspd) result(cd)
  !! Returns the drag coefficient given the input wind speed at 10 m,
  !! based on Deacon and Webb (1962).
  !!
  !! TODO add reference
  real(kind=rk),intent(in) :: wspd
    !! Wind speed at 10 m height [m/s]
  real(kind=rk) :: cd
    !! Drag coefficient
  cd = (1._rk+0.05_rk*wspd)*1e-3_rk
endfunction dragCoefficient_DW1962
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function dragCoefficient_LP1981(wspd) result(cd)
  !! Returns the drag coefficient given the input wind speed at 10 m,
  !! based on Large and Pond (1981).
  !!
  !! TODO add reference
  real(kind=rk),intent(in) :: wspd
    !! Wind speed at 10 m height [m/s]
  real(kind=rk) :: cd
    !! Drag coefficient
  if(wspd > 11)then
    cd = (0.49_rk+0.065_rk*wspd)*1e-3_rk
  else
    cd = 1.2e-3_rk
  endif
endfunction dragCoefficient_LP1981
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function dragCoefficient_LY2004(wspd) result(cd)
  !! Returns the drag coefficient given the input wind speed at 10 m,
  !! based on Large and Yeager (2004).
  !!
  !! TODO add reference
  real(kind=rk),intent(in) :: wspd
    !! Wind speed at 10 m height [m/s]
  real(kind=rk) :: cd
    !! Drag coefficient
  cd = (0.142_rk + 0.076_rk*wspd + 2.7_rk/wspd)*1e-3_rk
endfunction dragCoefficient_LY2004
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function dragCoefficient_TY2001(significant_height,&
  dominant_wavelength,vonkarman) result(cd)
  !! Returns the drag coefficient given the input significant wave height and
  !! dominant wavelength, based on Taylor and Yelland (2001).
  !!
  !! TODO add reference
  real(kind=rk),intent(in) :: significant_height
    !! Significant wave height [m]
  real(kind=rk),intent(in) :: dominant_wavelength
    !! Dominant wavelength [m]
  real(kind=rk),intent(in) :: vonkarman
    !! Von Karman constant
  real(kind=rk) :: cd
    !! Drag coefficient
  cd = dragCoefficientFromRoughnessLength(1.2e3_rk * significant_height&
     * (significant_height/dominant_wavelength)**4.5_rk,vonkarman)
endfunction dragCoefficient_TY2001
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function dragCoefficient_Y1988(wspd) result(cd)
  !! Returns the drag coefficient given the input wind speed at 10 m,
  !! based on Yelland et al. (1988).
  !!
  !! TODO add reference
  real(kind=rk),intent(in) :: wspd
    !! Wind speed at 10 m height [m/s]
  real(kind=rk) :: cd
    !! Drag coefficient
  cd = (0.5_rk+0.071_rk*wspd)*1e-3_rk
endfunction dragCoefficient_Y1988
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function dragCoefficient_YT1996(wspd) result(cd)
  !! Returns the drag coefficient given the input wind speed at 10 m,
  !! based on Yelland and Taylor (1996).
  !!
  !! TODO add reference
  real(kind=rk),intent(in) :: wspd
    !! Wind speed at 10 m height [m/s]
  real(kind=rk) :: cd
    !! Drag coefficient
  if(wspd < 11)then
    cd = (0.29_rk+3.1_rk/wspd+7.7_rk/wspd**2)*1e-3_rk
  else
    cd = (0.6_rk+0.070_rk*wspd)*1e-3_rk
  endif
endfunction dragCoefficient_YT1996
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function dragCoefficientFromFrictionVelocity(ustar,wspd)&
  result(cd)
  !! Returns the drag coefficient given the input friction velocity and wind
  !! speed.
  real(kind=rk),intent(in) :: ustar
    !! Friction velocity [m/s]
  real(kind=rk),intent(in) :: wspd
    !! Wind speed [m/s]
  real(kind=rk) :: cd
    !! Drag coefficient
  cd = (ustar/wspd)**2
endfunction dragCoefficientFromFrictionVelocity
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental function dragCoefficientFromRoughnessLength(z0,vonkarman)&
  result(cd)
  !! Returns the drag coefficient given the input von Karman constant
  !! and roughness length.
  !!
  !! TODO add reference
  real(kind=rk),intent(in) :: z0
    !! Aerodynamic roughness length [m]
  real(kind=rk),intent(in) :: vonkarman
    !! Von Karman constant
  real(kind=rk) :: cd
    !! Drag coefficient
  cd = (vonkarman/log(10/z0))**2
endfunction dragCoefficientFromRoughnessLength
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk)function skinFrictionVelocity(wspd,height,&
  viscosity,vonkarman) result(ustar)
  !! Computes the friction velocity for a smooth flow over the flat wall,
  !! given input wind speed [m/s], height of input wind [m], air viscosity
  !! [m^2/s], and Von Karman constant.
  real(kind=rk),intent(in) :: wspd
    !! Wind speed at input height above the surface [m/s]
  real(kind=rk),intent(in) :: height
    !! Height of input wind [m]
  real(kind=rk),intent(in) :: viscosity
    !! Viscosity [m^2/s]
  real(kind=rk),intent(in) :: vonkarman
    !! Von Karman constant
  real(kind=rk) :: z0 ! roughness length
  real(kind=rk) :: ustar_old
  integer :: n
  z0 = 0.001_rk
  ustar_old = 1e-2_rk
  do
    ustar = vonkarman*wspd/log(height/z0)
    if(abs(ustar-ustar_old) < 1e-3_rk)then
      return
    endif
    z0 = 0.132_rk*viscosity/ustar
  enddo
endfunction skinFrictionVelocity
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function windAtReferenceHeight(input_wspd,&
  input_height,reference_height,ustar,vonkarman) result(wspd)
  !! Computes the wind speed at reference height, given input wind speed [m/s],
  !! height of input wind [m] at input height, friction velocity [m/s], and
  !! Von Karman constant.
  real(kind=rk),intent(in) :: input_wspd
    !! Wind speed at input height above the surface [m/s]
  real(kind=rk),intent(in) :: input_height
    !! Height of input wind [m]
  real(kind=rk),intent(in) :: reference_height
    !! Reference height [m]
  real(kind=rk),intent(in) :: ustar
    !! Friction velocity [m/s]
  real(kind=rk),intent(in) :: vonkarman
    !! Von Karman constant
  wspd = input_wspd+ustar/vonkarman*log(reference_height/input_height)
endfunction windAtReferenceHeight
!-------------------------------------------------------------------------------
endmodule mod_aerodynamic_drag
