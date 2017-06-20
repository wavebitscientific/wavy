!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.

module mod_spectral_shapes

use mod_precision,only:ik => intkind,rk => realkind
use mod_const,only:twopi
use mod_nondimensional,only:nondimensionalFetch,nondimensionalFrequency

implicit none

private

public :: donelanHamiltonHui
public :: donelanHamiltonHuiDirectionalSpreading
public :: donelanHamiltonHuiDirectionalSpectrum
public :: jonswap
public :: jonswapPeakFrequency
public :: piersonMoskowitz
public :: piersonMoskowitzPeakFrequency
public :: phillips

!===============================================================================
contains



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function donelanHamiltonHui(f,fpeak,wspd,grav)&
  result(spec)

  !! The omnidirectional spectrum function based on the laboratory and field
  !! measurements by Donelan, Hamilton, and Hui (1985).
  !!
  !! References:
  !!
  !! Donelan, M. A., J. Hamilton, and W. H. Hui, 1985. Directional
  !! spectra of wind-generated waves. Phil. Trans. Royal Soc. London A.,
  !! 315, 509-562.

  real(kind=rk),intent(in) :: f !! Frequency [Hz]
  real(kind=rk),intent(in) :: fpeak !! Peak frequency [Hz]
  real(kind=rk),intent(in) :: wspd !! Wind speed at 10 m height [m/s]
  real(kind=rk),intent(in) :: grav !! Gravitational acceleration [m/s^2]

  real(kind=rk) :: beta
  real(kind=rk) :: gamma ! Peak enhancement factor
  real(kind=rk) :: nu ! Nondimensional frequency
  real(kind=rk) :: r ! Exponent in peak enhancement factor
  real(kind=rk) :: omega ! Radian frequency [rad/s]
  real(kind=rk) :: sigma ! Spread of the spectrum peak

  omega = twopi*f
  nu = nondimensionalFrequency(wspd,fpeak,grav)
  if(nu >= 0.159_rk)then
    gamma = 6.489_rk + 6*log10(nu)
    if(gamma < 0)gamma = 0
  else
    gamma = 1.7_rk
  endif
  sigma = 0.08_rk + 1.29e-3_rk*nu**(-3)
  r = exp(-0.5_rk*((f-fpeak)/(sigma*fpeak))**2)
  beta = 0.0165_rk*nu**0.55_rk
  spec = twopi*beta*grav**2/omega**4/fpeak*exp(-(fpeak/f)**4)*gamma**r

endfunction donelanHamiltonHui
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function donelanHamiltonHuiDirectionalSpreading(f,&
  wspd,fpeak,theta,theta_mean) result(spreading)

  !! Directional spreading function based on the laboratory and field
  !! measurements by Donelan, Hamilton, and Hui (1985). Includes the
  !! high-frequency form for beta_s found by Banner (1990).
  !!
  !! References:
  !!
  !! Donelan, M. A., J. Hamilton, and W. H. Hui, 1985. Directional
  !! spectra of wind-generated waves. *Phil. Trans. Royal Soc. London A.*,
  !! **315**, 509-562.
  !!
  !! Banner, M. L., 1990. Equilibrium spectra of wind waves. *J. Phys. 
  !! Oceanogr.*, **20**, 966-984.

  real(kind=rk),intent(in) :: f !! Frequency [Hz]
  real(kind=rk),intent(in) :: wspd !! Wind speed at 10 m height [m/s]
  real(kind=rk),intent(in) :: fpeak !! Peak frequency [Hz]
  real(kind=rk),intent(in) :: theta !! Wave direction [rad]
  real(kind=rk),intent(in) :: theta_mean !! Mean wave direction [rad]

  real(kind=rk) :: frel ! Frequency relative to peak frequency, f/fpeak
  real(kind=rk) :: nu ! Nondimensional frequency
  real(kind=rk) :: beta

  frel = f/fpeak

  if(frel < 0.56_rk)then
    beta = 1.24_rk
  elseif(frel >= 0.56_rk .and. frel < 0.95_rk)then
    beta = 2.61_rk*frel**1.3_rk
  elseif(frel >= 0.95_rk .and. frel < 1.6_rk)then
    beta = 2.28_rk*frel**(-1.3_rk)
  else
    beta = 10**(-0.4_rk+0.8393_rk*exp(-0.567_rk*log(frel**2)))
  endif

  spreading = 0.5_rk*beta/cosh(beta*(theta-theta_mean))**2

endfunction donelanHamiltonHuiDirectionalSpreading
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure function donelanHamiltonHuiDirectionalSpectrum(f,theta,wspd,fpeak,&
  theta_mean,grav) result(spec)

  !! Returns directional frequency spectrum based on the laboratory and field
  !! measurements by Donelan, Hamilton, and Hui (1985). Includes the high
  !! frequency form for beta_s found by Banner (1990). This function invokes the
  !!  DHH omnidirectional spectrum and the directional spreading functions to
  !! compute directional frequency spectrum:
  !!
  !! $$
  !!     F(f,\theta) = F'(f) * D(f,\theta)
  !! $$
  !!
  !! References:
  !!
  !! Donelan, M. A., J. Hamilton, and W. H. Hui, 1985. Directional
  !! spectra of wind-generated waves. Phil. Trans. Royal Soc. London A.,
  !! 315, 509-562.
  !!
  !! Banner, M. L., 1990. Equilibrium spectra of wind waves. J. Phys. Oceanogr.,
  !! 20, 966-984.

  real(kind=rk),dimension(:),intent(in) :: f !! Frequency [Hz]
  real(kind=rk),dimension(:),intent(in) :: theta !! Wave direction [rad]
  real(kind=rk),intent(in) :: wspd !! Wind speed at 10 m height [m/s]
  real(kind=rk),intent(in) :: fpeak !! Peak frequency [Hz]
  real(kind=rk),intent(in) :: theta_mean !! Mean wave direction [rad]
  real(kind=rk),intent(in) :: grav !! Gravitational acceleration [m/s^2]

  real(kind=rk),dimension(:,:),allocatable :: spec

  integer(kind=ik) :: ndir

  allocate(spec(size(f),size(theta)))

  do concurrent(ndir = 1:size(theta))
    spec(:,ndir) = donelanHamiltonHui(f,fpeak,wspd,grav)&
      * donelanHamiltonHuiDirectionalSpreading(f,wspd,fpeak,theta(ndir),&
                                               theta_mean)
  enddo

endfunction donelanHamiltonHuiDirectionalSpectrum
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function jonswap(f,wspd,fetch,grav) result(spec)

  !! Computes the JONSWAP equilibrium spectrum (Hasselmann et al. 1973) based on
  !!  input wind speed at the height of 10 m and fetch.
  !!
  !! References:
  !!
  !! Hasselmann, K. et al., 1973. Measurements of wind-wave growth and swell
  !! decay during the Joint North Sea Wave Project (JONSWAP). Dtsch. Hydrogh.
  !! Z., Suppl. A, 8, 12, 95pp.

  real(kind=rk),intent(in) :: f !! Frequency [Hz]
  real(kind=rk),intent(in) :: wspd !! Wind speed at 10 m height [m/s]
  real(kind=rk),intent(in) :: fetch !! Fetch [m]
  real(kind=rk),intent(in) :: grav !! Gravitational acceleration [m/s^2]

  real(kind=rk) :: alpha ! Phillips' constant
  real(kind=rk) :: sigma ! Spread of the spectrum peak
  real(kind=rk) :: r ! Exponent in peak enhancement factor

  real(kind=rk) :: fpeak ! Peak frequency [Hz]
  real(kind=rk) :: omega ! Radian frequency [rad/s]

  real(kind=rk),parameter :: beta = -1.25_rk
  real(kind=rk),parameter :: gamma = 3.3_rk ! peak enhancement

  omega = twopi*f
  alpha = 0.076_rk*nondimensionalFetch(wspd,fetch,grav)**(-0.22_rk)
  fpeak = jonswapPeakFrequency(wspd,fetch,grav)
  if(f > fpeak)then
    sigma = 0.09_rk
  else
    sigma = 0.07_rk
  endif
  r = exp(-0.5_rk*((f-fpeak)/(sigma*fpeak))**2)
  spec = twopi*alpha*grav**2/omega**5*exp(beta*(fpeak/f)**4)*gamma**r

endfunction jonswap
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function jonswapPeakFrequency(wspd,fetch,grav)&
  result(fpeak)

  !! Computes the JONSWAP equilibrium peak frequency [Hz] on the input
  !! based on the 10-m wind speed and fetch [km] (Hasselmann et al., 1973).
  !!
  !! References:
  !!
  !! Hasselmann, K. et al., 1973. Measurements of wind-wave growth and swell
  !! decay during the Joint North Sea Wave Project (JONSWAP). Dtsch. Hydrogh.
  !! Z., Suppl. A, 8, 12, 95pp.

  real(kind=rk),intent(in) :: wspd !! Wind speed at 10 m height [m/s]
  real(kind=rk),intent(in) :: fetch !! Fetch [m]
  real(kind=rk),intent(in) :: grav !! Gravitational acceleration [m/s^2]

  real(kind=rk),parameter :: a = 3.5_rk
  real(kind=rk),parameter :: b = 1._rk/3._rk

  fpeak = a*(grav**2/(wspd*fetch))**b

endfunction jonswapPeakFrequency
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function phillips(f,fpeak,grav) result(spec)

  !! Computes the Phillips (1958) equilibrium spectrum based on the input
  !! peak frequency [Hz].
  !!
  !! References:
  !!
  !! Phillips, O.M., 1958. The equilibrium range in the spectrum of
  !! wind-generated waves. J. Fluid Mech., 4, 426–434.
  !! doi:10.1017/S0022112058000550.

  real(kind=rk),intent(in) :: f !! Frequency [Hz]
  real(kind=rk),intent(in) :: fpeak !! Peak frequency [Hz]
  real(kind=rk),intent(in) :: grav !! Gravitational acceleration [m/s^2]

  real(kind=rk),parameter :: alpha = 8.13e-3 !! Phillips' parameter

  if(f < fpeak)then
    spec = 0
  else
    spec = alpha*grav**2*(twopi*f)**(-5)
  endif

endfunction phillips
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function piersonMoskowitz(f,wspd,grav) result(spec)

  !! Computes the Pierson-Moskowitz (1964) equilibrium spectrum based on input
  !! wind speed at the height of 10 m.
  !!
  !! References:
  !!
  !! Pierson Jr., W. J., and L. Moskowitz (1964), A proposed spectral form for
  !! fully developed wind seas based on the similarity theory of S. A.
  !! Kitaigorodskii, J. Geophys. Res., 69(24), 5181–5190,
  !! doi:10.1029/JZ069i024p05181.

  real(kind=rk),intent(in) :: f !! Frequency [Hz]
  real(kind=rk),intent(in) :: wspd !! Wind speed at 10 m height [m/s]
  real(kind=rk),intent(in) :: grav !! Gravitational acceleration [m/s^2]

  real(kind=rk) :: omega
  real(kind=rk) :: fpeak

  real(kind=rk),parameter :: alpha = 8.1e-3_rk ! Phillips' constant
  real(kind=rk),parameter :: beta = -1.25_rk

  omega = twopi*f
  fpeak = piersonMoskowitzPeakFrequency(wspd,grav)
  spec = twopi*alpha*grav**2/omega**5*exp(beta*(fpeak/f)**4)

endfunction piersonMoskowitz
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
pure elemental real(kind=rk) function piersonMoskowitzPeakFrequency(wspd,grav)&
  result(fpeak)

  !! Computes the Pierson-Moskowitz (1964) peak frequency based on input wind
  !! speed at the height of 10 m.
  !!
  !! References:
  !!
  !! Pierson Jr., W. J., and L. Moskowitz (1964), A proposed spectral form for
  !! fully developed wind seas based on the similarity theory of S. A.
  !! Kitaigorodskii, J. Geophys. Res., 69(24), 5181–5190,
  !! doi:10.1029/JZ069i024p05181.

  real(kind=rk),intent(in) :: wspd !! Wind speed at 10 m height [m/s]
  real(kind=rk),intent(in) :: grav !! Gravitational acceleration [m/s^2]
  real(kind=rk),parameter :: const = 0.1325_rk

  fpeak = const*grav/max(wspd,1e-2_rk)

endfunction piersonMoskowitzPeakFrequency
!-------------------------------------------------------------------------------
endmodule mod_spectral_shapes
