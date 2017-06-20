use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_spectral_shapes,only:jonswap

type(spectrum_type) :: spec

! initialize a spectrum instance
spec = spectrum_type(fmin=0.04_rk,fmax=2._rk,df=1.1_rk,ndirs=1,depth=1000._rk)

! assign a JONSWAP-shape spectrum to the instance
spec = jonswap(spec % getFrequency(),wspd=10._rk,fetch=1e5_rk,grav=9.8_rk)

write(*,*)'Significant wave height [m]: ',spec % significantWaveHeight()
write(*,*)'       Mean wave period [s]: ',spec % meanPeriod()
write(*,*)'      Mean square slope [-]: ',spec % meanSquareSlope()
write(*,*)'         Stokes drift [m/s]: ',spec % stokesDrift([0._rk])

end
