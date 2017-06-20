use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
type(spectrum_type) :: spec

! initialize a spectrum instance
spec = spectrum_type(fmin=0.04_rk,fmax=2._rk,df=1.1_rk,ndirs=1,depth=1000._rk)
end
