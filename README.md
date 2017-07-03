## wavy

A spectral ocean wave modeling framework.

* [Getting started](#getting-started)
  - [Building wavy](#building-wavy)
  - [Building a simple program with wavy](#building-a-simple-program-with-wavy)
  - [Examples](#examples)
* [Design principles](#design-principles)
* [Features](#features)
* [Development status](#development-status)

### Getting started

#### Building wavy

```
git clone --recursive https://github.com/wavebitscientific/wavy.git
FC=gfortran make
mkdir build
cd build
FC=gfortran cmake ..
make
ctest
```
Change the `FC` value if building with a different Fortran compiler.

By default, wavy will be built in single precision (32-bit reals).
To build it in double (64-bit reals) or quadruple precision (128-bit reals),
use the `-DREAL` argument when invoking `cmake`:
```
FC=gfortran cmake .. -DREAL=64 # for double precision
```
or:
```
FC=gfortran cmake .. -DREAL=128 # for quad precision
```

wavy needs gcc-6.3.0 or later to succesfully build and pass all tests.

#### Building a simple program with wavy

If you compiled wavy in `wavy/build`, then the module files and the library
are located in `wavy/build/include` and `wavy/build/lib`, respectively. 
For example, if we want to build a simple wavy hello world program from
the base directory, we could do it like this: 

```
gfortran hello.f90 -o hello -Ibuild/include -Lbuild/lib -lwavy
```

#### Examples

Initialize a omnidirectional spectrum instance in the frequency range from 
0.04 to 2 Hz with logarithmic increment of 1.1, in mean water depth of 1000 m:

```fortran
use mod_spectrum,only:spectrum_type
type(spectrum_type) :: spec

! initialize a spectrum instance
spec = spectrum_type(fmin=0.04,fmax=2.,df=1.1,ndirs=1,depth=1000.)
```

Same as above, but directional spectrum with 36 directional bins:

```fortran
spec = spectrum_type(fmin=0.04,fmax=2.,df=1.1,ndirs=36,depth=1000.)
```

Initialize omnidirectional spectrum with JONSWAP shape at wind speed of 10 m/s,
and fetch of 100 km:

```fortran
use mod_spectrum,only:spectrum_type
use mod_spectral_shapes,only:jonswap

type(spectrum_type) :: spec

! initialize a spectrum instance
spec = spectrum_type(fmin=0.04,fmax=2.,df=1.1,ndirs=1,depth=1000.)

! assign a JONSWAP-shape spectrum to the instance
spec = jonswap(spec % getFrequency(),wspd=10.,fetch=1e5,grav=9.8)
```

Above examples will work with default precision (`REAL32`). 
To write code that is always compatible with precision specified at 
build time, use `mod_precision` module:

```fortran
use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_spectral_shapes,only:jonswap

type(spectrum_type) :: spec

! initialize a spectrum instance
spec = spectrum_type(fmin=0.04_rk,fmax=2._rk,df=1.1_rk,ndirs=1,depth=1000._rk)

! assign a JONSWAP-shape spectrum to the instance
spec = jonswap(spec % getFrequency(),wspd=10._rk,fetch=1e5_rk,grav=9.8_rk)
``` 
There are many pre-built diagnostics that can be output from a `spectrum`
instance, here is a taste of a few:

```fortran
write(*,*)'Significant wave height [m]: ',spec % significantWaveHeight()
write(*,*)'       Mean wave period [s]: ',spec % meanPeriod()
write(*,*)'      Mean square slope [-]: ',spec % meanSquareSlope()
write(*,*)'         Stokes drift [m/s]: ',spec % stokesDrift([0._rk])
```
outputs:
```
 Significant wave height [m]:    2.13949418    
        Mean wave period [s]:    5.20506239    
       Mean square slope [-]:    2.39831898E-02
          Stokes drift [m/s]:    4.87080999E-02
```

### Design principles

  * Pure Fortran goodness
  * Object-oriented for high-level abstractions, functional for the computational kernels
  * Framework to construct arbitrary wave models
  * Provides a rich set of source functions, parametric shapes, numerical schemes, and wave diagnostics
  * Easy to use by existing atmosphere and ocean models
  * Supports single (32-bit), double (64-bit), and quadruple (128-bit)  precision floating point arithmetic
  * Self-documented
  * Fully unit-tested

### Features

* Classes:
    - [x] Spectrum class: supports omnidirectional and directional spectra
    - [x] Grid class: supports regular 1-d and 2-d grids
    - [x] Domain class: built from Grid and Spectrum instances.

* Source functions:
    - [x] Sin, Sds, Sdt, Snl, Sbf (Donelan et al., 2012)
    - [ ] Sin, Sds, WAM cycle 3 (Komen et al., 1984)
    - [ ] Sin, Sds (Tolman and Chalikov, 1996)
    - [ ] Sin, Sds, WAM cycle 4 (Janssen, 2004)
    - [ ] Sin, Sds, (Ardhuin et al., 2010) 
    - [ ] Snl, DIA (Hasselmann et al., 1985)

* Parametric spectra:
    - [x] Donelan, Hamilton & Hui (1985)
    - [x] JONSWAP (1973)
    - [x] Phillips (1956)
    - [x] Pierson-Moskowitz (1964)

* Grid projections:
    - [x] Cartesian
    - [ ] Spherical

* Time integration scheme:
    - [x] First order, forward Euler (explicit)
    - [x] First order, backward Euler (implicit)
    - [x] Exact exponential

* Advection schemes:
    - [x] Upwind differences, 1st order, 1-d
    - [x] Upwind differences, 1st order, 2-d
    - [x] Upwind differences, 2nd order, 1-d
    - [x] Upwind differences, 2nd order, 2-d
    - [x] Centered differences, 2nd order, 1-d
    - [x] Centered differences, 2nd order, 2-d
    - [ ] Flux corrected transport (Zalesak, 1979)
    - [ ] Smolarkiewitz anti-diffusive velocity scheme (Smolarkiewitz, 1981)
    - [ ] MPDATA (Smolarkiewitz, 1983)
    - [ ] ULTIMATE QUICKEST (Wavewatch III)

* Stokes drift:
    - [ ] Monochromatic
    - [x] Exact
    - [ ] Webb and Fox-Kemper (2011)
    - [ ] Breivik et al. (2015)

* Parallel execution:
    - [ ] Fortran Coarrays

* Input/Output:
    - [ ] NetCDF
    - [x] JSON

### Development status

wavy is in early development. Contributors are highly needed,
especially for implementing source functions. You can start 
contributing by opening an [issue](https://github.com/wavebitscientific/wavy/issues/new).
