# cmake version, project name, language
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(wavy)
enable_language(Fortran)

# set output paths for modules, archives, and executables
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# if build type not specified, default to release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "release")
endif()

# hande precision
if(REAL)
  message(STATUS "Configuring build for ${REAL}-bit reals")
  add_definitions(-DREAL${REAL})
  include_directories(src/external/json-fortran/build_real${REAL}/include)
  link_directories(src/external/json-fortran/build_real${REAL}/lib)
else()
  message(STATUS "Configuring build for 32-bit reals")
  add_definitions(-DREAL32)
  include_directories(src/external/json-fortran/build_real32/include)
  link_directories(src/external/json-fortran/build_real32/lib)
endif()

# compiler flags for gfortran
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -C -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

# compiler flags for ifort
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -assume realloc_lhs -heap-arrays")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -C -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

# compiler flags for Cray ftn
if(CMAKE_Fortran_COMPILER_ID MATCHES Cray)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -h noomp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

foreach(extlib datetime functional)
  include_directories(src/external/${extlib}-fortran/build/include)
  link_directories(src/external/${extlib}-fortran/build/lib)
endforeach()

file(GLOB wavy_lib "src/lib/mod_*.f90")

# library to archive (libwavy.a)
add_library(wavy ${wavy_lib})
add_library(testing src/tests/unit/mod_testing.f90)

# tests
enable_testing()

# unit tests
foreach(exec json utility grid spectrum domain)
  add_executable(test_${exec} src/tests/unit/test_${exec}.f90)
  target_link_libraries(test_${exec} wavy testing datetime functional jsonfortran)
  add_test(test_${exec} bin/test_${exec})
endforeach()

# integration tests
foreach(exec jonswap phillips pierson_moskowitz dhh dhh_directional advection advection_2d time_integration_spectrum time_integration_domain) 
  add_executable(test_${exec} src/tests/integration/test_${exec}.f90)
  target_link_libraries(test_${exec} wavy testing datetime functional jsonfortran)
  add_test(test_${exec} bin/test_${exec})
endforeach()

foreach(exec umwm_1d spectrum_omnidirectional spectrum_directional spectrum_jonswap)  
  add_executable(${exec} src/examples/${exec}.f90)
  target_link_libraries(${exec} wavy testing datetime functional jsonfortran)
endforeach()

install(DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}" DESTINATION "${CMAKE_INSTALL_PREFIX}")
install(DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" DESTINATION "${CMAKE_INSTALL_PREFIX}")
