cmake_minimum_required(VERSION 3.28)
project(spinsim VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(SPINSIM_BUILD_TESTS "Build the test suite" ON)
option(SPINSIM_WERROR "Treat warnings as errors" OFF)

add_library(spinsim_options INTERFACE)
target_compile_options(spinsim_options INTERFACE
  -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
  $<$<BOOL:${SPINSIM_WERROR}>:-Werror>
  # -ffp-contract=off keeps fused multiply-add out of the kernels, so a run is
  # reproducible across optimisation levels and comparable against the exact
  # matrix exponential at the tolerances the tests assert.
  -ffp-contract=off)

# LAPACK/BLAS: Accelerate on macOS, whatever CMake finds elsewhere. Used for
# reference eigendecompositions in tests and for the small-subsystem ground
# state (zheev in chinitPDDGnmr.f).
if(APPLE)
  find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
  add_library(spinsim_lapack INTERFACE)
  target_link_libraries(spinsim_lapack INTERFACE ${ACCELERATE_FRAMEWORK})
  target_compile_definitions(spinsim_lapack INTERFACE SPINSIM_ACCELERATE_LAPACK)
else()
  find_package(LAPACK REQUIRED)
  add_library(spinsim_lapack INTERFACE)
  target_link_libraries(spinsim_lapack INTERFACE LAPACK::LAPACK)
endif()

include(FetchContent)
FetchContent_Declare(
  tomlplusplus
  GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
  GIT_TAG v3.4.0
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(tomlplusplus)

add_library(spinsim_core
  src/core/StateVector.cpp
  src/core/SpinOperators.cpp
  src/core/Hamiltonian.cpp
  src/core/Bessel.cpp
  src/core/Chebyshev.cpp
  src/core/Propagator.cpp
  src/core/Observable.cpp
  src/core/Gate.cpp
  src/sequence/Sequence.cpp
  src/ensemble/Rng.cpp
  src/ensemble/Ensemble.cpp
  src/ensemble/EnsembleRunner.cpp
  src/io/Config.cpp
  src/io/Experiment.cpp)
target_include_directories(spinsim_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(spinsim_core PUBLIC tomlplusplus::tomlplusplus)
target_link_libraries(spinsim_core PRIVATE spinsim_options)

add_executable(bench_propagator tests/bench_propagator.cpp)
target_link_libraries(bench_propagator PRIVATE spinsim_core spinsim_options)

add_executable(spinsim apps/spinsim/main.cpp)
target_link_libraries(spinsim PRIVATE spinsim_core spinsim_options)

if(SPINSIM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
