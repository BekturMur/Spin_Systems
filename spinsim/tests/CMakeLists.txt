include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.7.1
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(Catch2)

add_executable(spinsim_tests
  DenseReference.cpp
  test_SpinOperators.cpp
  test_Hamiltonian.cpp
  test_StateVector.cpp
  test_Bessel.cpp
  test_Propagator.cpp
  test_Observable.cpp
  test_Gate.cpp
  test_Ensemble.cpp
  test_Config.cpp)

target_include_directories(spinsim_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(spinsim_tests PRIVATE
  SPINSIM_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_link_libraries(spinsim_tests PRIVATE
  spinsim_core spinsim_lapack spinsim_options Catch2::Catch2WithMain)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(spinsim_tests)
