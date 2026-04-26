cmake_minimum_required (VERSION 3.22)

project(cfextract
  LANGUAGES CXX
  VERSION 0.0.0
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(MAKE_PY "Compile python bindings" OFF)
option(MAKE_TEST "Compile and run tests binary" OFF)

option(ENABLE_DEBUG_FLAGS "Enable extra debug compile options" OFF)

# user paths to htslib if needed
set(HTSLIB_INCLUDE_DIR "" CACHE PATH "Path to directory that contains htslib/hts.h")
set(HTSLIB_LIBRARY "" CACHE FILEPATH "Absolute path to libhts.{so,dylib,a}")

# //--- shared deps ---//
find_package(PkgConfig QUIET)

set(HTSLIB_TARGET "")
set(_htslib_inc "")
set(_htslib_lib "")

if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
  message(STATUS "Using user-specified htslib install. Must be >=1.14, but this is not enforced.")
  set(_htslib_inc "${HTSLIB_INCLUDE_DIR}")
  set(_htslib_lib "${HTSLIB_LIBRARY}")
endif()

if(_htslib_inc AND _htslib_lib)
  add_library(htslib UNKNOWN IMPORTED GLOBAL)
  add_library(htslib::htslib ALIAS htslib)
  set_target_properties(htslib PROPERTIES
    IMPORTED_LOCATION "${_htslib_lib}"
    INTERFACE_INCLUDE_DIRECTORIES "${_htslib_inc}"
  )
  set(HTSLIB_TARGET htslib::htslib)
endif()

if(NOT HTSLIB_TARGET AND PkgConfig_FOUND)
  # assuming the pkg-config module is called "htslib"
  pkg_check_modules(HTSLIB QUIET IMPORTED_TARGET htslib>=1.14)
  if(HTSLIB_FOUND)
    message(STATUS "Located htslib>=1.14 install via pkg-config.")
    set(HTSLIB_TARGET PkgConfig::HTSLIB)
  endif()
endif()

if(NOT HTSLIB_TARGET)
  message(FATAL_ERROR
    "Could not find htslib directly or via pkg-config. "
    "Set both HTSLIB_INCLUDE_DIR and HTSLIB_LIBRARY to use user install."
  )
endif()
# //--- end deps ---//


# shared core lib to reduce problem surface
add_library(common_compile_opts INTERFACE)
target_compile_options(common_compile_opts INTERFACE
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
    -Wsign-conversion
    -Wfloat-equal
    -Wnon-virtual-dtor
    -Wold-style-cast
    -Woverloaded-virtual
    -Wnull-dereference
    -Wdouble-promotion
    -Wformat=2
)

if(ENABLE_DEBUG_FLAGS)
  target_compile_options(common_compile_opts INTERFACE
    -fsanitize=address,undefined
    -fno-omit-frame-pointer
    -O0 -g
  )
  target_link_options(common_compile_opts INTERFACE
    -fsanitize=address,undefined
  )
endif()


add_library(cfex-core-shared OBJECT
  src/core/access.cpp
  src/core/extract.cpp
  src/core/features.cpp
  src/core/stats.cpp
)
set_target_properties(cfex-core-shared PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(cfex-core-shared PUBLIC
  src
)
target_link_libraries(cfex-core-shared PUBLIC
  ${HTSLIB_TARGET}
)

# for target specifc deps
Include(FetchContent)

if(MAKE_TEST)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.8.1
  )

  FetchContent_MakeAvailable(Catch2)

  add_executable(test-cfextract
    tests/cpp/test_extract.cpp
  )
  target_include_directories(test-cfextract PRIVATE src)
  target_link_libraries(test-cfextract PRIVATE
    Catch2::Catch2WithMain
    cfex-core-shared
  )
endif(MAKE_TEST)


# TODO - possibly remove stub generation
# from install process. It is quite finicky
# for an interview project and would be a dreadful
# thing to fail on.
if(MAKE_PY)
  set(PYBIND11_FINDPYTHON NEW)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v3.0.3
  )
  FetchContent_MakeAvailable(pybind11)
  pybind11_add_module(cfextract
    src/bindings/python/bindings.cpp
  )
  target_link_libraries(cfextract PRIVATE
    common_compile_opts
    cfex-core-shared
  )
  set_target_properties(cfextract PROPERTIES
    OUTPUT_NAME cfextract
  )
  # janky, to revist
  find_program(PYBIND11_STUBGEN_EXECUTABLE pybind11-stubgen)
  if(PYBIND11_STUBGEN_EXECUTABLE)
    set(_stubgen_script "${CMAKE_CURRENT_BINARY_DIR}/run_stubgen.cmake")
    file(GENERATE
      OUTPUT "${_stubgen_script}"
      CONTENT
"set(ENV{PYTHONPATH} \"$<TARGET_FILE_DIR:cfextract>:$ENV{PYTHONPATH}\")
execute_process(
  COMMAND \"${PYBIND11_STUBGEN_EXECUTABLE}\" cfextract
    --output-dir \"${CMAKE_CURRENT_BINARY_DIR}/stubs\"
  RESULT_VARIABLE result
)
if(NOT result EQUAL 0)
  message(STATUS \"pybind11-stubgen failed (exit \${result}) -- cfextract.pyi will not be installed\")
endif()
"
    )
    add_custom_command(TARGET cfextract POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/stubs"
      COMMAND ${CMAKE_COMMAND} -P "${_stubgen_script}"
      COMMENT "Generating Python stubs for cfextract"
    )
  else()
    message(STATUS "pybind11-stubgen not found, skipping stub generation")
  endif()
  # for pip install
  install(TARGETS cfextract
    LIBRARY DESTINATION .
  )
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/stubs/cfextract.pyi
    DESTINATION .
    OPTIONAL
  )
endif()


