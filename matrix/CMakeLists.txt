cmake_minimum_required(VERSION 3.31)
project(matrix LANGUAGES C CXX ASM_NASM)
project(matrix VERSION 0.0.0)
enable_language(C CXX ASM_NASM)

set(CMAKE_C_STANDARD 23)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_LIST_DIR}/cmake/nasm-toolchain.cmake" CACHE STRING "")

include_directories(/usr/local/include inc)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(PROJECT_SOURCE_DIR "src")
set(PROJECT_INCLUDE_DIR "inc")
set(PROJECT_TEST_DIR "test")

set(CMAKE_CXX_FLAGS "-flto -funroll-loops -ffast-math -mavx -mavx2 -mfma -O3")
add_compile_options(
    -DAVX256
    -g
)

set(SOURCE_FILES
)

#add_library(matrix_static STATIC ${SOURCE_FILES})
#add_library(matrix_shared SHARED ${SOURCE_FILES})

add_executable(matrix_cpp_static_test ${PROJECT_TEST_DIR}/main.cpp)
add_executable(matrix_cpp_shared_test ${PROJECT_TEST_DIR}/main.cpp)

#target_link_libraries(matrix_cpp_static_test PRIVATE matrix_static)
#target_link_libraries(matrix_cpp_shared_test PRIVATE matrix_shared)
