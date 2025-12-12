cmake_minimum_required(VERSION 3.31)
project(matrix_c LANGUAGES C)
project(matrix_c VERSION 0.0.0)

set(CMAKE_C_STANDARD 23)

include_directories(/usr/local/include inc)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(PROJECT_SOURCE_DIR "src")
set(PROJECT_INCLUDE_DIR "inc")
set(PROJECT_TEST_DIR "test")


add_compile_options(
        # this flag adds pre processing to format strings in the case of NEGATIVE color effect (adds 1 whole copy)
        -DNEGATIVE_NEWLINE_CORRECTION
        -flto
)

set(SOURCE_FILES
    ${PROJECT_SOURCE_DIR}/matrix.c
)
add_library(matrix_static STATIC ${SOURCE_FILES})
add_library(matrix_shared SHARED ${SOURCE_FILES})

add_executable(matrix_c_static_test ${PROJECT_TEST_DIR}/main.c)
add_executable(matrix_c_shared_test ${PROJECT_TEST_DIR}/main.c)

target_link_libraries(matrix_c_static_test PRIVATE cprintf_static)
target_link_libraries(matrix_c_shared_test PRIVATE cprintf_shared)
