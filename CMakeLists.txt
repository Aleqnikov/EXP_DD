cmake_minimum_required(VERSION 3.10)
project(DD_Math)

set(CMAKE_CXX_STANDARD 17)

add_executable(dd_calc
        main.cpp
        DD.cpp
        low_level_math.cpp
)


if(NOT MSVC)
    target_compile_options(dd_calc PRIVATE -O3 -fno-fast-math)
else()
    target_compile_options(dd_calc PRIVATE /fp:strict)
endif()