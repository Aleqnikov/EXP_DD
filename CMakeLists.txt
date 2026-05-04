cmake_minimum_required(VERSION 3.14)
project(DD_Math)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

FetchContent_Declare(
        CLI11
        URL https://github.com/CLIUtils/CLI11/releases/download/v2.4.1/CLI11.hpp
        DOWNLOAD_NO_EXTRACT TRUE
)
FetchContent_MakeAvailable(CLI11)

add_executable(dd_calc
        main.cpp
        DD.cpp
        low_level_math.cpp
)


target_include_directories(dd_calc PRIVATE ${cli11_SOURCE_DIR})


if(NOT MSVC)
    target_compile_options(dd_calc PRIVATE -O3 -fno-fast-math)
else()
    target_compile_options(dd_calc PRIVATE /fp:strict)
endif()