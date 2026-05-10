cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0148 NEW)

project(eclipse LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)

# library
add_library(eclipse INTERFACE)
target_include_directories(eclipse INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(eclipse PROPERTIES POSITION_INDEPENDENT_CODE ON)


# compile main executable
add_executable(eclipse_main
    main.cxx
)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

target_include_directories(eclipse_main PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(eclipse_main PRIVATE ${Python3_LIBRARIES})

find_package(PkgConfig REQUIRED)
pkg_check_modules(ERFA REQUIRED erfa)

target_include_directories(eclipse_main PRIVATE ${ERFA_INCLUDE_DIRS})
target_link_libraries(eclipse_main PRIVATE ${ERFA_LIBRARIES})


