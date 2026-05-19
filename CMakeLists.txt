cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0148 NEW)

project(eclipse LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -DNDEBUG")

# include dirs
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)

# find libraries
find_package(PkgConfig REQUIRED)
pkg_check_modules(ERFA REQUIRED erfa)
find_package(Eigen3 REQUIRED)


# compile main executable
add_executable(eclipse_main main.cxx)

# Link ERFA and Eigen to the executable
target_include_directories(eclipse_main PRIVATE ${ERFA_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(eclipse_main PRIVATE ${ERFA_LIBRARIES} Eigen3::Eigen)


