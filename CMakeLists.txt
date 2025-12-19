cmake_minimum_required(VERSION 3.28)
project(LBM LANGUAGES CXX)

# ---------------- Build type ----------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ---------------- C++ standard ----------------
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------- Optimization (Release only) ----------------
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Configuring Release optimizations")

    # Explicit optimization flags
    add_compile_options(
        -O3
        -march=native
        -ffast-math
    )

    # Enable LTO (compile + link)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
endif()

# ---------------- OpenMP ----------------
find_package(OpenMP REQUIRED)

# ---------------- Source subdirectory ----------------
add_subdirectory(src/LBM)

# ---------------- Executables ----------------
add_executable(lbm3D src/main3D.cpp)
add_executable(lbm2D src/main2D.cpp)
add_executable(lbm3D_cluster src/cluster_main.cpp)

# ---------------- Link libraries ----------------
target_link_libraries(lbm3D PRIVATE LBM_lib OpenMP::OpenMP_CXX)
target_link_libraries(lbm2D PRIVATE LBM_lib OpenMP::OpenMP_CXX)
target_link_libraries(lbm3D_cluster PRIVATE LBM_lib OpenMP::OpenMP_CXX)

# ---------------- Tests ----------------
option(BUILD_TESTS "Build unit tests" OFF)

if(BUILD_TESTS)
    message(STATUS "Building unit tests with Catch2")
    enable_testing()
    add_subdirectory(test)
else()
    message(STATUS "Building without unit tests (enable with -DBUILD_TESTS=ON)")
endif()
