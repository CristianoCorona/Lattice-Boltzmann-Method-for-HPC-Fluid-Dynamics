cmake_minimum_required(VERSION 3.28)
project(LBM-01 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ----- Packages -----
find_package(OpenMP REQUIRED)

# --- Source files ----
set(SOURCE_FILES
    src/LBM/solver.cpp
    src/LBM/lattice.cpp
)

# ---- Executable ----
add_executable(${PROJECT_NAME}
    SOURCE_FILES
    src/main.cpp
)

# ---- Include dir ----
target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ----- Linking -----
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)

# ------ Tests -----
option(BUILD_TESTS "Build unit tests" OFF)

if(BUILD_TESTS)
    message(STATUS "Building unit tests with Catch2")
    enable_testing()
    add_subdirectory(test)
else()
    message(STATUS "Building withput unit tests, to enable them use -DBUILD_TESTS=ON")
endif()

