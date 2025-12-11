cmake_minimum_required(VERSION 3.28)
project(LBM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ----- Packages -----
find_package(OpenMP REQUIRED)

# ---- Source subdirectory ----
add_subdirectory(src/LBM)

# ---- Executable ----
add_executable(${PROJECT_NAME}
    src/main.cpp
)

# ---- Link LBM_lib ----
target_link_libraries(${PROJECT_NAME} PRIVATE LBM_lib)

# ---- Include dir ----
target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ------ Tests -----
option(BUILD_TESTS "Build unit tests" OFF)

if(BUILD_TESTS)
    message(STATUS "Building unit tests with Catch2")
    enable_testing()
    add_subdirectory(test)
else()
    message(STATUS "Building withput unit tests, to enable them use -DBUILD_TESTS=ON")
endif()

