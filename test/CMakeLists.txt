include(FetchContent)

# ---- Get Catch2 directly from GitHub ----
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.11.0
)
FetchContent_MakeAvailable(Catch2)

# ----- Test source files ------
set(SOURCE_FILES
    lattice_test_D2Q9.cpp 
)

# ---- Create test executable ----
add_executable(run_tests
    SOURCE_FILES
)

# ---- Link SFLA library and Catch2 ----
target_link_libraries(run_tests PRIVATE
    ${PROJECT_NAME}
    Catch2::Catch2WithMain
)

# ---- Automatically find tests for Catch2 ----
include(Catch)
catch_discover_tests(run_tests)
