# ----- Library sources ----
set(SOURCE_FILES
    solver.cpp
    lattice.cpp
    boundary.cpp
)

# ----- Create LBM_lib ----
add_library(LBM_lib ${SOURCE_FILES})

# ---- Include directory ----
target_include_directories(LBM_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

# ----- Link OpenMP -----
target_link_libraries(LBM_lib PUBLIC OpenMP::OpenMP_CXX)
