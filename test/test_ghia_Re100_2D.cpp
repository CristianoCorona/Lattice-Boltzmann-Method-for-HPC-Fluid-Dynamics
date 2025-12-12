#define CATCH_CONFIG_MAIN

#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "LBM/lattice.hpp"
#include "LBM/solver.hpp"

using float_type = double;
using Descriptor = D2Q9<float_type>;

TEST_CASE("Ghia 1982 Benchmark (Re=100, 2D lid-driven cavity flow)", "[ghia_re100_2d]") {

    // Parameters setup
    const int N = 129; //Odd grid size to have exact center -> 129 in the paper (81 might work and be faster)
    std::array<int, Descriptor:: d> dimensions = {N, N};

    std::array<float_type, Descriptor::d> lid_velocity = {1.0, 0.0}; // Lid velocity = 1.0 for normalized units
    float_type Re = 100.0;
    float_type L = static_cast<float_type>(N - 1); // Domain length in lattice units
    float_type nu = lid_velocity[0] * L / Re; // Kinematic viscosity (from Re definition) TODO -> norm of lid_velocity

    //float_type tau = 3.0 * nu + 0.5; // Relaxation time

    // Lattice creation
    Lattice<Descriptor, float_type> lattice(dimensions, 
                                            lid_velocity,
                                            nu,
                                            Direction<Descriptor::d>::TOP,
                                            "ghia_re100_2d.vtk");
    // Solver creation
    Solver<Descriptor, float_type> solver(lattice);
    // Initialize equilibrium
    lattice.initialize_equilibrium();

    // Simulation
    const unsigned long max_steps = 5000;
    solver.solve(max_steps, 10.0);

    // Validation at centerlines
    const int center_x = N / 2;
    const int center_y = N / 2;
    const int center_idx = lattice.idx(center_x, center_y); 

    // Simulation results at center
    float_type u_center = lattice.u[0][center_idx]/lid_velocity[0]; // u at (N/2, N/2), normalized by lid velocity 
    float_type v_center = lattice.u[1][center_idx]/lid_velocity[0]; // v at (N/2, N/2), normalized by lid velocity

    // Expected values from Ghia 1982 paper at Re=100
    const float_type u_ref = -0.20581;
    const float_type v_ref = 0.05454; 
    
    std::cout << "Re=100 Validation Results at Center (" << center_x << "," << center_y << "):\n";
    std::cout << "U_sim (norm): " << u_center << " | Ghia: " << u_ref << "\n";
    std::cout << "V_sim (norm): " << v_center << " | Ghia: " << v_ref << "\n";

    // Assertions with a tolerance
    REQUIRE_THAT(u_center, Catch::Matchers::WithinAbs(u_ref, 0.05));
    REQUIRE_THAT(v_center, Catch::Matchers::WithinAbs(v_ref, 0.05));
    }