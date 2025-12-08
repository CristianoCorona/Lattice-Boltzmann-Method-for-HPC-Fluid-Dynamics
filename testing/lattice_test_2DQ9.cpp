#define CATCH_CONFIG_MAIN

#include "catch2.hpp"
#include "lattice.hpp"

using float_type = double;
using Descriptor = D2Q9<float_type>;


TEST_CASE("Lattice initialization and indexing", "[lattice]") {

    // test parameters
    std::array<int, Descriptor::d> dimensions = {10, 5}; 
    float_type lid_velocity = 0.1;
    float_type nu = 0.01;

    // Create lattice
    Lattice<Descriptor, float_type> lattice(
        dimensions, 
        lid_velocity, 
        nu,
        "test_output.vtk");

    SECTION("Dimension correctness") {
        REQUIRE(lattice.sizes[0] == 10);
        REQUIRE(lattice.sizes[1] == 5);
        REQUIRE(lattice.u_lid == lid_velocity);
        REQUIRE(lattice.nu == nu);
        REQUIRE(lattice.total_cells == 50);
    }

    SECTION("Buffer sizes") {
        for(int i = 0; i < Descriptor::q; ++i) {
            REQUIRE(lattice.f[i].size() == lattice.total_cells);
        }
        for(int k = 0; k < Descriptor::d; ++k) {
            REQUIRE(lattice.u[k].size() == lattice.total_cells);
        }
        REQUIRE(lattice.rho.size() == lattice.total_cells);
    }

    SECTION("Strides correctness") {
        // x stride = 1, y stride = nx
        REQUIRE(lattice.strides[0] == 1);       // x stride
        REQUIRE(lattice.strides[1] == 10);      // y stride
    }
}



