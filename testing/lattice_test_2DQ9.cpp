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

    SECTION("Indexing function") {
        // Testing some known indices
        // (0,0) -> 0 (first cell)
        REQUIRE(lattice.idx(0, 0) == 0);
        // (1,0) -> 1
        REQUIRE(lattice.idx(1, 0) == 1);
        // (0,1) -> 10
        REQUIRE(lattice.idx(0, 1) == 10);
        // (5,3) -> 35
        REQUIRE(lattice.idx(5, 3) == 35);
        // (9,4) -> 49 (last cell)
        REQUIRE(lattice.idx(9, 4) == 49);
    }

    SECTION("Precomputed neighbor offsets correctness") {
        // D2Q9 expected offsets:
        // c_0 -> 0
        REQUIRE(lattice.neighbor_offsets[0] == 0);    // c_0
        // c_1 -> 1
        REQUIRE(lattice.neighbor_offsets[1] == 1);    // c_1
        // c_2 -> nx = 10
        REQUIRE(lattice.neighbor_offsets[2] == dimensions[0]);   // c_2
        // c_3 -> -1
        REQUIRE(lattice.neighbor_offsets[3] == -1);   // c_3
        // c_4 -> -nx = -10 
        REQUIRE(lattice.neighbor_offsets[4] == -dimensions[0]);  // c_4
        // c_5 -> nx + 1 = 11
        REQUIRE(lattice.neighbor_offsets[5] == dimensions[0] + 1);   // c_5
        // c_6 -> nx - 1 = 9
        REQUIRE(lattice.neighbor_offsets[6] == dimensions[0] - 1);    // c_6
        // c_7 -> -nx - 1 = -11
        REQUIRE(lattice.neighbor_offsets[7] == -dimensions[0] - 1);  // c_7
        // c_8 -> -nx + 1 = -9
        REQUIRE(lattice.neighbor_offsets[8] == -dimensions[0] + 1);   // c_8
        
    
    }
}



