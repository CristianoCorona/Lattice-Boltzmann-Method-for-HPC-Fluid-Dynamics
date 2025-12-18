#include <iostream>
#include <array>
#include <LBM/LBM.hpp>

using float_type = double;
int main(int argc, char *argv[]) {
    std::array<int, D2Q9<float_type>::d> lattice_dim = {257, 257};
    std::array<float_type, D2Q9<float_type>::d> lid_velocity = {0.1, 0.0};
    float_type Re = 1000.0;
    float_type L = static_cast<float_type>(lattice_dim[0] - 1);
    float_type nu = lid_velocity[0] * L / Re;

    Lattice<D2Q9<float_type>, float_type> lattice(
            lattice_dim,
            lid_velocity,
            nu,
            Direction<D2Q9<float_type>::d>::TOP,
            "lbm_cluster_output.vtk");
    
    Solver<D2Q9<float_type>, float_type> solver(lattice);

    unsigned long n_iter = 200000;
    int output_interval = 5000;
    solver.solve(n_iter, output_interval);

    return 0;
}