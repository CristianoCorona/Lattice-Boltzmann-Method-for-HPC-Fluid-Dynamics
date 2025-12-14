#include <iostream>
#include <array>
#include <LBM/LBM.hpp>

int main(int argc, char *argv[]) {
    std::array<int, D2Q9<float>::d> lattice_dim = {81, 81};
    std::array<float, D2Q9<float>::d> lid_velocity = {0.1f, 0.0f};
    float nu = 0.1f*81.0f/100.0f; // from Re=100

    Lattice<D2Q9<float>, float> lattice(
            lattice_dim,
            lid_velocity,
            nu,
            Direction<D2Q9<float>::d>::TOP);
    
    Solver<D2Q9<float>, float> solver(lattice);

    unsigned long n_iter = 10000;
    float delta_t = 1.0f;
    solver.solve(n_iter, delta_t);

    return 0;
}
