#include <iostream>
#include <array>
#include <LBM/LBM.hpp>

int main(int argc, char *argv[]) {
    std::array<int, D2Q9<float>::d> lattice_dim = {100, 100};
    std::array<float, D2Q9<float>::d> lid_velocity = {1.0f, 0.0f};
    float nu = 1.0f;

    Lattice<D2Q9<float>, float> lattice(
            lattice_dim,
            lid_velocity,
            nu,
            Direction<D2Q9<float>::d>::TOP);

    Solver<D2Q9<float>, float> solver(lattice);

    unsigned long n_iter = 1000;
    float delta_t = 0.1f;
    solver.solve(n_iter, delta_t);

    return 0;
}
