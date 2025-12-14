#include <iostream>
#include <array>
#include <LBM/LBM.hpp>

int main(int argc, char *argv[]) {
    std::array<int, D3Q19<float>::d> lattice_dim = {19, 19, 19};
    std::array<float, D3Q19<float>::d> lid_velocity = {1.0f, 0.0f, 0.0f};
    float nu = 1.0f;

    Lattice<D3Q19<float>, float> lattice(
            lattice_dim,
            lid_velocity,
            nu,
            Direction<D3Q19<float>::d>::BOTTOM);
    
    Solver<D3Q19<float>, float> solver(lattice);

    unsigned long n_iter = 10000;
    float delta_t = 1.0f;
    solver.solve(n_iter, delta_t);

    return 0;
}
