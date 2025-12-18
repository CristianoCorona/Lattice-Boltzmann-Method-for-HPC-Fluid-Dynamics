#include <iostream>
#include <array>
#include <LBM/LBM.hpp>

using namespace std;

int main(int argc, char *argv[]) {


    array<int, D3Q19<float>::d> lattice_dim = {81, 81, 81};
    array<float, D3Q19<float>::d> lid_velocity = {0.1f, 0.0f, 0.0f};
    float nu = 1.0f;

    Lattice<D3Q19<float>, float> lattice(
            lattice_dim,
            lid_velocity,
            nu,
            Direction<D3Q19<float>::d>::TOP);
    
    Solver<D3Q19<float>, float> solver(lattice);

    unsigned long n_iter = 1000;
    
    solver.solve(n_iter, 10);

    return 0;
}
