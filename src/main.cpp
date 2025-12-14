#include <iostream>
#include <array>
#include <LBM/LBM.hpp>

int main(int argc, char *argv[]) {
<<<<<<< Updated upstream
    std::array<int, D2Q9<float>::d> lattice_dim = {129, 129};
    std::array<float, D2Q9<float>::d> lid_velocity = {0.1f, 0.0f};
    float Re = 100.0f;
    float L = static_cast<float>(lattice_dim[0] - 1);
    float nu = lid_velocity[0] * L / Re;
=======
    std::array<int, D2Q9<float>::d> lattice_dim = {1000, 1000};
    std::array<float, D2Q9<float>::d> lid_velocity = {0.1f, 0.0f};
    float nu = 1.0f;
>>>>>>> Stashed changes

    Lattice<D2Q9<float>, float> lattice(
            lattice_dim,
            lid_velocity,
            nu,
            Direction<D2Q9<float>::d>::TOP);
    
    Solver<D2Q9<float>, float> solver(lattice);

<<<<<<< Updated upstream
    unsigned long n_iter = 5000;
=======
    unsigned long n_iter = 10000;
>>>>>>> Stashed changes
    float delta_t = 1.0f;
    solver.solve(n_iter, delta_t);

    return 0;
}
