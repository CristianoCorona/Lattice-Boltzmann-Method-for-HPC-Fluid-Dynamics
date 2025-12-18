#include <iostream>
#include <array>
#include <LBM/LBM.hpp>
#include <omp.h>

using namespace std;

int main(int argc, char *argv[]) {


    array<int, D3Q19<float>::d> lattice_dim = {19, 19, 19};
    array<float, D3Q19<float>::d> lid_velocity = {0.1f, 0.0f, 0.0f};
    float nu = 1.0f;

    Lattice<D3Q19<float>, float> lattice(
            lattice_dim,
            lid_velocity,
            nu,
            Direction<D3Q19<float>::d>::TOP);
    
    Solver<D3Q19<float>, float> solver(lattice);

    unsigned long n_iter = 1000;
    float delta_t = 1.0f;

    cout << "Lattice Descriptor: D" << lattice.d << "Q" << lattice.q << endl;

    cout << "Dimension of the Lattice: " ;
    for(int i=0; i<lattice.sizes.size(); ++i){
        cout << lattice.sizes[i];
        if (i<lattice.d -1) cout << "x";
    }
    cout << endl;

    cout << "Lid velocity on x: " << lattice.u_lid[0] << endl;
    cout << "Lid velocity on y: " << lattice.u_lid[1] << endl;
    if(lattice.d == 3)
        cout << "Lid velocity on z: " << lattice.u_lid[2] << endl;
    

    double initial_time;
    double finish_time;
    double solve_time;

    initial_time = omp_get_wtime();
    solver.solve(n_iter, delta_t);
    finish_time = omp_get_wtime();

    solve_time = finish_time - initial_time;

    cout << "Total Time Elapsed: " << solve_time << "sec" << endl;

    return 0;
}
