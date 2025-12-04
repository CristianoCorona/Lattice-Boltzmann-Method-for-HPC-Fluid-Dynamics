// Solver for LBM simulation
// Handles:
// - Stream and collide operations (merged to avoid race conditions)
// - Boundary conditions
// - Moment computation
// - Output writing
//
// Notes:
// - feq (equilibrium) computed as register variable
// - fi* (collision and propagation are merged)
// - Each cell pulls from its neighbors (no push) to avoid race conditions

#include "lattice.hpp"

#ifndef SOLVER_HPP
#define SOLVER_HPP

template <isDescriptor Descriptor>
class Solver {

private:

    double tau;
    Lattice<Descriptor> &lattice;

    // helper function to compute equilibrium
    inline double feq(int i, double rho, double ux, double uy); 

public:

    // constructor
    Solver(double kinematic_viscosity, Lattice<Descriptor> &lattice) : lattice(lattice); 

    void compute_moments(Lattice<Descriptor>& grid);
    void stream_collide_inner(Lattice<Descriptor>& grid);
    
    void apply_bcs(Lattice<Descriptor>& grid, double u_lid);

    // function for the output, maybe change second parameter, 
    // use flag to set what output you want
    void write_output(Lattice<Descriptor>& grid, char* out_buffer, char flag);

};

#endif
