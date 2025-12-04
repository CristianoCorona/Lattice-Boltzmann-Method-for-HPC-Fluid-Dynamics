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

// assumption: Lattice has a function to which we pass coordinates, direction, rho_w, u_w and it gives us a boolean (to indicate if it's a 
// boundary or not) and, if true, also rho_w, u_w of the wall/lid through parameters
// so boundary_values is a bool: bool boundary_values(const std::vector<int> &coords, int i, float_type &rho_b, float_type &u_b);

#include "lattice.hpp"

#ifndef SOLVER_HPP
#define SOLVER_HPP

template <isDescriptor Descriptor, std::floating_point float_type = double>
class Solver {

private:

    constexpr float_type tau_coefficient = 0.9; // to compute tau from delta_t
    Lattice<Descriptor, float_type> &lattice;

    // helper function to compute equilibrium
    inline double feq(int ci, float_type wi, float_type rho, std::array<float_type, Descriptor::d> &u); 
    inline double f_star(float_type fi, float_type feq_i, float_type tau, float_type delta_t);

public:

    // constructor
    Solver(double kinematic_viscosity, Lattice<Descriptor, float_type> &lattice) : lattice(lattice); 

    void solve(int n_iterations, float_type delta_t);

    void compute_moments(Lattice<Descriptor>& grid);
    void stream_collide_inner(Lattice<Descriptor>& grid);
    
    void apply_bcs(Lattice<Descriptor>& grid, double u_lid);

    // function for the output, maybe change second parameter, 
    // use flag to set what output you want
    void write_output(Lattice<Descriptor>& grid, char* out_buffer, char flag);

};

#endif
